
#include "gdalcubes/src/gdalcubes.h"

// [[Rcpp::plugins("cpp11")]]
// [[Rcpp::depends(RcppProgress)]]
#include <Rcpp.h>
#include <progress.hpp>
#include <progress_bar.hpp>
#include <memory>
#include <thread>
#include <algorithm>


using namespace Rcpp;
using namespace gdalcubes;



/**
 * @brief Implementation of the chunk_processor class for multithreaded parallel chunk processing, interruptible by R
 */
class chunk_processor_multithread_interruptible : public chunk_processor {
public:
  /**
   * @brief Construct a multithreaded chunk processor
   * @param nthreads number of threads
   */
  chunk_processor_multithread_interruptible(uint16_t nthreads) : _nthreads(nthreads) {}
  
  /**
   * @copydoc chunk_processor::max_threads
   */
  uint32_t max_threads() override {
    return _nthreads;
  }
  
  /**
   * @copydoc chunk_processor::apply
   */
  void apply(std::shared_ptr<cube> c,
             std::function<void(chunkid_t, std::shared_ptr<chunk_data>, std::mutex &)> f) override;
  
  /**
   * Query the number of threads to be used in parallel chunk processing
   * @return the number of threads
   */
  inline uint16_t get_threads() { return _nthreads; }
  
private:
  uint16_t _nthreads;
};

void chunk_processor_multithread_interruptible::apply(std::shared_ptr<cube> c,
                                        std::function<void(chunkid_t, std::shared_ptr<chunk_data>, std::mutex &)> f) {
  std::mutex mutex;
  bool interrupted = false;
  std::vector<std::thread> workers;
  std::vector<bool> finished(_nthreads, false);
  for (uint16_t it = 0; it < _nthreads; ++it) {
    workers.push_back(std::thread([this, &c, f, it, &mutex, &finished, &interrupted](void) {
      for (uint32_t i = it; i < c->count_chunks(); i += _nthreads) {
        try {
          if (!interrupted) {
            std::shared_ptr<chunk_data> dat = c->read_chunk(i);
            f(i, dat, mutex);
          }
        } catch (std::string s) {
          GCBS_ERROR(s);
          continue;
        } catch (...) {
          GCBS_ERROR("unexpected exception while processing chunk " + std::to_string(i));
          continue;
        }
      }
      finished[it] = true;
    }));
  }
  
  bool done = false;
  int i=0;
  while (!done) {
    done = true;
    for (uint16_t it = 0; it < _nthreads; ++it) {
      done = done && finished[it];
    }
    if (!done) {
      if (Progress::check_abort()) {
        interrupted = true; // still need to wait for threads to finish current chunk
      }
      uint32_t cur_ms = std::min(750, 100 + i*50);
      std::this_thread::sleep_for(std::chrono::milliseconds(cur_ms));
    }
    ++i;
  }
  for (uint16_t it = 0; it < _nthreads; ++it) {
    workers[it].join();
  }
  if (interrupted) {
    throw std::string("computations have been interrupted by the user");
  }

}









struct error_handling_r {
  static std::mutex _m_errhandl;
  static std::stringstream _err_stream;
  static bool _defer;
  
  static void defer_output() {
    _m_errhandl.lock();
    _defer = true;
    _m_errhandl.unlock();
  }
  
  static void do_output() {
    _m_errhandl.lock();
    _defer = false;
    Rcpp::Rcerr << _err_stream.str() << std::endl;
    _err_stream.str(""); 
    _m_errhandl.unlock();
  }
  
  static void debug(error_level type, std::string msg, std::string where, int error_code) {
    _m_errhandl.lock();
    std::string code = (error_code != 0) ? " (" + std::to_string(error_code) + ")" : "";
    std::string where_str = (where.empty()) ? "" : " [in " + where + "]";
    if (type == error_level::ERRLVL_ERROR || type == error_level::ERRLVL_FATAL ) {
      _err_stream << "Error  message: "  << msg << where_str << std::endl;
    } else if (type == error_level::ERRLVL_WARNING) {
      _err_stream << "Warning  message: " << msg << where_str << std::endl;
    } else if (type == error_level::ERRLVL_INFO) {
      _err_stream << "Info message: " << msg << where_str << std::endl;
    } else if (type == error_level::ERRLVL_DEBUG) {
      _err_stream << "Debug message: "  << msg << where_str << std::endl;
    }
    if (!_defer) {
      Rcpp::Rcerr << _err_stream.str() ;
      _err_stream.str(""); 
    }
    _m_errhandl.unlock();
  }
  
  static void standard(error_level type, std::string msg, std::string where, int error_code) {
    _m_errhandl.lock();
    std::string code = (error_code != 0) ? " (" + std::to_string(error_code) + ")" : "";
    if (type == error_level::ERRLVL_ERROR || type == error_level::ERRLVL_FATAL) {
      _err_stream << "Error: " << msg << std::endl;
    } else if (type == error_level::ERRLVL_WARNING) {
      _err_stream << "Warning: " << msg << std::endl;
    } else if (type == error_level::ERRLVL_INFO) {
      _err_stream << "## " << msg << std::endl;
    }
    if (!_defer) {
      Rcpp::Rcerr << _err_stream.str();
      _err_stream.str(""); 
    }
    
    _m_errhandl.unlock();
  }
};
std::mutex error_handling_r::_m_errhandl;
std::stringstream error_handling_r::_err_stream;
bool error_handling_r::_defer = false;



struct progress_simple_R : public progress {
  std::shared_ptr<progress> get() override { return std::make_shared<progress_simple_R>(); }
  
  void set(double p) override {
    _m.lock();
    _set(p);
    _m.unlock();
  };
  
  void increment(double dp) override {
    _m.lock();
    _set(_p + dp);
    _m.unlock();
  }
  virtual void finalize() override {
    _m.lock();
    _rp->update(100);
    error_handling_r::do_output();
    _m.unlock();
  }


  progress_simple_R() : _p(0), _rp(nullptr) {}

  ~progress_simple_R(){
    if (_rp) {
      delete _rp;
    }
  }
  
 
  
  

private:
  
  std::mutex _m;
  double _p;
  Progress *_rp;
  
  void _set(double p) { // call this function only with a lock on _m
    //Rcpp::checkUserInterrupt();
    // if (Progress::check_abort()) {
    //   throw std::string("Operation has been interrupted by user");
    // }
    if (!_rp) {
      error_handling_r::defer_output();
      _rp = new Progress(100,true);
    }
    _p = p;
    _rp->update((int)(_p*100));
  }
};

// see https://stackoverflow.com/questions/26666614/how-do-i-check-if-an-externalptr-is-null-from-within-r
// [[Rcpp::export]]
Rcpp::LogicalVector libgdalcubes_is_null(SEXP pointer) {
  return Rcpp::LogicalVector(!R_ExternalPtrAddr(pointer));
}


// [[Rcpp::export]]
Rcpp::List libgdalcubes_version() {
  version_info v = config::instance()->get_version_info();
  return(Rcpp::List::create(
      Rcpp::Named("VERSION_MAJOR") = v.VERSION_MAJOR ,
      Rcpp::Named("VERSION_MINOR") = v.VERSION_MINOR ,
      Rcpp::Named("VERSION_PATCH") = v.VERSION_PATCH ,
      Rcpp::Named("BUILD_DATE") = v.BUILD_DATE,
      Rcpp::Named("BUILD_TIME") = v.BUILD_TIME,
      Rcpp::Named("GIT_DESC") = v.GIT_DESC,
      Rcpp::Named("GIT_COMMIT") = v.GIT_COMMIT));
}

// [[Rcpp::export]]
std::vector<std::string> libgdalcubes_gdalformats() {
  return config::instance()->gdal_formats();
}

// [[Rcpp::export]]
std::string libgdalcubes_gdalversion() {
  return config::instance()->gdal_version_info();
}

// [[Rcpp::export]]
bool libgdalcubes_gdal_has_geos() {
  return config::instance()->gdal_has_geos();
}


// [[Rcpp::export]]
void libgdalcubes_add_format_dir(std::string dir) {
  config::instance()->add_collection_format_preset_dir(dir);
}

// [[Rcpp::export]]
void libgdalcubes_init() {
  config::instance()->gdalcubes_init();
  config::instance()->set_default_progress_bar(std::make_shared<progress_simple_R>());
  //config::instance()->set_default_progress_bar(std::make_shared<progress_none>());
  
  config::instance()->set_error_handler(error_handling_r::standard); 
  
  // Interruptible chunk processor
  config::instance()->set_default_chunk_processor(std::dynamic_pointer_cast<chunk_processor>(std::make_shared<chunk_processor_multithread_interruptible>(1)));
  
}

// [[Rcpp::export]]
void libgdalcubes_cleanup() {
  config::instance()->gdalcubes_cleanup();
}

// This function is also covered by libgdalcubes_dimension_values() and should be deprecated
// [[Rcpp::export]]
Rcpp::StringVector libgdalcubes_datetime_values(SEXP pin) {
  Rcpp::XPtr<std::shared_ptr<cube>> aa = Rcpp::as<Rcpp::XPtr<std::shared_ptr<cube>>>(pin);
  std::shared_ptr<cube> x = *aa;
  
  Rcpp::CharacterVector out(x->size_t());
  
  for (uint32_t i = 0; i < x->size_t(); ++i) {
    out[i] = x->st_reference()->datetime_at_index(i).to_string();
  }
  return out;
}

// [[Rcpp::export]]
Rcpp::List libgdalcubes_cube_info( SEXP pin) {

  try {
    Rcpp::XPtr<std::shared_ptr<cube>> aa = Rcpp::as<Rcpp::XPtr<std::shared_ptr<cube>>>(pin);
    
    std::shared_ptr<cube> x = *aa;
    
    Rcpp::List tdim_list;
    if (x->st_reference()->has_regular_time()) {
      std::shared_ptr<cube_stref_regular> stref = std::dynamic_pointer_cast<cube_stref_regular>( x->st_reference());
      tdim_list = Rcpp::List::create(
        Rcpp::Named("low") = stref->t0().to_string(),
        Rcpp::Named("high") = stref->t1().to_string(),
        Rcpp::Named("count") = stref->nt(),
        Rcpp::Named("pixel_size") = stref->dt().to_string(),
        Rcpp::Named("chunk_size") = x->chunk_size()[0]);
    }
    else {
      std::shared_ptr<cube_stref_labeled_time> stref = std::dynamic_pointer_cast<cube_stref_labeled_time>( x->st_reference());
      tdim_list = Rcpp::List::create(
        Rcpp::Named("low") = stref->t0().to_string(),
        Rcpp::Named("high") = stref->t1().to_string(),
        Rcpp::Named("count") = stref->nt(),
        Rcpp::Named("values") = stref->get_time_labels_as_string(), 
        Rcpp::Named("pixel_size") = stref->dt().to_string(),
        Rcpp::Named("chunk_size") = x->chunk_size()[0]);
    }
    
    Rcpp::List dims =
      Rcpp::List::create(
        Rcpp::Named("t")= tdim_list,
        Rcpp::Named("y") =  
          Rcpp::List::create(
            Rcpp::Named("low") = x->st_reference()->bottom(),
            Rcpp::Named("high") = x->st_reference()->top(),
            Rcpp::Named("count") = x->st_reference()->ny(),
            Rcpp::Named("pixel_size") = x->st_reference()->dy(),
            Rcpp::Named("chunk_size") = x->chunk_size()[1]),
        Rcpp::Named("x")=
          Rcpp::List::create(
            Rcpp::Named("low") = x->st_reference()->left(),
            Rcpp::Named("high") = x->st_reference()->right(),
            Rcpp::Named("count") = x->st_reference()->nx(),
            Rcpp::Named("pixel_size") = x->st_reference()->dx(),
            Rcpp::Named("chunk_size") = x->chunk_size()[2])
      );
                         
                         
                      
                       
    
    
    
    Rcpp::CharacterVector b_names(x->bands().count(), "");
    //Rcpp::CharacterVector b_type(x->bands().count(), "");
    Rcpp::NumericVector b_offset(x->bands().count(), NA_REAL);
    Rcpp::NumericVector b_scale(x->bands().count(), NA_REAL);
    Rcpp::NumericVector b_nodata(x->bands().count(), NA_REAL);
    Rcpp::CharacterVector b_unit(x->bands().count(), "");
    
    for (uint16_t i=0; i<x->bands().count(); ++i) {
      b_names[i] = x->bands().get(i).name;
     // b_type[i] = x->bands().get(i).type;
      b_offset[i] = x->bands().get(i).offset;
      b_scale[i] = x->bands().get(i).scale;
      b_nodata[i] = (x->bands().get(i).no_data_value.empty())? NA_REAL : std::stod(x->bands().get(i).no_data_value);
      b_unit[i] = x->bands().get(i).unit;
    }
    
    Rcpp::DataFrame bands =
      Rcpp::DataFrame::create(Rcpp::Named("name")=b_names,
                              //Rcpp::Named("type")=b_type,
                              Rcpp::Named("offset")=b_offset,
                              Rcpp::Named("scale")=b_scale,
                              Rcpp::Named("nodata")=b_nodata,
                              Rcpp::Named("unit")=b_unit);
    
    // extract proj4 string 
    char* proj4 = NULL;
    std::string sproj4 = "";
    if (x->st_reference()->srs_ogr().exportToProj4(&proj4) == OGRERR_NONE) {
      sproj4 = proj4;  
      CPLFree(proj4);
    }
    
    return Rcpp::List::create(Rcpp::Named("bands") = bands,
                              Rcpp::Named("dimensions") = dims,
                              Rcpp::Named("srs") = x->st_reference()->srs(),
                              Rcpp::Named("proj4") = sproj4,
                              Rcpp::Named("graph") = x->make_constructible_json().dump(),
                              Rcpp::Named("size") = Rcpp::IntegerVector::create(x->size()[0], x->size()[1], x->size()[2], x->size()[3])); // TODO: remove size element
    
  }
  catch (std::string s) {
    Rcpp::stop(s);
  }
  
  
  
}

// [[Rcpp::export]]
Rcpp::List libgdalcubes_dimension_values_from_view(Rcpp::List view, std::string dt_unit="") {
  
  cube_view cv;
  
  if (Rcpp::as<Rcpp::List>(view["space"])["right"] != R_NilValue) {
    cv.right(Rcpp::as<Rcpp::List>(view["space"])["right"]);
  }
  if (Rcpp::as<Rcpp::List>(view["space"])["left"] != R_NilValue) {
    cv.left(Rcpp::as<Rcpp::List>(view["space"])["left"]);
  }
  if (Rcpp::as<Rcpp::List>(view["space"])["top"] != R_NilValue) {
    cv.top(Rcpp::as<Rcpp::List>(view["space"])["top"]);
  }
  if (Rcpp::as<Rcpp::List>(view["space"])["bottom"] != R_NilValue) {
    cv.bottom(Rcpp::as<Rcpp::List>(view["space"])["bottom"]);
  }
  if (Rcpp::as<Rcpp::List>(view["space"])["dx"] != R_NilValue) {
    cv.dx(Rcpp::as<Rcpp::List>(view["space"])["dx"]);
  }
  if (Rcpp::as<Rcpp::List>(view["space"])["nx"] != R_NilValue) {
    cv.nx(Rcpp::as<Rcpp::List>(view["space"])["nx"]);
  }
  if (Rcpp::as<Rcpp::List>(view["space"])["dy"] != R_NilValue) {
    cv.dy(Rcpp::as<Rcpp::List>(view["space"])["dy"]);
  }
  if (Rcpp::as<Rcpp::List>(view["space"])["ny"] != R_NilValue) {
    cv.ny(Rcpp::as<Rcpp::List>(view["space"])["ny"]);
  }
  if (Rcpp::as<Rcpp::List>(view["space"])["srs"] != R_NilValue) {
      std::string srs = Rcpp::as<Rcpp::String>(Rcpp::as<Rcpp::List>(view["space"])["srs"]);
      cv.srs(srs);
  }
  if (Rcpp::as<Rcpp::List>(view["time"])["t0"] != R_NilValue) {
    std::string tmp = Rcpp::as<Rcpp::String>(Rcpp::as<Rcpp::List>(view["time"])["t0"]);
    cv.t0(datetime::from_string(tmp));
  }
  if (Rcpp::as<Rcpp::List>(view["time"])["t1"] != R_NilValue) {
    std::string tmp = Rcpp::as<Rcpp::String>(Rcpp::as<Rcpp::List>(view["time"])["t1"]);
    cv.t1(datetime::from_string(tmp));
  }
  if (Rcpp::as<Rcpp::List>(view["time"])["nt"] != R_NilValue) {
    cv.nt(Rcpp::as<Rcpp::List>(view["time"])["nt"]);
  }
  if (Rcpp::as<Rcpp::List>(view["time"])["dt"] != R_NilValue) {
    std::string tmp = Rcpp::as<Rcpp::String>(Rcpp::as<Rcpp::List>(view["time"])["dt"]);
    cv.dt(duration::from_string(tmp));
    cv.t0().unit() = cv.dt().dt_unit; 
    cv.t1().unit() = cv.dt().dt_unit; 
  }
  
  if (view["aggregation"] != R_NilValue) {
    std::string tmp = Rcpp::as<Rcpp::String>(view["aggregation"]);
    cv.aggregation_method() = aggregation::from_string(tmp);
  }
  if (view["resampling"] != R_NilValue) {
    std::string tmp = Rcpp::as<Rcpp::String>(view["resampling"]);
    cv.resampling_method() = resampling::from_string(tmp);
  }
  
  Rcpp::CharacterVector dimt(cv.nt());
  Rcpp::NumericVector dimx(cv.nx());
  Rcpp::NumericVector dimy(cv.ny());
  
  
  datetime_unit u = cv.dt_unit();
  if (dt_unit == "Y") {
    u = datetime_unit::YEAR;
  }
  else if (dt_unit == "m") {
    u = datetime_unit::MONTH;
  }
  else if (dt_unit == "d") {
    u = datetime_unit::DAY;
  }
  else if (dt_unit == "H") {
    u = datetime_unit::HOUR;
  }
  else if (dt_unit == "M") {
    u = datetime_unit::MINUTE;
  }
  else if (dt_unit == "S") {
    u = datetime_unit::SECOND;
  }
  
  for (uint32_t i = 0; i < cv.nt(); ++i) {
    dimt[i] = (cv.t0() + cv.dt() * i).to_string(u); 
  }
  for (uint32_t i = 0; i < cv.ny(); ++i) {
    dimy[i] = (cv.win().bottom +cv.dy() * i); 
  }
  for (uint32_t i = 0; i < cv.nx(); ++i) {
    dimx[i] = (cv.win().left + cv.dx() * i); 
  }
  
  return Rcpp::List::create(Rcpp::Named("t") = dimt,
                            Rcpp::Named("y") = dimy,
                            Rcpp::Named("x") = dimx);
  
  
}




// [[Rcpp::export]]
Rcpp::List libgdalcubes_dimension_bounds(SEXP pin, std::string dt_unit="") {
  
  Rcpp::XPtr<std::shared_ptr<cube>> aa = Rcpp::as<Rcpp::XPtr<std::shared_ptr<cube>>>(pin);
  
  std::shared_ptr<cube> x = *aa;
  Rcpp::CharacterVector dimt(2*x->st_reference()->nt());
  Rcpp::NumericVector dimx(2*x->st_reference()->nx());
  Rcpp::NumericVector dimy(2*x->st_reference()->ny());
  
  if (!x->st_reference()->has_regular_space()) {
    Rcpp::stop("Irregular spatial dimensions are currently not supprted");
  }
  
  // NOTE: the following will only work as long as all cube st reference types with regular spatial dimensions inherit from  cube_stref_regular class
  std::shared_ptr<cube_stref_regular> stref = std::dynamic_pointer_cast<cube_stref_regular>(x->st_reference());
  
  
  datetime_unit u = stref->dt_unit();
  if (dt_unit == "Y") {
    u = datetime_unit::YEAR;
  }
  else if (dt_unit == "m") {
    u = datetime_unit::MONTH;
  }
  else if (dt_unit == "d") {
    u = datetime_unit::DAY;
  }
  else if (dt_unit == "H") {
    u = datetime_unit::HOUR;
  }
  else if (dt_unit == "M") {  
    u = datetime_unit::MINUTE;
  }
  else if (dt_unit == "S") {
    u = datetime_unit::SECOND;
  }
  
  for (uint32_t i = 0; i < x->st_reference()->nt(); ++i) {
    dimt[2*i] = stref->datetime_at_index(i).to_string(u);
    dimt[2*i+1] = stref->datetime_at_index(i+1).to_string(u);
  }
  for (uint32_t i = 0; i < x->st_reference()->ny(); ++i) {
    dimy[2*i] = (stref->win().bottom + stref->dy() * i); 
    dimy[2*i+1] = (stref->win().bottom + stref->dy() * (i+1)); 
  }
  for (uint32_t i = 0; i < x->st_reference()->nx(); ++i) {
    dimx[2*i] = (stref->win().left + stref->dx() * i); 
    dimx[2*i+1] = (stref->win().left + stref->dx() * (i+1)); 
  }
  
  return Rcpp::List::create(Rcpp::Named("t") = dimt,
                            Rcpp::Named("y") = dimy,
                            Rcpp::Named("x") = dimx);
  
}

// [[Rcpp::export]]
Rcpp::List libgdalcubes_dimension_values(SEXP pin, std::string dt_unit="") {
  
  Rcpp::XPtr<std::shared_ptr<cube>> aa = Rcpp::as<Rcpp::XPtr<std::shared_ptr<cube>>>(pin);
  
  std::shared_ptr<cube> x = *aa;
  Rcpp::CharacterVector dimt(x->st_reference()->nt());
  Rcpp::NumericVector dimx(x->st_reference()->nx());
  Rcpp::NumericVector dimy(x->st_reference()->ny());
  
  if (!x->st_reference()->has_regular_space()) {
    Rcpp::stop("Irregular spatial dimensions are currently not supprted");
  }
  
  // NOTE: the following will only work as long as all cube st reference types with regular spatial dimensions inherit from  cube_stref_regular class
  std::shared_ptr<cube_stref_regular> stref = std::dynamic_pointer_cast<cube_stref_regular>(x->st_reference());
  
  
  datetime_unit u = stref->dt_unit();
  if (dt_unit == "Y") {
    u = datetime_unit::YEAR;
  }
  else if (dt_unit == "m") {
    u = datetime_unit::MONTH;
  }
  else if (dt_unit == "d") {
    u = datetime_unit::DAY;
  }
  else if (dt_unit == "H") {
    u = datetime_unit::HOUR;
  }
  else if (dt_unit == "M") {  
    u = datetime_unit::MINUTE;
  }
  else if (dt_unit == "S") {
    u = datetime_unit::SECOND;
  }
  
  for (uint32_t i = 0; i < x->st_reference()->nt(); ++i) {
    dimt[i] = stref->datetime_at_index(i).to_string(u);
  }
  for (uint32_t i = 0; i < x->st_reference()->ny(); ++i) {
    dimy[i] = (stref->win().bottom + stref->dy() * i); 
  }
  for (uint32_t i = 0; i < x->st_reference()->nx(); ++i) {
    dimx[i] = (stref->win().left + stref->dx() * i); 
  }
  
  return Rcpp::List::create(Rcpp::Named("t") = dimt,
                            Rcpp::Named("y") = dimy,
                            Rcpp::Named("x") = dimx);
  
}



// [[Rcpp::export]]
SEXP libgdalcubes_open_image_collection(std::string filename) {
  
  try {
    std::shared_ptr<image_collection>* x = new std::shared_ptr<image_collection>( std::make_shared<image_collection>(filename));
    Rcpp::XPtr< std::shared_ptr<image_collection> > p(x, true) ;
    return p;
  }
  catch (std::string s) {
    Rcpp::stop(s);
  }
}



// [[Rcpp::export]]
Rcpp::List libgdalcubes_image_collection_info( SEXP pin) {
  
  try {
    Rcpp::XPtr<std::shared_ptr<image_collection>> aa = Rcpp::as<Rcpp::XPtr<std::shared_ptr<image_collection>>>(pin);
    std::shared_ptr<image_collection> ic = *aa;
    
    std::vector<image_collection::images_row> img = ic->get_images();
    
    Rcpp::CharacterVector images_name(img.size());
    Rcpp::NumericVector images_left(img.size());
    Rcpp::NumericVector images_top(img.size());
    Rcpp::NumericVector images_right(img.size());
    Rcpp::NumericVector images_bottom(img.size());
    Rcpp::CharacterVector images_datetime(img.size());
    Rcpp::CharacterVector images_proj(img.size());
    
    
    for (uint32_t i=0; i<img.size(); ++i) {
      images_name[i] = img[i].name;
      images_left[i] = img[i].left;
      images_right[i] = img[i].right;
      images_top[i] = img[i].top;
      images_bottom[i] = img[i].bottom;
      images_proj[i] = img[i].proj;
      images_datetime[i] = img[i].datetime;
    }
    
    
    Rcpp::DataFrame images_df =
      Rcpp::DataFrame::create(Rcpp::Named("name")=images_name,
                              Rcpp::Named("left")=images_left,
                              Rcpp::Named("top")=images_top,
                              Rcpp::Named("bottom")=images_bottom,
                              Rcpp::Named("right")=images_right,
                              Rcpp::Named("datetime")=images_datetime,
                              Rcpp::Named("srs")=images_proj);
    
    
    
    
    std::vector<image_collection::bands_row> bands = ic->get_available_bands();
    
    
    Rcpp::CharacterVector bands_name( bands.size());
    //Rcpp::CharacterVector bands_type(bands.size());
    Rcpp::NumericVector bands_offset( bands.size());
    Rcpp::NumericVector bands_scale( bands.size());
    Rcpp::CharacterVector bands_unit( bands.size());
    Rcpp::CharacterVector bands_nodata( bands.size());
    Rcpp::IntegerVector bands_image_count( bands.size());
    
    for (uint32_t i=0; i<bands.size(); ++i) {
      bands_name[i] = bands[i].name;
      //bands_type[i] = bands[i].type;
      bands_offset[i] = bands[i].offset;
      bands_scale[i] = bands[i].scale;
      bands_unit[i] = bands[i].unit;
      bands_nodata[i] = bands[i].nodata;
      bands_image_count[i] = bands[i].image_count;
    }
    
    Rcpp::DataFrame bands_df =
      Rcpp::DataFrame::create(Rcpp::Named("name")=bands_name,
                              //Rcpp::Named("type")=bands_type,
                              Rcpp::Named("offset")=bands_offset,
                              Rcpp::Named("scale")=bands_scale,
                              Rcpp::Named("unit")=bands_unit,
                              Rcpp::Named("nodata")=bands_nodata,
                              Rcpp::Named("image_count")=bands_image_count);
    
    
    std::vector<image_collection::gdalrefs_row> gdalrefs = ic->get_gdalrefs();
    
    Rcpp::IntegerVector gdalrefs_imageid(gdalrefs.size());
    Rcpp::IntegerVector gdalrefs_bandid(gdalrefs.size());
    Rcpp::CharacterVector gdalrefs_descriptor(gdalrefs.size());
    Rcpp::IntegerVector gdalrefs_bandnum(gdalrefs.size());
    
    for (uint32_t i=0; i<gdalrefs.size(); ++i) {
      gdalrefs_imageid[i] = gdalrefs[i].image_id;
      gdalrefs_bandid[i] = gdalrefs[i].band_id;
      gdalrefs_descriptor[i] = gdalrefs[i].descriptor;
      gdalrefs_bandnum[i] = gdalrefs[i].band_num;
    }
    
    Rcpp::DataFrame gdalrefs_df =
      Rcpp::DataFrame::create(Rcpp::Named("image_id")=gdalrefs_imageid,
                              Rcpp::Named("band_id")=gdalrefs_bandid,
                              Rcpp::Named("descriptor")=gdalrefs_descriptor,
                              Rcpp::Named("band_num")=gdalrefs_bandnum);
    
    
    return Rcpp::List::create(Rcpp::Named("images") = images_df,
                              Rcpp::Named("bands") = bands_df,
                              Rcpp::Named("gdalrefs") = gdalrefs_df);
    
  }
  catch (std::string s) {
    Rcpp::stop(s);
  }
}





// [[Rcpp::export]]
Rcpp::List libgdalcubes_image_collection_extent( SEXP pin, std::string srs) {
  
  try {
    Rcpp::XPtr<std::shared_ptr<image_collection>> aa = Rcpp::as<Rcpp::XPtr<std::shared_ptr<image_collection>>>(pin);
    std::shared_ptr<image_collection> ic = *aa;
    
    bounds_st ext = ic->extent();
    ext.s = ext.s.transform("EPSG:4326", srs.c_str());
    
    return Rcpp::List::create(Rcpp::Named("left") = ext.s.left,
                              Rcpp::Named("right") = ext.s.right,
                              Rcpp::Named("top") = ext.s.top,
                              Rcpp::Named("bottom") = ext.s.bottom,
                              Rcpp::Named("t0") = ext.t0.to_string(),
                              Rcpp::Named("t1") = ext.t1.to_string());
    
  }
  catch (std::string s) {
    Rcpp::stop(s);
  }
}







// [[Rcpp::export]]
void libgdalcubes_create_image_collection_from_format(std::vector<std::string> files, std::string format_file, std::string outfile, bool unroll_archives=true) {
  try {
    collection_format cfmt(format_file);
    if (unroll_archives) {
      files = image_collection::unroll_archives(files);
    }
    image_collection::create(cfmt, files)->write(outfile);
  }
  catch (std::string s) {
    Rcpp::stop(s);
  }
}

// [[Rcpp::export]]
void libgdalcubes_create_image_collection_from_datetime(std::string outfile, std::vector<std::string> files, 
                                                        std::vector<std::string> date_time, bool use_subdatasets, 
                                                        std::vector<std::string> band_names) {
  try {
    image_collection::create(files, date_time, band_names, use_subdatasets)->write(outfile);
  }
  catch (std::string s) {
    Rcpp::stop(s);
  }
}


// [[Rcpp::export]]
void libgdalcubes_add_images(SEXP pin, std::vector<std::string> files, bool unroll_archives=true, std::string outfile = "") {
  
  try {
    Rcpp::XPtr<std::shared_ptr<image_collection>> aa = Rcpp::as<Rcpp::XPtr<std::shared_ptr<image_collection>>>(pin);
    if (!outfile.empty()) {
      (*aa)->write(outfile);
    }
    if (unroll_archives) {
      files = image_collection::unroll_archives(files);
    }
    (*aa)->add_with_collection_format(files);
  }
  catch (std::string s) {
    Rcpp::stop(s);
  }
}

// [[Rcpp::export]]
SEXP libgdalcubes_list_collection_formats() {
  try {
    
    // get package directory and add to presets... (file.path(system.file(package="gdalcubes"),"formats"))
    Rcpp::Environment base("package:base"); 
    Rcpp::Function sfile = base["system.file"];  
    Rcpp::Function fpath = base["file.path"];
    Rcpp::CharacterVector preset_dir = fpath(sfile(Rcpp::_["package"] = "gdalcubes"), "formats");
    std::string temp = Rcpp::as<std::string>(preset_dir[0]);
    config::instance()->add_collection_format_preset_dir(temp);
    
    
    std::map<std::string,std::string> fmts = collection_format::list_presets();
    Rcpp::CharacterVector out_keys(fmts.size());
    Rcpp::CharacterVector out_values(fmts.size());
    uint32_t i=0;
    for (auto it=fmts.begin(); it != fmts.end(); ++it) {
      out_values[i] = it->second;
      out_keys[i++] = it->first;
      
    }
    return  Rcpp::DataFrame::create(Rcpp::Named("name")=out_keys,
                                    Rcpp::Named("path")=out_values);
  }
  catch (std::string s) {
    Rcpp::stop(s);
  }
}

// [[Rcpp::export]]
SEXP libgdalcubes_create_view(SEXP v) {
  Rcpp::List view = Rcpp::as<Rcpp::List>(v);
  cube_view cv;
  
  if (Rcpp::as<Rcpp::List>(view["space"])["right"] != R_NilValue) {
    cv.right(Rcpp::as<Rcpp::List>(view["space"])["right"]);
  }
  if (Rcpp::as<Rcpp::List>(view["space"])["left"] != R_NilValue) {
    cv.left(Rcpp::as<Rcpp::List>(view["space"])["left"]);
  }
  if (Rcpp::as<Rcpp::List>(view["space"])["top"] != R_NilValue) {
    cv.top(Rcpp::as<Rcpp::List>(view["space"])["top"]);
  }
  if (Rcpp::as<Rcpp::List>(view["space"])["bottom"] != R_NilValue) {
    cv.bottom(Rcpp::as<Rcpp::List>(view["space"])["bottom"]);
  }
  if (Rcpp::as<Rcpp::List>(view["space"])["dx"] != R_NilValue) {
    cv.dx(Rcpp::as<Rcpp::List>(view["space"])["dx"]);
  }
  if (Rcpp::as<Rcpp::List>(view["space"])["nx"] != R_NilValue) {
    cv.nx(Rcpp::as<Rcpp::List>(view["space"])["nx"]);
  }
  if (Rcpp::as<Rcpp::List>(view["space"])["dy"] != R_NilValue) {
    cv.dy(Rcpp::as<Rcpp::List>(view["space"])["dy"]);
  }
  if (Rcpp::as<Rcpp::List>(view["space"])["ny"] != R_NilValue) {
    cv.ny(Rcpp::as<Rcpp::List>(view["space"])["ny"]);
  }
  if (Rcpp::as<Rcpp::List>(view["space"])["srs"] != R_NilValue) {
    std::string srs = Rcpp::as<Rcpp::String>(Rcpp::as<Rcpp::List>(view["space"])["srs"]);
    cv.srs(srs);
  }
  if (Rcpp::as<Rcpp::List>(view["time"])["t0"] != R_NilValue) {
    std::string tmp = Rcpp::as<Rcpp::String>(Rcpp::as<Rcpp::List>(view["time"])["t0"]);
    cv.t0(datetime::from_string(tmp));
  }
  if (Rcpp::as<Rcpp::List>(view["time"])["t1"] != R_NilValue) {
    std::string tmp = Rcpp::as<Rcpp::String>(Rcpp::as<Rcpp::List>(view["time"])["t1"]);
    cv.t1(datetime::from_string(tmp));
  }
  if (Rcpp::as<Rcpp::List>(view["time"])["nt"] != R_NilValue) {
    cv.nt(Rcpp::as<Rcpp::List>(view["time"])["nt"]);
  }
  if (Rcpp::as<Rcpp::List>(view["time"])["dt"] != R_NilValue) {
    std::string tmp = Rcpp::as<Rcpp::String>(Rcpp::as<Rcpp::List>(view["time"])["dt"]);
    cv.dt(duration::from_string(tmp));
    cv.t0().unit() = cv.dt().dt_unit; 
    cv.t1().unit() = cv.dt().dt_unit; 
  }
  
  if (view["aggregation"] != R_NilValue) {
    std::string tmp = Rcpp::as<Rcpp::String>(view["aggregation"]);
    cv.aggregation_method() = aggregation::from_string(tmp);
  }
  if (view["resampling"] != R_NilValue) {
    std::string tmp = Rcpp::as<Rcpp::String>(view["resampling"]);
    cv.resampling_method() = resampling::from_string(tmp);
  }
  
  
  return Rcpp::List::create(
    Rcpp::Named("space") = Rcpp::List::create(
      Rcpp::Named("right") = cv.right(),
      Rcpp::Named("left") = cv.left(),
      Rcpp::Named("top") = cv.top(),
      Rcpp::Named("bottom") = cv.bottom(),
      Rcpp::Named("nx") = cv.nx(),
      Rcpp::Named("ny") = cv.ny(),
      Rcpp::Named("srs") = cv.srs(),
      Rcpp::Named("dx") = cv.dx(),
      Rcpp::Named("dy") = cv.dy()
    ),
    Rcpp::Named("time") = Rcpp::List::create(
      Rcpp::Named("t0") = cv.t0().to_string(),
      Rcpp::Named("t1") = cv.t1().to_string(),
      Rcpp::Named("dt") = cv.dt().to_string(),
      Rcpp::Named("nt") = cv.nt()
    ),
    Rcpp::Named("aggregation") = aggregation::to_string(cv.aggregation_method()) ,
    Rcpp::Named("resampling") = resampling::to_string(cv.resampling_method())
  );
 
}


// [[Rcpp::export]]
SEXP libgdalcubes_create_image_collection_cube(SEXP pin, Rcpp::IntegerVector chunk_sizes, SEXP mask, SEXP v = R_NilValue) {

  try {
    Rcpp::XPtr<std::shared_ptr<image_collection>> aa = Rcpp::as<Rcpp::XPtr<std::shared_ptr<image_collection>>>(pin);
    
    std::shared_ptr<image_collection_cube>* x;
    if (v == R_NilValue) {
      x = new std::shared_ptr<image_collection_cube>( image_collection_cube::create(*aa));
    }
    else {
      Rcpp::List view = Rcpp::as<Rcpp::List>(v);
      cube_view cv;
      
      if (Rcpp::as<Rcpp::List>(view["space"])["right"] != R_NilValue) {
        cv.right(Rcpp::as<Rcpp::List>(view["space"])["right"]);
      }
      if (Rcpp::as<Rcpp::List>(view["space"])["left"] != R_NilValue) {
        cv.left(Rcpp::as<Rcpp::List>(view["space"])["left"]);
      }
      if (Rcpp::as<Rcpp::List>(view["space"])["top"] != R_NilValue) {
        cv.top(Rcpp::as<Rcpp::List>(view["space"])["top"]);
      }
      if (Rcpp::as<Rcpp::List>(view["space"])["bottom"] != R_NilValue) {
        cv.bottom(Rcpp::as<Rcpp::List>(view["space"])["bottom"]);
      }
      if (Rcpp::as<Rcpp::List>(view["space"])["dx"] != R_NilValue) {
        cv.dx(Rcpp::as<Rcpp::List>(view["space"])["dx"]);
      }
      if (Rcpp::as<Rcpp::List>(view["space"])["nx"] != R_NilValue) {
        cv.nx(Rcpp::as<Rcpp::List>(view["space"])["nx"]);
      }
      if (Rcpp::as<Rcpp::List>(view["space"])["dy"] != R_NilValue) {
        cv.dy(Rcpp::as<Rcpp::List>(view["space"])["dy"]);
      }
      if (Rcpp::as<Rcpp::List>(view["space"])["ny"] != R_NilValue) {
        cv.ny(Rcpp::as<Rcpp::List>(view["space"])["ny"]);
      }
      if (Rcpp::as<Rcpp::List>(view["space"])["srs"] != R_NilValue) {
        std::string srs = Rcpp::as<Rcpp::String>(Rcpp::as<Rcpp::List>(view["space"])["srs"]);
        cv.srs(srs);
      }
      if (Rcpp::as<Rcpp::List>(view["time"])["t0"] != R_NilValue) {
        std::string tmp = Rcpp::as<Rcpp::String>(Rcpp::as<Rcpp::List>(view["time"])["t0"]);
        cv.t0(datetime::from_string(tmp));
      }
      if (Rcpp::as<Rcpp::List>(view["time"])["t1"] != R_NilValue) {
        std::string tmp = Rcpp::as<Rcpp::String>(Rcpp::as<Rcpp::List>(view["time"])["t1"]);
        cv.t1(datetime::from_string(tmp));
      }
      if (Rcpp::as<Rcpp::List>(view["time"])["nt"] != R_NilValue) {
        cv.nt(Rcpp::as<Rcpp::List>(view["time"])["nt"]);
      }
      if (Rcpp::as<Rcpp::List>(view["time"])["dt"] != R_NilValue) {
        std::string tmp = Rcpp::as<Rcpp::String>(Rcpp::as<Rcpp::List>(view["time"])["dt"]);
        cv.dt(duration::from_string(tmp));
      }
      
      if (view["aggregation"] != R_NilValue) {
        std::string tmp = Rcpp::as<Rcpp::String>(view["aggregation"]);
        cv.aggregation_method() = aggregation::from_string(tmp);
      }
      if (view["resampling"] != R_NilValue) {
        std::string tmp = Rcpp::as<Rcpp::String>(view["resampling"]);
        cv.resampling_method() = resampling::from_string(tmp);
      }
      x = new std::shared_ptr<image_collection_cube>( image_collection_cube::create(*aa, cv));
    }
    (*x)->set_chunk_size(chunk_sizes[0], chunk_sizes[1], chunk_sizes[2]);
    
    
    if (mask != R_NilValue) {
      std::string band_name = Rcpp::as<Rcpp::List>(mask)["band"]; 
      bool invert = Rcpp::as<Rcpp::List>(mask)["invert"];
      
      if (Rcpp::as<Rcpp::List>(mask)["values"] != R_NilValue) {
        std::vector<double> values = Rcpp::as<std::vector<double>>(Rcpp::as<Rcpp::List>(mask)["values"]);
        std::vector<uint8_t> bits;
        if (Rcpp::as<Rcpp::List>(mask)["bits"] != R_NilValue)
          bits =  Rcpp::as<std::vector<uint8_t>>(Rcpp::as<Rcpp::List>(mask)["bits"]);
        (*x)->set_mask(band_name, std::make_shared<value_mask>(std::unordered_set<double>(values.begin(), values.end()), invert, bits));
      }
      else {
        double min = Rcpp::as<Rcpp::List>(mask)["min"];
        double max = Rcpp::as<Rcpp::List>(mask)["max"];
        std::vector<uint8_t> bits;
        if (Rcpp::as<Rcpp::List>(mask)["bits"] != R_NilValue)
          bits =  Rcpp::as<std::vector<uint8_t>>(Rcpp::as<Rcpp::List>(mask)["bits"]);
        (*x)->set_mask(band_name, std::make_shared<range_mask>(min, max, invert, bits));
      }
    }
    
    
    Rcpp::XPtr< std::shared_ptr<image_collection_cube> > p(x, true) ;
    
    return p;
    
  }
  catch (std::string s) {
    Rcpp::stop(s);
  }
 
  
}











// [[Rcpp::export]]
SEXP libgdalcubes_create_dummy_cube(SEXP v, uint16_t nbands, double fill, Rcpp::IntegerVector chunk_sizes) {
  
  try {
    
      Rcpp::List view = Rcpp::as<Rcpp::List>(v);
      cube_view cv;
      
      if (Rcpp::as<Rcpp::List>(view["space"])["right"] != R_NilValue) {
        cv.right(Rcpp::as<Rcpp::List>(view["space"])["right"]);
      }
      if (Rcpp::as<Rcpp::List>(view["space"])["left"] != R_NilValue) {
        cv.left(Rcpp::as<Rcpp::List>(view["space"])["left"]);
      }
      if (Rcpp::as<Rcpp::List>(view["space"])["top"] != R_NilValue) {
        cv.top(Rcpp::as<Rcpp::List>(view["space"])["top"]);
      }
      if (Rcpp::as<Rcpp::List>(view["space"])["bottom"] != R_NilValue) {
        cv.bottom(Rcpp::as<Rcpp::List>(view["space"])["bottom"]);
      }
      if (Rcpp::as<Rcpp::List>(view["space"])["dx"] != R_NilValue) {
        cv.dx(Rcpp::as<Rcpp::List>(view["space"])["dx"]);
      }
      if (Rcpp::as<Rcpp::List>(view["space"])["nx"] != R_NilValue) {
        cv.nx(Rcpp::as<Rcpp::List>(view["space"])["nx"]);
      }
      if (Rcpp::as<Rcpp::List>(view["space"])["dy"] != R_NilValue) {
        cv.dy(Rcpp::as<Rcpp::List>(view["space"])["dy"]);
      }
      if (Rcpp::as<Rcpp::List>(view["space"])["ny"] != R_NilValue) {
        cv.ny(Rcpp::as<Rcpp::List>(view["space"])["ny"]);
      }
      if (Rcpp::as<Rcpp::List>(view["space"])["srs"] != R_NilValue) {
        std::string srs = Rcpp::as<Rcpp::String>(Rcpp::as<Rcpp::List>(view["space"])["srs"]);
        cv.srs(srs);  
      }
      if (Rcpp::as<Rcpp::List>(view["time"])["t0"] != R_NilValue) {
        std::string tmp = Rcpp::as<Rcpp::String>(Rcpp::as<Rcpp::List>(view["time"])["t0"]);
        cv.t0(datetime::from_string(tmp));
      }
      if (Rcpp::as<Rcpp::List>(view["time"])["t1"] != R_NilValue) {
        std::string tmp = Rcpp::as<Rcpp::String>(Rcpp::as<Rcpp::List>(view["time"])["t1"]);
        cv.t1(datetime::from_string(tmp));
      }
      if (Rcpp::as<Rcpp::List>(view["time"])["nt"] != R_NilValue) {
        cv.nt(Rcpp::as<Rcpp::List>(view["time"])["nt"]);
      }
      if (Rcpp::as<Rcpp::List>(view["time"])["dt"] != R_NilValue) {
        std::string tmp = Rcpp::as<Rcpp::String>(Rcpp::as<Rcpp::List>(view["time"])["dt"]);
        cv.dt(duration::from_string(tmp));
        cv.t0().unit() = cv.dt().dt_unit; 
        cv.t1().unit() = cv.dt().dt_unit; 
      }
      
      if (view["aggregation"] != R_NilValue) {
        std::string tmp = Rcpp::as<Rcpp::String>(view["aggregation"]);
        cv.aggregation_method() = aggregation::from_string(tmp);
      }
      if (view["resampling"] != R_NilValue) {
        std::string tmp = Rcpp::as<Rcpp::String>(view["resampling"]);
        cv.resampling_method() = resampling::from_string(tmp);
      }
    
    std::shared_ptr<dummy_cube>* x = new std::shared_ptr<dummy_cube>( dummy_cube::create(cv, nbands, fill));
    (*x)->set_chunk_size(chunk_sizes[0], chunk_sizes[1], chunk_sizes[2]);
    Rcpp::XPtr< std::shared_ptr<dummy_cube> > p(x, true) ;
    return p;
  }
  catch (std::string s) {
    Rcpp::stop(s);
  }
}











// [[Rcpp::export]]
SEXP libgdalcubes_create_reduce_time_cube(SEXP pin, std::vector<std::string> reducers, std::vector<std::string> bands) {
  try {
    Rcpp::XPtr< std::shared_ptr<cube> > aa = Rcpp::as<Rcpp::XPtr<std::shared_ptr<cube>>>(pin);
    
    std::vector<std::pair<std::string, std::string>> reducer_bands;
    for (uint16_t i=0; i<reducers.size(); ++i) {
      // assuming reducers.size() == bands.size(), this is checked in R code calling this function
      reducer_bands.push_back(std::make_pair(reducers[i], bands[i]));
    }
    
    std::shared_ptr<reduce_time_cube>* x = new std::shared_ptr<reduce_time_cube>(reduce_time_cube::create(*aa, reducer_bands));
    Rcpp::XPtr< std::shared_ptr<reduce_time_cube> > p(x, true) ;
    
    return p;
    
  }
  catch (std::string s) {
    Rcpp::stop(s);
  }
}


// [[Rcpp::export]]
SEXP libgdalcubes_create_stream_reduce_time_cube(SEXP pin, std::string cmd, uint16_t nbands, std::vector<std::string> names) {
  try {
    Rcpp::XPtr< std::shared_ptr<cube> > aa = Rcpp::as<Rcpp::XPtr<std::shared_ptr<cube>>>(pin);
    std::shared_ptr<stream_reduce_time_cube>* x = new std::shared_ptr<stream_reduce_time_cube>(stream_reduce_time_cube::create(*aa, cmd, nbands, names));
    Rcpp::XPtr< std::shared_ptr<stream_reduce_time_cube> > p(x, true) ;
    return p;
  }
  catch (std::string s) {
    Rcpp::stop(s);
  }
}

// [[Rcpp::export]]
SEXP libgdalcubes_create_stream_reduce_space_cube(SEXP pin, std::string cmd, uint16_t nbands, std::vector<std::string> names) {
  try {
    Rcpp::XPtr< std::shared_ptr<cube> > aa = Rcpp::as<Rcpp::XPtr<std::shared_ptr<cube>>>(pin);
    std::shared_ptr<stream_reduce_space_cube>* x = new std::shared_ptr<stream_reduce_space_cube>(stream_reduce_space_cube::create(*aa, cmd, nbands, names));
    Rcpp::XPtr< std::shared_ptr<stream_reduce_space_cube> > p(x, true) ;
    return p;
  }
  catch (std::string s) {
    Rcpp::stop(s);
  }
}


// [[Rcpp::export]]
SEXP libgdalcubes_create_reduce_space_cube(SEXP pin, std::vector<std::string> reducers, std::vector<std::string> bands) {
  try {
    Rcpp::XPtr< std::shared_ptr<cube> > aa = Rcpp::as<Rcpp::XPtr<std::shared_ptr<cube>>>(pin);
    
    std::vector<std::pair<std::string, std::string>> reducer_bands;
    for (uint16_t i=0; i<reducers.size(); ++i) {
      // assuming reducers.size() == bands.size(), this is checked in R code calling this function
      reducer_bands.push_back(std::make_pair(reducers[i], bands[i]));
    }
    
    std::shared_ptr<reduce_space_cube>* x = new std::shared_ptr<reduce_space_cube>(reduce_space_cube::create(*aa, reducer_bands));
    Rcpp::XPtr< std::shared_ptr<reduce_space_cube> > p(x, true) ;
    
    return p;
    
  }
  catch (std::string s) {
    Rcpp::stop(s);
  }
}


// [[Rcpp::export]]
SEXP libgdalcubes_create_window_time_cube_reduce(SEXP pin, std::vector<int> window, std::vector<std::string> reducers, std::vector<std::string> bands) {
  try {
    Rcpp::XPtr< std::shared_ptr<cube> > aa = Rcpp::as<Rcpp::XPtr<std::shared_ptr<cube>>>(pin);
    
    std::vector<std::pair<std::string, std::string>> reducer_bands;
    for (uint16_t i=0; i<reducers.size(); ++i) {
      // assuming reducers.size() == bands.size(), this is checked in R code calling this function
      reducer_bands.push_back(std::make_pair(reducers[i], bands[i]));
    }
    
    std::shared_ptr<window_time_cube>* x = new std::shared_ptr<window_time_cube>(window_time_cube::create(*aa, reducer_bands, window[0], window[1]));
    Rcpp::XPtr< std::shared_ptr<window_time_cube> > p(x, true) ;
    
    return p;
    
  }
  catch (std::string s) {
    Rcpp::stop(s);
  }
}

// [[Rcpp::export]]
SEXP libgdalcubes_create_window_time_cube_kernel(SEXP pin, std::vector<int> window, std::vector<double> kernel) {
  try {
    Rcpp::XPtr< std::shared_ptr<cube> > aa = Rcpp::as<Rcpp::XPtr<std::shared_ptr<cube>>>(pin);
    
    std::shared_ptr<window_time_cube>* x = new std::shared_ptr<window_time_cube>(window_time_cube::create(*aa, kernel, window[0], window[1]));
    Rcpp::XPtr< std::shared_ptr<window_time_cube> > p(x, true) ;
    return p;
    
  }
  catch (std::string s) {
    Rcpp::stop(s);
  }
}


// [[Rcpp::export]]
SEXP libgdalcubes_create_join_bands_cube(Rcpp::List pin_list,  std::vector<std::string> cube_names) {
  try {
    std::vector< std::shared_ptr<cube> > cube_list;
    for (uint16_t i=0; i<pin_list.size(); ++i) {
      cube_list.push_back(*(Rcpp::as<Rcpp::XPtr<std::shared_ptr<cube>>>(pin_list[i])));
    }
  
    std::shared_ptr<join_bands_cube>* x = new std::shared_ptr<join_bands_cube>(join_bands_cube::create(cube_list, cube_names));
    Rcpp::XPtr< std::shared_ptr<join_bands_cube> > p(x, true) ;
    
    return p;
    
  }
  catch (std::string s) {
    Rcpp::stop(s);
  }
}

// [[Rcpp::export]]
SEXP libgdalcubes_create_select_bands_cube(SEXP pin, std::vector<std::string> bands) {
  try {
    Rcpp::XPtr< std::shared_ptr<cube> > aa = Rcpp::as<Rcpp::XPtr<std::shared_ptr<cube>>>(pin);
    
    std::shared_ptr<select_bands_cube>* x = new std::shared_ptr<select_bands_cube>(select_bands_cube::create(*aa, bands));
    Rcpp::XPtr< std::shared_ptr<select_bands_cube> > p(x, true) ;
    
    return p;
    
  }
  catch (std::string s) {
    Rcpp::stop(s);
  }
}


// [[Rcpp::export]]
SEXP libgdalcubes_create_select_time_cube(SEXP pin, std::vector<std::string> t) {
  try {
    Rcpp::XPtr< std::shared_ptr<cube> > aa = Rcpp::as<Rcpp::XPtr<std::shared_ptr<cube>>>(pin);
    
      std::shared_ptr<select_time_cube>* x = new std::shared_ptr<select_time_cube>(select_time_cube::create(*aa, t));
    Rcpp::XPtr< std::shared_ptr<select_time_cube> > p(x, true) ;
    
    return p;
  }
  catch (std::string s) {
    Rcpp::stop(s);
  }
}

// [[Rcpp::export]]
SEXP libgdalcubes_create_apply_pixel_cube(SEXP pin, std::vector<std::string> expr, std::vector<std::string> names, bool keep_bands=false) {
  try {
    Rcpp::XPtr< std::shared_ptr<cube> > aa = Rcpp::as<Rcpp::XPtr<std::shared_ptr<cube>>>(pin);
    
    std::shared_ptr<apply_pixel_cube>* x = new std::shared_ptr<apply_pixel_cube>(apply_pixel_cube::create(*aa, expr, names, keep_bands));
    Rcpp::XPtr< std::shared_ptr<apply_pixel_cube> > p(x, true) ;
    
    return p;
    
  }
  catch (std::string s) {
    Rcpp::stop(s);
  }
}

// [[Rcpp::export]]
SEXP libgdalcubes_create_stream_apply_pixel_cube(SEXP pin, std::string cmd, uint16_t nbands, std::vector<std::string> names, bool keep_bands = false) {
  try {
    Rcpp::XPtr< std::shared_ptr<cube> > aa = Rcpp::as<Rcpp::XPtr<std::shared_ptr<cube>>>(pin);
    std::shared_ptr<stream_apply_pixel_cube>* x = new std::shared_ptr<stream_apply_pixel_cube>(stream_apply_pixel_cube::create(*aa, cmd, nbands, names, keep_bands));
    Rcpp::XPtr< std::shared_ptr<stream_apply_pixel_cube> > p(x, true) ;
    return p;
  }
  catch (std::string s) {
    Rcpp::stop(s);
  }
}

// [[Rcpp::export]]
SEXP libgdalcubes_create_stream_apply_time_cube(SEXP pin, std::string cmd, uint16_t nbands, std::vector<std::string> names, bool keep_bands = false) {
  try {
    Rcpp::XPtr< std::shared_ptr<cube> > aa = Rcpp::as<Rcpp::XPtr<std::shared_ptr<cube>>>(pin);
    std::shared_ptr<stream_apply_time_cube>* x = new std::shared_ptr<stream_apply_time_cube>(stream_apply_time_cube::create(*aa, cmd, nbands, names, keep_bands));
    Rcpp::XPtr< std::shared_ptr<stream_apply_time_cube> > p(x, true) ;
    return p;
  }
  catch (std::string s) {
    Rcpp::stop(s);
  }
}



// [[Rcpp::export]]
SEXP libgdalcubes_create_filter_predicate_cube(SEXP pin, std::string pred) {
  try {
    Rcpp::XPtr< std::shared_ptr<cube> > aa = Rcpp::as<Rcpp::XPtr<std::shared_ptr<cube>>>(pin);
    std::shared_ptr<filter_pixel_cube>* x = new std::shared_ptr<filter_pixel_cube>(filter_pixel_cube::create(*aa, pred));
    Rcpp::XPtr< std::shared_ptr<filter_pixel_cube> > p(x, true) ;
    return p;
  }
  catch (std::string s) {
    Rcpp::stop(s);
  }
}

// [[Rcpp::export]]
SEXP libgdalcubes_create_filter_geom_cube(SEXP pin, std::string wkt, std::string srs) {
  try {
    Rcpp::XPtr< std::shared_ptr<cube> > aa = Rcpp::as<Rcpp::XPtr<std::shared_ptr<cube>>>(pin);
    std::shared_ptr<filter_geom_cube>* x = new std::shared_ptr<filter_geom_cube>(filter_geom_cube::create(*aa, wkt, srs));
    Rcpp::XPtr< std::shared_ptr<filter_geom_cube> > p(x, true) ;
    return p;
  }
  catch (std::string s) {
    Rcpp::stop(s);
  }
}


// [[Rcpp::export]]
void libgdalcubes_debug_output( bool debug) {
  try {
    if (debug) {
      config::instance()->set_error_handler(error_handling_r::debug); 
    }
    else {
      config::instance()->set_error_handler(error_handling_r::standard); 
    }
  }
  catch (std::string s) {
    Rcpp::stop(s);
  }
}




// [[Rcpp::export]]
void libgdalcubes_eval_cube( SEXP pin, std::string outfile, uint8_t compression_level=0, bool with_VRT=false, 
                             bool write_bounds = true,  SEXP packing = R_NilValue) {
  try {
    Rcpp::XPtr< std::shared_ptr<cube> > aa = Rcpp::as<Rcpp::XPtr< std::shared_ptr<cube> >>(pin);
    packed_export p = packed_export::make_none();
    if (packing != R_NilValue) {
      
      std::string type = Rcpp::as<Rcpp::List>(packing)["type"];
      
      if (type == "uint8") {
        p.type = packed_export::packing_type::PACK_UINT8;
      }
      else if (type == "uint16") {
        p.type = packed_export::packing_type::PACK_UINT16;
      }
      else if (type == "uint32") {
        p.type = packed_export::packing_type::PACK_UINT32;
      }
      else if (type == "int16") {
        p.type = packed_export::packing_type::PACK_INT16;
      }
      else if (type == "int32") {
        p.type = packed_export::packing_type::PACK_INT32;
      }
      else {
        // TODO: NO PACKING-> WARNING
      }
      p.offset = Rcpp::as<std::vector<double>>(Rcpp::as<Rcpp::List>(packing)["offset"]);
      p.scale = Rcpp::as<std::vector<double>>(Rcpp::as<Rcpp::List>(packing)["scale"]);
      p.nodata = Rcpp::as<std::vector<double>>(Rcpp::as<Rcpp::List>(packing)["nodata"]);
    }
    (*aa)->write_netcdf_file(outfile, compression_level, with_VRT, write_bounds, p);
  }
  catch (std::string s) {
    Rcpp::stop(s);
  }
}

// [[Rcpp::export]]
void libgdalcubes_write_chunks_ncdf( SEXP pin, std::string dir, std::string name, uint8_t compression_level=0) {
  try {
    Rcpp::XPtr< std::shared_ptr<cube> > aa = Rcpp::as<Rcpp::XPtr< std::shared_ptr<cube> >>(pin);
    (*aa)->write_chunks_netcdf(dir, name, compression_level);
  }
  catch (std::string s) {
    Rcpp::stop(s);
  }
}


// [[Rcpp::export]]
void libgdalcubes_write_tif( SEXP pin, std::string dir, std::string prefix="", 
                             bool overviews = false, bool cog = false, 
                             SEXP creation_options = R_NilValue,
                             std::string rsmpl_overview = "nearest", 
                             SEXP packing = R_NilValue) {
  try {
    Rcpp::XPtr< std::shared_ptr<cube> > aa = Rcpp::as<Rcpp::XPtr< std::shared_ptr<cube> >>(pin);
    std::map<std::string, std::string> co;
    
    if (creation_options != R_NilValue) {
      Rcpp::List colist = Rcpp::as<Rcpp::List>(creation_options);
      Rcpp::CharacterVector names = colist.names();
      for(uint16_t i=0; i< names.size(); ++i) {
        std::string key = Rcpp::as<std::string>(names[i]);
        std::string value =Rcpp::as<std::string>(Rcpp::as<Rcpp::CharacterVector>(colist[i]));
        co[key] = value;
      }
    }
    
    packed_export p = packed_export::make_none();
    if (packing != R_NilValue) {
      
      std::string type = Rcpp::as<Rcpp::List>(packing)["type"];
      
      if (type == "uint8") {
        p.type = packed_export::packing_type::PACK_UINT8;
      }
      else if (type == "uint16") {
        p.type = packed_export::packing_type::PACK_UINT16;
      }
      else if (type == "uint32") {
        p.type = packed_export::packing_type::PACK_UINT32;
      }
      else if (type == "int16") {
        p.type = packed_export::packing_type::PACK_INT16;
      }
      else if (type == "int32") {
        p.type = packed_export::packing_type::PACK_INT32;
      }
      else {
        // TODO: NO PACKING-> WARNING
      }
      
      p.offset = Rcpp::as<std::vector<double>>(Rcpp::as<Rcpp::List>(packing)["offset"]);
      p.scale = Rcpp::as<std::vector<double>>(Rcpp::as<Rcpp::List>(packing)["scale"]);
      p.nodata = Rcpp::as<std::vector<double>>(Rcpp::as<Rcpp::List>(packing)["nodata"]);
    }
    
    (*aa)->write_tif_collection(dir, prefix, overviews, cog, co, rsmpl_overview, p);
    
    
  }
  catch (std::string s) {
    Rcpp::stop(s);
  }
}


// [[Rcpp::export]]
SEXP libgdalcubes_create_stream_cube(SEXP pin, std::string cmd) {
  try {
    Rcpp::XPtr< std::shared_ptr<cube> > aa = Rcpp::as<Rcpp::XPtr< std::shared_ptr<cube> >>(pin);
    
    std::shared_ptr<stream_cube>* x = new std::shared_ptr<stream_cube>( stream_cube::create(*aa, cmd, true));
    
    Rcpp::XPtr< std::shared_ptr<stream_cube> > p(x, true) ;
  
    return p;
  }
  catch (std::string s) {
    Rcpp::stop(s);
  }
}

// [[Rcpp::export]]
SEXP libgdalcubes_create_fill_time_cube(SEXP pin, std::string method) {
  try {
    Rcpp::XPtr< std::shared_ptr<cube> > aa = Rcpp::as<Rcpp::XPtr< std::shared_ptr<cube> >>(pin);
    std::shared_ptr<fill_time_cube>* x = new std::shared_ptr<fill_time_cube>( fill_time_cube::create(*aa, method));
    Rcpp::XPtr< std::shared_ptr<fill_time_cube> > p(x, true) ;
    return p;
  } 
  catch (std::string s) {
    Rcpp::stop(s);
  }
}



// [[Rcpp::export]]
SEXP libgdalcubes_query_points(SEXP pin, std::vector<double> px, std::vector<double> py, std::vector<std::string> pt, std::string srs) {
  try {
    Rcpp::XPtr< std::shared_ptr<cube> > aa = Rcpp::as<Rcpp::XPtr< std::shared_ptr<cube> >>(pin);
    std::vector<std::vector<double>> res = vector_queries::query_points(*aa, px, py, pt, srs);
    Rcpp::List df(res.size());
  
    for (uint16_t i=0; i<res.size(); ++i) {
     df[i] = res[i];
    }
    return df;
  }
  catch (std::string s) {
    Rcpp::stop(s);
  }
}


// [[Rcpp::export]]
SEXP libgdalcubes_query_timeseries(SEXP pin, std::vector<double> px, std::vector<double> py, std::string srs) {
  try {
    Rcpp::XPtr< std::shared_ptr<cube> > aa = Rcpp::as<Rcpp::XPtr< std::shared_ptr<cube> >>(pin);
    std::vector<std::vector<std::vector<double>>> res = vector_queries::query_timeseries(*aa, px, py, srs);
    Rcpp::List dflist(res.size());
    
    for (uint16_t i=0; i<res.size(); ++i) {
      dflist[i] = res[i];
    }
    return dflist;
  }
  catch (std::string s) {
    Rcpp::stop(s);
  }
}

// [[Rcpp::export]]
void libgdalcubes_zonal_statistics(SEXP pin, std::string ogr_dataset, std::vector<std::string> agg_funcs, std::vector<std::string> agg_bands, std::string out_path, bool overwrite, std::string ogr_layer) {
 try {
    Rcpp::XPtr< std::shared_ptr<cube> > aa = Rcpp::as<Rcpp::XPtr< std::shared_ptr<cube> >>(pin);

    std::vector<std::pair<std::string, std::string>> agg;
    // assumption: agg_funcs.size() == agg_bands.size 
    for (uint32_t i=0; i<agg_funcs.size(); ++i) {
      agg.push_back(std::make_pair(agg_funcs[i],agg_bands[i]));
    }
   vector_queries::zonal_statistics(*aa, ogr_dataset, agg, out_path, overwrite, ogr_layer);
  }
  catch (std::string s) {
    Rcpp::stop(s);
  } 
}


// [[Rcpp::export]]
void libgdalcubes_set_threads(IntegerVector n) {
  config::instance()->set_default_chunk_processor(std::dynamic_pointer_cast<chunk_processor>(std::make_shared<chunk_processor_multithread_interruptible>(n[0])));
}


// [[Rcpp::export]]
void libgdalcubes_set_use_overviews(bool use_overviews) {
  config::instance()->set_gdal_use_overviews(use_overviews);
}


// [[Rcpp::export]]
std::string libgdalcubes_translate_cog(SEXP collection, std::string out_dir, uint16_t nthreads, bool overwrite, std::vector<std::string> creation_options) {
  try {
    Rcpp::XPtr<std::shared_ptr<image_collection>> aa = Rcpp::as<Rcpp::XPtr<std::shared_ptr<image_collection>>>(collection);
    image_collection_ops::translate_cog(*aa, out_dir, nthreads, overwrite, creation_options);
    return filesystem::join(out_dir, filesystem::filename((*aa)->get_filename()));
  }
  catch (std::string s) {
    Rcpp::stop(s);
  } 
  
}


// [[Rcpp::export]]
std::string libgdalcubes_translate_gtiff(SEXP collection, std::string out_dir, uint16_t nthreads, bool overwrite, std::vector<std::string> creation_options) {
  try {
    Rcpp::XPtr<std::shared_ptr<image_collection>> aa = Rcpp::as<Rcpp::XPtr<std::shared_ptr<image_collection>>>(collection);
    image_collection_ops::translate_gtiff(*aa, out_dir, nthreads, overwrite, creation_options);
    return filesystem::join(out_dir, filesystem::filename((*aa)->get_filename()));
  }
  catch (std::string s) {
    Rcpp::stop(s);
  } 
  
}


// [[Rcpp::export]]
void libgdalcubes_set_swarm(std::vector<std::string> swarm) {
  auto p = gdalcubes_swarm::from_urls(swarm);
  //p->set_threads(nthreads);
  config::instance()->set_default_chunk_processor(p);
  // TODO: how to make swarm interruptible
}

// [[Rcpp::export]]
std::string libgdalcubes_simple_hash(std::string instr) {
  return utils::hash(instr);
}



