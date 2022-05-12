#include "gdalcubes/src/gdalcubes.h"
#include "gdalcubes/src/cube_factory.h"
#include "multiprocess.h"
#include "error.h"

// [[Rcpp::plugins("cpp11")]]
#include <Rcpp.h>
#include <memory>
#include <thread>
#include <algorithm>
//#include <fstream>


using namespace Rcpp;
using namespace gdalcubes;



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
    _set(1);
    std::stringstream ss;
    ss << std::endl;
    r_stderr_buf::print(ss.str());
    error_handling_r::do_output();
    _m.unlock();
  }
  

  
  progress_simple_R() : _p(0) {}
  
  ~progress_simple_R(){}
  
private:
  
  std::mutex _m;
  double _p;
  
  void _set(double p) { // call this function only with a lock on _m
    error_handling_r::defer_output();
    _p = p;
    
    std::stringstream s;
    s << "[";
    int pp = 50 * _p;
    int i = 0;
    while (i < pp) {
      if (i < pp) s << "=";
      ++i;
    }
    s << ">";
    ++i;
    while (i < 50) {
      s << " ";
      ++i;
    }
    s << "] " << int(p * 100.0) << " %\r";
    r_stderr_buf::print(s.str());
  }
};

struct progress_none_R : public progress {
  std::shared_ptr<progress> get() override { return std::make_shared<progress_none_R>(); }
  void set(double p) override {};
  void increment(double dp) override {}
  virtual void finalize() override {}
  progress_none_R() {}
};




cube_view cube_view_from_list(SEXP v) {
  Rcpp::List view = Rcpp::as<Rcpp::List>(v);
  cube_view cv;
  
  // x
  Rcpp::RObject dx = Rcpp::as<Rcpp::List>(view["space"])["dx"];
  Rcpp::RObject nx = Rcpp::as<Rcpp::List>(view["space"])["nx"];
  Rcpp::RObject left = Rcpp::as<Rcpp::List>(view["space"])["left"];
  Rcpp::RObject right = Rcpp::as<Rcpp::List>(view["space"])["right"];
  if (!left.isNULL() && !right.isNULL() && !dx.isNULL() && !nx.isNULL()) {
    //GCBS_WARN("Expected only three out of [left, right, dx, nx], ignoring dx");
    cv.set_x_axis(double(Rcpp::as<Rcpp::NumericVector>(left)[0]), 
                  double(Rcpp::as<Rcpp::NumericVector>(right)[0]), 
                  uint32_t(Rcpp::as<Rcpp::IntegerVector>(nx)[0]));
  }
  else if (!left.isNULL() && !right.isNULL() && !dx.isNULL() && nx.isNULL()) {  
    cv.set_x_axis(double(Rcpp::as<Rcpp::NumericVector>(left)[0]), 
                  double(Rcpp::as<Rcpp::NumericVector>(right)[0]), 
                  double(Rcpp::as<Rcpp::NumericVector>(dx)[0]));
  }
  else if (!left.isNULL()  && !right.isNULL() && dx.isNULL() && !nx.isNULL()) {
    cv.set_x_axis(double(Rcpp::as<Rcpp::NumericVector>(left)[0]), 
                  double(Rcpp::as<Rcpp::NumericVector>(right)[0]), 
                  uint32_t(Rcpp::as<Rcpp::IntegerVector>(nx)[0]));
  }
  else if (!left.isNULL()  && right.isNULL() && !dx.isNULL() && !nx.isNULL()) {
    cv.set_x_axis(double(Rcpp::as<Rcpp::NumericVector>(left)[0]), 
                  uint32_t(Rcpp::as<Rcpp::IntegerVector>(nx)[0]), 
                  double(Rcpp::as<Rcpp::NumericVector>(dx)[0]));
  }
  else if (left.isNULL()  && !right.isNULL() && !dx.isNULL() && !nx.isNULL()) {
    cv.set_x_axis(uint32_t(Rcpp::as<Rcpp::IntegerVector>(nx)[0]), 
                  double(Rcpp::as<Rcpp::NumericVector>(right)[0]), 
                  double(Rcpp::as<Rcpp::NumericVector>(dx)[0]));
  }
  else {
    GCBS_ERROR("Specification of x dimension is incomplete");
    throw std::string("Specification of x dimension is incomplete");
  }
  
  
  // y
  Rcpp::RObject dy = Rcpp::as<Rcpp::List>(view["space"])["dy"];
  Rcpp::RObject ny = Rcpp::as<Rcpp::List>(view["space"])["ny"];
  Rcpp::RObject bottom = Rcpp::as<Rcpp::List>(view["space"])["bottom"];
  Rcpp::RObject top = Rcpp::as<Rcpp::List>(view["space"])["top"];
  if (!bottom.isNULL() && !top.isNULL() && !dy.isNULL() && !ny.isNULL()) {
    //GCBS_WARN("Expected only three out of [bottom, top, ny, dy], ignoring dy");
    cv.set_y_axis(double(Rcpp::as<Rcpp::NumericVector>(bottom)[0]), 
                  double(Rcpp::as<Rcpp::NumericVector>(top)[0]), 
                  uint32_t(Rcpp::as<Rcpp::IntegerVector>(ny)[0]));
  }
  else if (!bottom.isNULL() && !top.isNULL() && !dy.isNULL() && ny.isNULL()) {  
    cv.set_y_axis(double(Rcpp::as<Rcpp::NumericVector>(bottom)[0]), 
                  double(Rcpp::as<Rcpp::NumericVector>(top)[0]), 
                  double(Rcpp::as<Rcpp::NumericVector>(dy)[0]));
  }
  else if (!bottom.isNULL()  && !top.isNULL() && dy.isNULL() && !ny.isNULL()) {
    cv.set_y_axis(double(Rcpp::as<Rcpp::NumericVector>(bottom)[0]), 
                  double(Rcpp::as<Rcpp::NumericVector>(top)[0]), 
                  uint32_t(Rcpp::as<Rcpp::IntegerVector>(ny)[0]));
  }
  else if (!bottom.isNULL()  && top.isNULL() && !dy.isNULL() && !ny.isNULL()) {
    cv.set_y_axis(double(Rcpp::as<Rcpp::NumericVector>(bottom)[0]), 
                  uint32_t(Rcpp::as<Rcpp::IntegerVector>(ny)[0]), 
                  double(Rcpp::as<Rcpp::NumericVector>(dy)[0]));
  }
  else if (bottom.isNULL()  && !top.isNULL() && !dy.isNULL() && !ny.isNULL()) {
    cv.set_y_axis(uint32_t(Rcpp::as<Rcpp::IntegerVector>(ny)[0]), 
                  double(Rcpp::as<Rcpp::NumericVector>(top)[0]), 
                  double(Rcpp::as<Rcpp::NumericVector>(dy)[0]));
  }
  else {
    GCBS_ERROR("Specification of y dimension is incomplete");
    throw std::string("Specification of y dimension is incomplete");
  }
  
  
  // t
  Rcpp::RObject t0 = Rcpp::as<Rcpp::List>(view["time"])["t0"];
  Rcpp::RObject t1 = Rcpp::as<Rcpp::List>(view["time"])["t1"];
  Rcpp::RObject nt = Rcpp::as<Rcpp::List>(view["time"])["nt"];
  Rcpp::RObject dt = Rcpp::as<Rcpp::List>(view["time"])["dt"];
  
  if (!t0.isNULL() && !t1.isNULL() && !dt.isNULL() && !nt.isNULL()) {
    //GCBS_WARN("Expected only three out of [t0, t1, nt, dt], ignoring nt");
    cv.set_t_axis(datetime::from_string(Rcpp::as<Rcpp::String>(t0)), 
                  datetime::from_string(Rcpp::as<Rcpp::String>(t1)), 
                  duration::from_string(Rcpp::as<Rcpp::String>(dt)));
  }
  else if (!t0.isNULL() && !t1.isNULL() && !dt.isNULL() && nt.isNULL()) {  
    cv.set_t_axis(datetime::from_string(Rcpp::as<Rcpp::String>(t0)), 
                  datetime::from_string(Rcpp::as<Rcpp::String>(t1)), 
                  duration::from_string(Rcpp::as<Rcpp::String>(dt)));
  }
  else if (!t0.isNULL()  && !t1.isNULL() && dt.isNULL() && !nt.isNULL()) {
    cv.set_t_axis(datetime::from_string(Rcpp::as<Rcpp::String>(t0)), 
                  datetime::from_string(Rcpp::as<Rcpp::String>(t1)), 
                  uint32_t(Rcpp::as<Rcpp::IntegerVector>(nt)[0]));
  }
  // else if (!t0.isNULL()  && t1.isNULL() && !dt.isNULL() && !nt.isNULL()) {
  //   cv.set_t_axis(double(Rcpp::as<Rcpp::NumericVector>(t0)[0]), 
  //                 uint32_t(Rcpp::as<Rcpp::IntegerVector>(nt)[0]), 
  //                 double(Rcpp::as<Rcpp::NumericVector>(dt)[0]));
  // }
  // else if (t0.isNULL()  && !t1.isNULL() && !dt.isNULL() && !nt.isNULL()) {
  //   cv.set_t_axis(uint32_t(Rcpp::as<Rcpp::IntegerVector>(nt)[0]), 
  //                 double(Rcpp::as<Rcpp::NumericVector>(t1)[0]), 
  //                 double(Rcpp::as<Rcpp::NumericVector>(dt)[0]));
  // }
  else {
    GCBS_ERROR("Specification of t dimension is incomplete");
    throw std::string("Specification of t dimension is incomplete");
  }
  
  
  if (Rcpp::as<Rcpp::List>(view["space"])["srs"] != R_NilValue) {
    std::string srs = Rcpp::as<Rcpp::String>(Rcpp::as<Rcpp::List>(view["space"])["srs"]);
    cv.srs(srs);  
  }
  else {
    GCBS_ERROR("Missing spatial reference system");
    throw std::string("Missing spatial reference system");
  }
  if (view["aggregation"] != R_NilValue) {
    std::string tmp = Rcpp::as<Rcpp::String>(view["aggregation"]);
    cv.aggregation_method() = aggregation::from_string(tmp);
  }
  if (view["resampling"] != R_NilValue) {
    std::string tmp = Rcpp::as<Rcpp::String>(view["resampling"]);
    cv.resampling_method() = resampling::from_string(tmp);
  }
  return cv;
}






// see https://stackoverflow.com/questions/26666614/how-do-i-check-if-an-externalptr-is-null-from-within-r
// [[Rcpp::export]]
Rcpp::LogicalVector gc_is_null(SEXP pointer) {
  return Rcpp::LogicalVector(!R_ExternalPtrAddr(pointer));
}


// [[Rcpp::export]]
Rcpp::List gc_version() {
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
std::vector<std::string> gc_gdalformats() {
  return config::instance()->gdal_formats();
}

// [[Rcpp::export]]
void gc_set_gdal_config(std::string k, std::string v) {
  config::instance()->set_gdal_option(k, v);
}

// [[Rcpp::export]]
void gc_set_streamining_dir(std::string dir) {
  config::instance()->set_streaming_dir(dir);
}


// [[Rcpp::export]]
std::string gc_gdalversion() {
  return config::instance()->gdal_version_info();
}

// [[Rcpp::export]]
bool gc_gdal_has_geos() {
  return config::instance()->gdal_has_geos();
}


// [[Rcpp::export]]
void gc_add_format_dir(std::string dir) {
  config::instance()->add_collection_format_preset_dir(dir);
}

// [[Rcpp::export]]
void gc_init() {
  config::instance()->gdalcubes_init();
  config::instance()->set_default_progress_bar(std::make_shared<progress_simple_R>());
  config::instance()->set_error_handler(error_handling_r::standard); 
  
  // Interruptible chunk processor
  config::instance()->set_gdal_option("GDAL_NUM_THREADS", "ALL_CPUS");
  
}

// [[Rcpp::export]]
void gc_cleanup() {
  config::instance()->gdalcubes_cleanup();
}












// This function is also covered by gc_dimension_values() and should be deprecated
// [[Rcpp::export]]
Rcpp::StringVector gc_datetime_values(SEXP pin) {
  Rcpp::XPtr<std::shared_ptr<cube>> aa = Rcpp::as<Rcpp::XPtr<std::shared_ptr<cube>>>(pin);
  std::shared_ptr<cube> x = *aa;
  
  Rcpp::CharacterVector out(x->size_t());
  
  for (uint32_t i = 0; i < x->size_t(); ++i) {
    out[i] = x->st_reference()->datetime_at_index(i).to_string();
  }
  return out;
}

// [[Rcpp::export]]
Rcpp::List gc_cube_info( SEXP pin) {

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
Rcpp::List gc_dimension_values_from_view(Rcpp::List view, std::string dt_unit="") {
  
  cube_view cv = cube_view_from_list(view);
  
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
Rcpp::List gc_dimension_bounds(SEXP pin, std::string dt_unit="") {
  
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
Rcpp::List gc_dimension_values(SEXP pin, std::string dt_unit="") {
  
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
SEXP gc_open_image_collection(std::string filename) {
  
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
Rcpp::List gc_image_collection_info( SEXP pin) {
  
  try {
    Rcpp::XPtr<std::shared_ptr<image_collection>> aa = Rcpp::as<Rcpp::XPtr<std::shared_ptr<image_collection>>>(pin);
    std::shared_ptr<image_collection> ic = *aa;
    
    if ((*aa)->is_empty()) {
      return Rcpp::List::create();
    }
    
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
Rcpp::List gc_image_collection_extent( SEXP pin, std::string srs) {
  
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
void gc_create_image_collection_from_format(std::vector<std::string> files, std::string format_file, std::string outfile, bool unroll_archives=true) {
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
void gc_create_image_collection_from_datetime(std::string outfile, std::vector<std::string> files, 
                                                        std::vector<std::string> date_time, bool use_subdatasets, 
                                                        std::vector<std::string> band_names, bool one_band_per_file) {
  try {
    image_collection::create(files, date_time, band_names, use_subdatasets, one_band_per_file)->write(outfile);
  }
  catch (std::string s) {
    Rcpp::stop(s);
  }
}


// [[Rcpp::export]]
void gc_add_images(SEXP pin, std::vector<std::string> files, bool unroll_archives=true, std::string outfile = "") {
  
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
SEXP gc_list_collection_formats() {
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
SEXP gc_create_view(SEXP v) {
  cube_view cv = cube_view_from_list(v);
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
SEXP gc_create_image_collection_cube(SEXP pin, Rcpp::IntegerVector chunk_sizes, SEXP mask, SEXP v = R_NilValue) {

  try {
    Rcpp::XPtr<std::shared_ptr<image_collection>> aa = Rcpp::as<Rcpp::XPtr<std::shared_ptr<image_collection>>>(pin);
    
    std::shared_ptr<image_collection_cube>* x;
    if (v == R_NilValue) {
      x = new std::shared_ptr<image_collection_cube>( image_collection_cube::create(*aa));
    }
    else {
      cube_view cv = cube_view_from_list(v);
      x = new std::shared_ptr<image_collection_cube>( image_collection_cube::create(*aa, cv));
    }
    (*x)->set_chunk_size(chunk_sizes[0], chunk_sizes[1], chunk_sizes[2]);
    
    
    if (mask != R_NilValue) {
      std::string band_name = Rcpp::as<Rcpp::List>(mask)["band"]; 
      bool invert = Rcpp::as<Rcpp::List>(mask)["invert"];
      
      if (Rcpp::as<Rcpp::List>(mask).containsElementNamed("values") && Rcpp::as<Rcpp::List>(mask)["values"] != R_NilValue) {
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
SEXP gc_create_ncdf_cube(std::string path, Rcpp::IntegerVector chunk_sizes, bool auto_unpack) {
  
  try {
    std::shared_ptr<ncdf_cube>* x  = new std::shared_ptr<ncdf_cube>( ncdf_cube::create(path, auto_unpack));
    if (chunk_sizes.size() == 3) {
      (*x)->set_chunk_size(chunk_sizes[0], chunk_sizes[1], chunk_sizes[2]);
    }
    Rcpp::XPtr< std::shared_ptr<ncdf_cube> > p(x, true) ;
    return p;
  }
  catch (std::string s) {
    Rcpp::stop(s);
  }
}









// [[Rcpp::export]]
SEXP gc_create_dummy_cube(SEXP v, uint16_t nbands, double fill, Rcpp::IntegerVector chunk_sizes) {
  try {
    Rcpp::List view = Rcpp::as<Rcpp::List>(v);
    cube_view cv = cube_view_from_list(v);
    
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
SEXP gc_create_empty_cube(SEXP v, uint16_t nbands, Rcpp::IntegerVector chunk_sizes) {
  try {
    Rcpp::List view = Rcpp::as<Rcpp::List>(v);
    cube_view cv = cube_view_from_list(v);
    
    std::shared_ptr<empty_cube>* x = new std::shared_ptr<empty_cube>(empty_cube::create(cv, nbands));
    (*x)->set_chunk_size(chunk_sizes[0], chunk_sizes[1], chunk_sizes[2]);
    Rcpp::XPtr< std::shared_ptr<empty_cube> > p(x, true) ;
    return p;
  }
  catch (std::string s) {
    Rcpp::stop(s);
  }
}



// [[Rcpp::export]]
SEXP gc_copy_cube(SEXP pin) {
  try {
    Rcpp::XPtr< std::shared_ptr<cube> > aa = Rcpp::as<Rcpp::XPtr<std::shared_ptr<cube>>>(pin);
    std::shared_ptr<cube>* x  = new std::shared_ptr<cube>(cube_factory::instance()->create_from_json((*aa)->make_constructible_json()));
    Rcpp::XPtr< std::shared_ptr<cube> > p(x, true) ;
    return p;
  }
  catch (std::string s) {
    Rcpp::stop(s);
  }
}


// [[Rcpp::export]]
SEXP gc_from_json_file(std::string path) {
  try {
    std::shared_ptr<cube>* x  = new std::shared_ptr<cube>(cube_factory::instance()->create_from_json_file(path));
    Rcpp::XPtr< std::shared_ptr<cube> > p(x, true) ;
    return p;
  }
  catch (std::string s) {
    Rcpp::stop(s);
  }
}

// [[Rcpp::export]]
SEXP gc_from_json_string(std::string json) {
  try {
    std::string err;
    json11::Json j = json11::Json::parse(json, err);
    if (!err.empty()) {
      Rcpp::stop(err);
    }
    std::shared_ptr<cube>* x  = new std::shared_ptr<cube>(cube_factory::instance()->create_from_json(j));
    Rcpp::XPtr< std::shared_ptr<cube> > p(x, true);
    return p;
  }
  catch (std::string s) {
    Rcpp::stop(s);
  }
}







// [[Rcpp::export]]
SEXP gc_create_rename_bands_cube(SEXP pin, std::vector<std::string> names_old, std::vector<std::string> names_new) {
  try {
    Rcpp::XPtr< std::shared_ptr<cube> > aa = Rcpp::as<Rcpp::XPtr<std::shared_ptr<cube>>>(pin);
    
    std::map<std::string, std::string> bandnames;
    for (uint16_t i=0; i<names_old.size(); ++i) {
      bandnames[names_old[i]] = names_new[i];
    }

    std::shared_ptr<rename_bands_cube>* x = new std::shared_ptr<rename_bands_cube>(rename_bands_cube::create(*aa, bandnames));
    Rcpp::XPtr< std::shared_ptr<rename_bands_cube> > p(x, true) ;
    
    return p;
  }
  catch (std::string s) {
    Rcpp::stop(s);
  }
}
  


// [[Rcpp::export]]
SEXP gc_create_reduce_time_cube(SEXP pin, std::vector<std::string> reducers, std::vector<std::string> bands) {
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
SEXP gc_create_stream_reduce_time_cube(SEXP pin, std::string cmd, uint16_t nbands, std::vector<std::string> names) {
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
SEXP gc_create_stream_reduce_space_cube(SEXP pin, std::string cmd, uint16_t nbands, std::vector<std::string> names) {
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
SEXP gc_create_reduce_space_cube(SEXP pin, std::vector<std::string> reducers, std::vector<std::string> bands) {
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
SEXP gc_create_window_time_cube_reduce(SEXP pin, std::vector<int> window, std::vector<std::string> reducers, std::vector<std::string> bands) {
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
SEXP gc_create_window_time_cube_kernel(SEXP pin, std::vector<int> window, std::vector<double> kernel) {
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
SEXP gc_create_join_bands_cube(Rcpp::List pin_list,  std::vector<std::string> cube_names) {
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
SEXP gc_create_select_bands_cube(SEXP pin, std::vector<std::string> bands) {
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
SEXP gc_create_select_time_cube(SEXP pin, std::vector<std::string> t) {
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
SEXP gc_create_apply_pixel_cube(SEXP pin, std::vector<std::string> expr, std::vector<std::string> names, bool keep_bands=false) {
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
SEXP gc_create_stream_apply_pixel_cube(SEXP pin, std::string cmd, uint16_t nbands, std::vector<std::string> names, bool keep_bands = false) {
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
SEXP gc_create_stream_apply_time_cube(SEXP pin, std::string cmd, uint16_t nbands, std::vector<std::string> names, bool keep_bands = false) {
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
SEXP gc_create_filter_predicate_cube(SEXP pin, std::string pred) {
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
SEXP gc_create_filter_geom_cube(SEXP pin, std::string wkt, std::string srs) {
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
void gc_set_err_handler(bool debug, std::string log_to_file = "") {
  try {
    if (log_to_file.empty()) {
      if (debug) {
        config::instance()->set_error_handler(error_handling_r::debug); 
      }
      else {
        config::instance()->set_error_handler(error_handling_r::standard); 
      }
    }
    else {
      error_handling_r::_logfile = log_to_file;
      if (debug) {
        config::instance()->set_error_handler(error_handling_r::debug_file); 
      }
      else {
        config::instance()->set_error_handler(error_handling_r::standard_file); 
      }
    }
  }
  catch (std::string s) {
    Rcpp::stop(s);
  }
}




// [[Rcpp::export]]
void gc_eval_cube( SEXP pin, std::string outfile, uint8_t compression_level=0, bool with_VRT=false, 
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
void gc_write_chunks_ncdf( SEXP pin, std::string dir, std::string name, uint8_t compression_level=0) {
  try {
    Rcpp::XPtr< std::shared_ptr<cube> > aa = Rcpp::as<Rcpp::XPtr< std::shared_ptr<cube> >>(pin);
    (*aa)->write_chunks_netcdf(dir, name, compression_level);
  }
  catch (std::string s) {
    Rcpp::stop(s);
  }
}

// [[Rcpp::export]]
void gc_write_tif( SEXP pin, std::string dir, std::string prefix="", 
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
SEXP gc_create_stream_cube(SEXP pin, std::string cmd) {
  try {
    Rcpp::XPtr< std::shared_ptr<cube> > aa = Rcpp::as<Rcpp::XPtr< std::shared_ptr<cube> >>(pin);
    
    std::shared_ptr<stream_cube>* x = new std::shared_ptr<stream_cube>( stream_cube::create(*aa, cmd));
    
    Rcpp::XPtr< std::shared_ptr<stream_cube> > p(x, true) ;
  
    return p;
  }
  catch (std::string s) {
    Rcpp::stop(s);
  }
}



// [[Rcpp::export]]
SEXP gc_create_simple_cube(std::vector<std::string> files,std::vector<std::string> datetime_values, 
                                     std::vector<std::string> bands,  std::vector<std::string> band_names, 
                                     double dx, double dy,  Rcpp::IntegerVector chunk_sizes) {
  try {
    std::shared_ptr<simple_cube>* x = new std::shared_ptr<simple_cube>( simple_cube::create(files, datetime_values,
                                                                                            bands, band_names,
                                                                                            dx, dy));
    (*x)->set_chunk_size(chunk_sizes[0], chunk_sizes[1], chunk_sizes[2]);
    Rcpp::XPtr< std::shared_ptr<simple_cube> > p(x, true) ;
    return p;
  }
  catch (std::string s) {
    Rcpp::stop(s);
  }
}



// [[Rcpp::export]]
SEXP gc_create_fill_time_cube(SEXP pin, std::string method) {
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
SEXP gc_create_aggregate_time_cube(SEXP pin, std::string dt, std::string method, uint32_t fact=0) {
  try {
    Rcpp::XPtr< std::shared_ptr<cube> > aa = Rcpp::as<Rcpp::XPtr< std::shared_ptr<cube> >>(pin);
    
    std::shared_ptr<aggregate_time_cube>* x;
    if (fact >= 1) {
      x = new std::shared_ptr<aggregate_time_cube>(aggregate_time_cube::create(*aa, fact, method));
    }
    else {
      x = new std::shared_ptr<aggregate_time_cube>(aggregate_time_cube::create(*aa, dt, method));
    }
    Rcpp::XPtr< std::shared_ptr<aggregate_time_cube> > p(x, true) ;
    return p;
  } 
  catch (std::string s) {
    Rcpp::stop(s);
  }
}


// [[Rcpp::export]]
SEXP gc_create_aggregate_space_cube(SEXP pin, double dx, double dy, std::string method, uint32_t fact=0) {
  try {
    Rcpp::XPtr< std::shared_ptr<cube> > aa = Rcpp::as<Rcpp::XPtr< std::shared_ptr<cube> >>(pin);
    
    std::shared_ptr<aggregate_space_cube>* x;
    if (fact >= 1) {
      x = new std::shared_ptr<aggregate_space_cube>(aggregate_space_cube::create(*aa, fact, method));
    }
    else {
      x = new std::shared_ptr<aggregate_space_cube>(aggregate_space_cube::create(*aa, dx, dy, method));
    }
    Rcpp::XPtr< std::shared_ptr<aggregate_space_cube> > p(x, true) ;
    return p;
  } 
  catch (std::string s) {
    Rcpp::stop(s);
  }
}


// [[Rcpp::export]]
SEXP gc_create_slice_time_cube(SEXP pin, std::string dt, int32_t it=0) {
  try {
    Rcpp::XPtr< std::shared_ptr<cube> > aa = Rcpp::as<Rcpp::XPtr< std::shared_ptr<cube> >>(pin);
    
    std::shared_ptr<slice_time_cube>* x;
    if (dt.empty()) {
      x = new std::shared_ptr<slice_time_cube>(slice_time_cube::create(*aa, it));
    }
    else {
      x = new std::shared_ptr<slice_time_cube>(slice_time_cube::create(*aa, dt));
    }
    Rcpp::XPtr< std::shared_ptr<slice_time_cube> > p(x, true) ;
    return p;
  } 
  catch (std::string s) {
    Rcpp::stop(s);
  }
}


// [[Rcpp::export]]
SEXP gc_create_slice_space_cube(SEXP pin, std::vector<double> loc, std::vector<int32_t> i) {
  try {
    Rcpp::XPtr< std::shared_ptr<cube> > aa = Rcpp::as<Rcpp::XPtr< std::shared_ptr<cube> >>(pin);
    
    std::shared_ptr<slice_space_cube>* x;
    if (loc.empty()) {
      x = new std::shared_ptr<slice_space_cube>(slice_space_cube::create(*aa, i[0], i[1]));
    }
    else {
      x = new std::shared_ptr<slice_space_cube>(slice_space_cube::create(*aa, loc[0], loc[1]));
    }
    Rcpp::XPtr< std::shared_ptr<slice_space_cube> > p(x, true) ;
    return p;
  } 
  catch (std::string s) {
    Rcpp::stop(s);
  }
}



// [[Rcpp::export]]
SEXP gc_create_crop_cube(SEXP pin, Rcpp::List extent, std::vector<int32_t> iextent, std::string snap) {
  try {
    Rcpp::XPtr< std::shared_ptr<cube> > aa = Rcpp::as<Rcpp::XPtr< std::shared_ptr<cube> >>(pin);
    std::shared_ptr<crop_cube>* x;
    
    if (iextent.empty()) {
      // x = new std::shared_ptr<crop_cube>(crop_cube::create(*aa, Rcpp::as<double>(extent["left"]), 
      //                                                      Rcpp::as<double>(extent["right"]), Rcpp::as<double>(extent["bottom"]), 
      //                                                      Rcpp::as<double>(extent["top"]), Rcpp::as<std::string>(extent["t0"]),  
      //                                                      Rcpp::as<std::string>(extent["t1"]), Rcpp::as<std::string>(snap)));
      x = new std::shared_ptr<crop_cube>(crop_cube::create(*aa, (extent["left"]), 
                                                           (extent["right"]), (extent["bottom"]), 
                                                           (extent["top"]), (extent["t0"]),  
                                                           (extent["t1"]), (snap)));
    }
    else {
      x = new std::shared_ptr<crop_cube>(crop_cube::create(*aa, iextent[0], iextent[1], iextent[2], iextent[3], iextent[4], iextent[5]));
    }
    Rcpp::XPtr< std::shared_ptr<crop_cube> > p(x, true) ;
    return p;
  } 
  catch (std::string s) {
    Rcpp::stop(s);
  }
}


// [[Rcpp::export]]
Rcpp::DataFrame gc_extract(SEXP pin, std::string ogr_dataset, std::string time_column = "") {
  try {
    CPLPushErrorHandler(config::gdal_err_handler_default);
    Rcpp::XPtr< std::shared_ptr<cube> > aa = Rcpp::as<Rcpp::XPtr< std::shared_ptr<cube> >>(pin);
    auto x = std::shared_ptr<extract_geom>(extract_geom::create(*aa, ogr_dataset, time_column));
    
    std::vector<std::vector<double>> out;
    out.resize(x->size_bands());
    
    auto p = config::instance()->get_default_chunk_processor();
    
    std::shared_ptr<progress> prg = config::instance()->get_default_progress_bar()->get();
    prg->set(0);  // explicitly set to zero to show progress bar immediately
    std::function<void(chunkid_t, std::shared_ptr<chunk_data>, std::mutex &)> f = [prg, &out, x](chunkid_t id, std::shared_ptr<chunk_data> dat, std::mutex &m) {
      if (!dat->empty()) {
        for (uint32_t i=0; i<out.size(); ++i) {
          out[i].insert(out[i].end(), &((double*)(dat->buf()))[i*dat->size()[1]],&((double*)(dat->buf()))[(i+1)*dat->size()[1]]);
        }
      }
      prg->increment((double)1 / (double)x->count_chunks());
    };
    p->apply(x, f);
    prg->finalize();
    
    
    uint32_t ncol = out.size();
    uint32_t nrow = out[0].size();
    
    Rcpp::List df; 
    IntegerVector col_FID(nrow);
    for (uint32_t i=0; i<nrow; ++i) {
      col_FID[i] = (int)out[0][i];
    }
    df.push_back(col_FID);
    CharacterVector col_time(nrow);
    for (uint32_t i=0; i<nrow; ++i) {
      col_time[i] = (*aa)->st_reference()->datetime_at_index((int)out[1][i]).to_string();
    }
    df.push_back(col_time);
    
    for (uint32_t j=2; j<ncol; ++j) {
      df.push_back(out[j]);
    }
     

    StringVector row_names(nrow);
    for (uint32_t i = 0; i < nrow; ++i) {
      row_names(i) = std::to_string(i+1);
    }
    df.attr("row.names") = row_names;
    
    
    StringVector col_names(ncol);
    col_names(0) = "FID";
    col_names(1) = "time";
    for (uint32_t i = 2; i < ncol; ++i) {
      col_names(i) = (*aa)->bands().get(i-2).name;
    }
    df.attr("names") = col_names;
    df.attr("class") = "data.frame";
    
    CPLPopErrorHandler();
    return df;
  }
  catch (std::string s) {
    CPLPopErrorHandler();
    Rcpp::stop(s);
  } 
}








// [[Rcpp::export]]
void gc_exec_worker(std::string json_path, uint32_t pid, uint32_t nworker, std::string work_dir, int compression = 0) {
  chunk_processor_multiprocess::exec(json_path, pid, nworker, work_dir, compression);
}


// [[Rcpp::export]]
void gc_set_process_execution(IntegerVector n_worker, std::string cmd, bool debug, int ncdf_compression_level, 
                              bool use_overviews, Rcpp::List gdal_options) {
  auto p = std::make_shared<chunk_processor_multiprocess>();
  p->set_cmd(cmd);
  p->set_nworker(n_worker[0]);
  p->set_debug(debug);
  p->set_ncdf_compression_level(ncdf_compression_level);
  p->set_use_overviews(use_overviews);
  
  std::unordered_map<std::string, std::string> opt;
  if (gdal_options.size() > 0) {
    std::vector<std::string> nms = gdal_options.names();
    if ((size_t)gdal_options.size() == (size_t)nms.size()) {
      for (int i=0; i<gdal_options.size(); ++i) {
        std::string key = nms[i];
        std::string value = gdal_options[key];
        opt[key] = value;
      }
      p->set_gdal_options(opt);
    }
  }
  config::instance()->set_default_chunk_processor(std::dynamic_pointer_cast<chunk_processor>(p));
}

// [[Rcpp::export]]
void gc_set_progress(bool show_progress) {
  if (show_progress) {
    config::instance()->set_default_progress_bar(std::make_shared<progress_simple_R>());
  }
  else {
    config::instance()->set_default_progress_bar(std::make_shared<progress_none_R>());
  }
}


// [[Rcpp::export]]
void gc_set_use_overviews(bool use_overviews) {
  config::instance()->set_gdal_use_overviews(use_overviews);
}

// [[Rcpp::export]]
int gc_detect_cores() {
  return std::thread::hardware_concurrency();
}


// [[Rcpp::export]]
std::string gc_simple_hash(std::string instr) {
  return utils::hash(instr);
}

// [[Rcpp::export]]
void gc_create_stac_collection(Rcpp::DataFrame bands, Rcpp::DataFrame images, Rcpp::DataFrame gdalrefs, std::string outfile, Rcpp::DataFrame image_md) {
  
  try {
    //std::shared_ptr<image_collection>* x = new std::shared_ptr<image_collection>();
    std::shared_ptr<image_collection> x = image_collection::create();
    
    x->transaction_start();
    
    Rcpp::CharacterVector band_name =  bands["name"];
    Rcpp::IntegerVector band_id = bands["id"];
    for (int32_t i=0; i<bands.nrows(); ++i) {
      x->insert_band(band_id[i], Rcpp::as<std::string>(band_name[i])); // TODO add further data if available
    }

    Rcpp::IntegerVector image_id = images["id"];
    Rcpp::CharacterVector image_name = images["name"];
    Rcpp::NumericVector image_left = images["left"];
    Rcpp::NumericVector image_top = images["top"];
    Rcpp::NumericVector image_bottom =images["bottom"];
    Rcpp::NumericVector image_right =images["right"];
    Rcpp::CharacterVector image_datetime = images["datetime"];
    Rcpp::CharacterVector image_proj = images["proj"];
    for (int32_t i=0; i<images.nrows(); ++i) {
      x->insert_image(image_id[i], Rcpp::as<std::string>(image_name[i]), image_left[i], image_top[i],
                         image_bottom[i], image_right[i], Rcpp::as<std::string>(image_datetime[i]), Rcpp::as<std::string>(image_proj[i]));
    }

    Rcpp::IntegerVector gdalrefs_image_id = gdalrefs["image_id"];
    Rcpp::IntegerVector gdalrefs_band_id = gdalrefs["band_id"];
    Rcpp::CharacterVector gdalrefs_descriptor = gdalrefs["descriptor"];
    Rcpp::IntegerVector gdalrefs_band_num = gdalrefs["band_num"];
    for (int32_t i=0; i<gdalrefs.nrows(); ++i) {
      x->insert_dataset(gdalrefs_image_id[i],gdalrefs_band_id[i], Rcpp::as<std::string>(gdalrefs_descriptor[i]), gdalrefs_band_num[i]);
    }
    
    
    if (image_md.nrows() > 0) {
      for (int32_t i=0; i<image_md.nrows(); ++i) {
        Rcpp::IntegerVector image_md_image_id = image_md["image_id"];
        Rcpp::CharacterVector image_md_key = image_md["key"];
        Rcpp::CharacterVector image_md_value = image_md["value"];
        
        x->insert_image_md(image_md_image_id[i], Rcpp::as<std::string>(image_md_key[i]), Rcpp::as<std::string>(image_md_value[i]));
      }
    }
    
    // TODO: add collection, image, and band metadata
    x->transaction_end();
    
    x->write(outfile);
  }
  catch (std::string s) {
    Rcpp::stop(s);
  }
}

