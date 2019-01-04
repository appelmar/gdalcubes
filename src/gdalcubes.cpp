
#include <gdalcubes.h>

// [[Rcpp::plugins("cpp11")]]
// [[Rcpp::depends(RcppProgress)]]
#include <Rcpp.h>
#include <progress.hpp>
#include <progress_bar.hpp>
#include <memory>

using namespace Rcpp;


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
    Rcpp::Rcout << _err_stream.str() << std::endl;
    _err_stream.str(""); 
    _m_errhandl.unlock();
  }
  
  static void debug(error_level type, std::string msg, std::string where, int error_code) {
    _m_errhandl.lock();
    std::string code = (error_code != 0) ? " (" + std::to_string(error_code) + ")" : "";
    std::string where_str = (where.empty()) ? "" : " [in " + where + "]";
    if (type == error_level::ERRLVL_ERROR || type == error_level::ERRLVL_FATAL ) {
      _err_stream << "ERROR" << code << ": " << msg << where_str << std::endl;
    } else if (type == error_level::ERRLVL_WARNING) {
      _err_stream << "WARNING" << code << ": " << msg << where_str << std::endl;
    } else if (type == error_level::ERRLVL_INFO) {
      _err_stream << "INFO" << code << ": " << msg << where_str << std::endl;
    } else if (type == error_level::ERRLVL_DEBUG) {
      _err_stream << "DEBUG" << code << ": " << msg << where_str << std::endl;
    }
    if (!_defer) {
      Rcpp::Rcout << _err_stream.str() << std::endl;
      _err_stream.str(""); 
    }
    _m_errhandl.unlock();
  }
  
  static void standard(error_level type, std::string msg, std::string where, int error_code) {
    _m_errhandl.lock();
    std::string code = (error_code != 0) ? " (" + std::to_string(error_code) + ")" : "";
    std::string where_str = (where.empty()) ? "" : " [in " + where + "]";
    if (type == error_level::ERRLVL_ERROR || type == error_level::ERRLVL_FATAL) {
      _err_stream << "ERROR" << code << ": " << msg << where_str << std::endl;
    } else if (type == error_level::ERRLVL_WARNING) {
      _err_stream << "WARNING" << code << ": " << msg << where_str << std::endl;
    } else if (type == error_level::ERRLVL_INFO) {
      _err_stream << "INFO" << code << ": " << msg << where_str << std::endl;
    }
    if (!_defer) {
      Rcpp::Rcout << _err_stream.str() << std::endl;
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
  
  void _set(double p) {
    //Rcpp::checkUserInterrupt();
    // if (Progress::check_abort()) {
    //   throw std::string("Operation has been interrupted by user");
    // }
    if (!_rp) {
      error_handling_r::defer_output();
      _rp = new Progress(100,true);
    }
    double p_old = _p;
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
void libgdalcubes_init() {
  config::instance()->gdalcubes_init();
  config::instance()->set_default_progress_bar(std::make_shared<progress_simple_R>());
  //config::instance()->set_default_progress_bar(std::make_shared<progress_none>());
  
  config::instance()->set_error_handler(error_handling_r::standard); 
}

// [[Rcpp::export]]
void libgdalcubes_cleanup() {
  config::instance()->gdalcubes_cleanup();
}

// [[Rcpp::export]]
Rcpp::StringVector libgdalcubes_datetime_values(SEXP pin) {
  Rcpp::XPtr<std::shared_ptr<cube>> aa = Rcpp::as<Rcpp::XPtr<std::shared_ptr<cube>>>(pin);
  std::shared_ptr<cube> x = *aa;
  
  Rcpp::CharacterVector out(x->size_t());
  
  for (uint32_t i = 0; i < x->size_t(); ++i) {
    out[i] = (x->st_reference()->t0() + ( x->st_reference()->dt() * i)).to_string();
  }
  return out;
}

// [[Rcpp::export]]
Rcpp::List libgdalcubes_cube_info( SEXP pin) {

  Rcpp::XPtr<std::shared_ptr<cube>> aa = Rcpp::as<Rcpp::XPtr<std::shared_ptr<cube>>>(pin);
  
  std::shared_ptr<cube> x = *aa;
  
  Rcpp::CharacterVector d_name(3);
  Rcpp::NumericVector d_low(3);
  Rcpp::NumericVector d_high(3);
  Rcpp::IntegerVector d_n(3);
  Rcpp::IntegerVector d_chunk(3);

  d_name[0] = "t";
  d_low[0] = x->st_reference()->t0().to_double();
  d_high[0] = x->st_reference()->t1().to_double();
  d_n[0] = x->st_reference()->nt();
  d_chunk[0] = x->chunk_size()[0];

  d_name[1] = "y";
  d_low[1] = x->st_reference()->bottom();
  d_high[1] = x->st_reference()->top();
  d_n[1] = x->st_reference()->ny();
  d_chunk[1] = x->chunk_size()[1];

  d_name[2] = "x";
  d_low[2] = x->st_reference()->left();
  d_high[2] = x->st_reference()->right();
  d_n[2] = x->st_reference()->nx();
  d_chunk[2] = x->chunk_size()[2];

  Rcpp::DataFrame dims =
    Rcpp::DataFrame::create(Rcpp::Named("name")=d_name,
                            Rcpp::Named("low")=d_low,
                            Rcpp::Named("high")=d_high,
                            Rcpp::Named("size")=d_n,
                            Rcpp::Named("chunk_size")=d_chunk);

  Rcpp::CharacterVector b_names(x->bands().count(), "");
  Rcpp::CharacterVector b_type(x->bands().count(), "");
  Rcpp::NumericVector b_offset(x->bands().count(), NA_REAL);
  Rcpp::NumericVector b_scale(x->bands().count(), NA_REAL);
  Rcpp::NumericVector b_nodata(x->bands().count(), NA_REAL);
  Rcpp::CharacterVector b_unit(x->bands().count(), "");

  for (uint16_t i=0; i<x->bands().count(); ++i) {
    b_names[i] = x->bands().get(i).name;
    b_type[i] = x->bands().get(i).type;
    b_offset[i] = x->bands().get(i).offset;
    b_scale[i] = x->bands().get(i).scale;
    b_nodata[i] = (x->bands().get(i).no_data_value.empty())? NA_REAL : std::stod(x->bands().get(i).no_data_value);
    b_unit[i] = x->bands().get(i).unit;
  }

  Rcpp::DataFrame bands =
    Rcpp::DataFrame::create(Rcpp::Named("name")=b_names,
                           // Rcpp::Named("type")=b_type,
                            Rcpp::Named("offset")=b_offset,
                            Rcpp::Named("scale")=b_scale,
                            Rcpp::Named("nodata")=b_nodata,
                            Rcpp::Named("unit")=b_unit);

  return Rcpp::List::create(Rcpp::Named("bands") = bands,
                            Rcpp::Named("dimensions") = dims,
                            Rcpp::Named("proj") = x->st_reference()->proj(),
                            Rcpp::Named("graph") = x->make_constructible_json().dump(2),
                            Rcpp::Named("size") = Rcpp::IntegerVector::create(x->size()[0], x->size()[1], x->size()[2], x->size()[3]));

}

// [[Rcpp::export]]
Rcpp::List libgdalcubes_get_cube_view( SEXP pin) {
  
  Rcpp::XPtr<std::shared_ptr<cube>> aa = Rcpp::as<Rcpp::XPtr<std::shared_ptr<cube>>>(pin);
  
  std::shared_ptr<cube> x = *aa;
  
  std::shared_ptr<cube_view> v = std::dynamic_pointer_cast<cube_view>(x->st_reference());
  
  Rcpp::List view = Rcpp::List::create(
    Rcpp::Named("space") = Rcpp::List::create(
      Rcpp::Named("right") = x->st_reference()->right(),
      Rcpp::Named("left") = x->st_reference()->left(),
      Rcpp::Named("top") = x->st_reference()->top(),
      Rcpp::Named("bottom") = x->st_reference()->bottom(),
      Rcpp::Named("nx") = x->st_reference()->nx(),
      Rcpp::Named("ny") = x->st_reference()->ny(),
      Rcpp::Named("proj") = x->st_reference()->proj(),
      Rcpp::Named("dx") = R_NilValue,
      Rcpp::Named("dy") = R_NilValue
    ),
    Rcpp::Named("time") = Rcpp::List::create(
      Rcpp::Named("t0") = x->st_reference()->t0().to_string(),
      Rcpp::Named("t1") = x->st_reference()->t1().to_string(),
      Rcpp::Named("dt") = x->st_reference()->dt().to_string(),
      Rcpp::Named("nt") = R_NilValue
    ),
    Rcpp::Named("aggregation") = v ? aggregation::to_string(v->aggregation_method()) : "none",
    Rcpp::Named("resampling") = v ? resampling::to_string(v->resampling_method()) : "near"
  );
  return view;
}
 

// TODO: do not update the view of pin directly
// but update view for all image_collection_cubes that are somehow connected....
// [[Rcpp::export]]
void libgdalcubes_update_cube_view( SEXP pin, SEXP v) {
  
  Rcpp::XPtr<std::shared_ptr<cube>> aa = Rcpp::as<Rcpp::XPtr<std::shared_ptr<cube>>>(pin);
  
  std::shared_ptr<cube> x = *aa;
  
  std::shared_ptr<cube_view> vv = std::dynamic_pointer_cast<cube_view>(x->st_reference());
  
  // 1. Create cube view shared_ptr that is default from pin and update in place...
  std::shared_ptr<cube_view> s = std::make_shared<cube_view>();
  s->win() = x->st_reference()->win(); 
  s->proj() = x->st_reference()->proj();
  s->nx() = x->st_reference()->nx();
  s->ny() = x->st_reference()->ny();
  s->t0() =  x->st_reference()->t0();
  s->t1() =  x->st_reference()->t1();
  s->dt() =  x->st_reference()->dt();
  
  Rcpp::List view = Rcpp::as<Rcpp::List>(v);
  if (Rcpp::as<Rcpp::List>(view["space"])["right"] != R_NilValue) {
    s->right() = Rcpp::as<Rcpp::List>(view["space"])["right"];
  }
  if (Rcpp::as<Rcpp::List>(view["space"])["left"] != R_NilValue) {
    s->left() = Rcpp::as<Rcpp::List>(view["space"])["left"];
  }
  if (Rcpp::as<Rcpp::List>(view["space"])["top"] != R_NilValue) {
    s->top() = Rcpp::as<Rcpp::List>(view["space"])["top"];
  }
  if (Rcpp::as<Rcpp::List>(view["space"])["bottom"] != R_NilValue) {
   s->bottom() = Rcpp::as<Rcpp::List>(view["space"])["bottom"];
  }
  // nx overwrites dx
  if (Rcpp::as<Rcpp::List>(view["space"])["dx"] != R_NilValue) {
   s->dx(Rcpp::as<Rcpp::List>(view["space"])["dx"]);
  }
  if (Rcpp::as<Rcpp::List>(view["space"])["nx"] != R_NilValue) {
    s->nx() = Rcpp::as<Rcpp::List>(view["space"])["nx"];
  }
  // ny overwrites dy
  if (Rcpp::as<Rcpp::List>(view["space"])["dy"] != R_NilValue) {
   s->dy(Rcpp::as<Rcpp::List>(view["space"])["dy"]);
  }
  if (Rcpp::as<Rcpp::List>(view["space"])["ny"] != R_NilValue) {
    s->ny() = Rcpp::as<Rcpp::List>(view["space"])["ny"];
  }
  if (Rcpp::as<Rcpp::List>(view["space"])["proj"] != R_NilValue) {
    s->proj() = Rcpp::as<Rcpp::CharacterVector>(Rcpp::as<Rcpp::List>(view["space"])["proj"])[0];
  }
  if (Rcpp::as<Rcpp::List>(view["time"])["t0"] != R_NilValue) {
    std::string tmp = Rcpp::as<Rcpp::String>(Rcpp::as<Rcpp::List>(view["time"])["t0"]);
    s->t0() = datetime::from_string(tmp);
  }
  if (Rcpp::as<Rcpp::List>(view["time"])["t1"] != R_NilValue) {
    std::string tmp = Rcpp::as<Rcpp::String>(Rcpp::as<Rcpp::List>(view["time"])["t1"]);
    s->t1() = datetime::from_string(tmp);
  }
  
  // dt overwrites nt
  if (Rcpp::as<Rcpp::List>(view["time"])["nt"] != R_NilValue) {
    s->nt(Rcpp::as<Rcpp::List>(view["time"])["nt"]);
  }
  if (Rcpp::as<Rcpp::List>(view["time"])["dt"] != R_NilValue) {
    std::string tmp = Rcpp::as<Rcpp::String>(Rcpp::as<Rcpp::List>(view["time"])["dt"]);
    s->dt() = duration::from_string(tmp);
    s->t0().unit() = x->st_reference()->dt().dt_unit; 
    s->t1().unit() = x->st_reference()->dt().dt_unit; 
  }
  
  if (vv) {
    if (view["aggregation"] != R_NilValue) {
      std::string tmp = Rcpp::as<Rcpp::String>(view["aggregation"]);
     s->aggregation_method() = aggregation::from_string(tmp);
    }
    if (view["resampling"] != R_NilValue) {
      std::string tmp = Rcpp::as<Rcpp::String>(view["resampling"]);
      s->resampling_method() = resampling::from_string(tmp);
    }
  }
  
  // 2. Call x->update_st_reference()...
  x->update_st_reference(s);
  
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
                              Rcpp::Named("proj")=images_proj);
    
    
    
    
    std::vector<image_collection::bands_row> bands = ic->get_bands();
    
    Rcpp::CharacterVector bands_name(bands.size());
    Rcpp::CharacterVector bands_type(bands.size());
    Rcpp::NumericVector bands_offset(bands.size());
    Rcpp::NumericVector bands_scale(bands.size());
    Rcpp::CharacterVector bands_unit(bands.size());
    Rcpp::CharacterVector bands_nodata(bands.size());
    
    for (uint32_t i=0; i<bands.size(); ++i) {
      bands_name[i] = bands[i].name;
      bands_type[i] = bands[i].type;
      bands_offset[i] = bands[i].offset;
      bands_scale[i] = bands[i].scale;
      bands_unit[i] = bands[i].unit;
      bands_nodata[i] = bands[i].nodata;
    }
    
    Rcpp::DataFrame bands_df =
      Rcpp::DataFrame::create(Rcpp::Named("name")=bands_name,
                              Rcpp::Named("type")=bands_type,
                              Rcpp::Named("offset")=bands_offset,
                              Rcpp::Named("scale")=bands_scale,
                              Rcpp::Named("unit")=bands_unit,
                              Rcpp::Named("nodata")=bands_nodata);
    
    
    
    
    
    
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
SEXP libgdalcubes_create_image_collection(std::vector<std::string> files, std::string format_file, std::string outfile, bool unroll_archives=true) {

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
SEXP libgdalcubes_list_collection_formats() {
  try {
  
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
SEXP libgdalcubes_create_image_collection_cube(SEXP pin, Rcpp::IntegerVector chunk_sizes, SEXP v = R_NilValue) {

  try {
    
    Rcpp::XPtr<std::shared_ptr<image_collection>> aa = Rcpp::as<Rcpp::XPtr<std::shared_ptr<image_collection>>>(pin);
    
   
    
    std::shared_ptr<image_collection_cube>* x = new std::shared_ptr<image_collection_cube>( image_collection_cube::create(*aa));

    if (v != R_NilValue) {
      Rcpp::List view = Rcpp::as<Rcpp::List>(v);
      if (Rcpp::as<Rcpp::List>(view["space"])["right"] != R_NilValue) {
        (*x)->st_reference()->right() = Rcpp::as<Rcpp::List>(view["space"])["right"];
      }
      if (Rcpp::as<Rcpp::List>(view["space"])["left"] != R_NilValue) {
        (*x)->st_reference()->left() = Rcpp::as<Rcpp::List>(view["space"])["left"];
      }
      if (Rcpp::as<Rcpp::List>(view["space"])["top"] != R_NilValue) {
        (*x)->st_reference()->top() = Rcpp::as<Rcpp::List>(view["space"])["top"];
      }
      if (Rcpp::as<Rcpp::List>(view["space"])["bottom"] != R_NilValue) {
        (*x)->st_reference()->bottom() = Rcpp::as<Rcpp::List>(view["space"])["bottom"];
      }
      // nx overwrites dx
      if (Rcpp::as<Rcpp::List>(view["space"])["dx"] != R_NilValue) {
        (*x)->st_reference()->dx(Rcpp::as<Rcpp::List>(view["space"])["dx"]);
      }
      if (Rcpp::as<Rcpp::List>(view["space"])["nx"] != R_NilValue) {
        (*x)->st_reference()->nx() = Rcpp::as<Rcpp::List>(view["space"])["nx"];
      }
      // ny overwrites dy
      if (Rcpp::as<Rcpp::List>(view["space"])["dy"] != R_NilValue) {
        (*x)->st_reference()->dy(Rcpp::as<Rcpp::List>(view["space"])["dy"]);
      }
      if (Rcpp::as<Rcpp::List>(view["space"])["ny"] != R_NilValue) {
        (*x)->st_reference()->ny() = Rcpp::as<Rcpp::List>(view["space"])["ny"];
      }
      if (Rcpp::as<Rcpp::List>(view["space"])["proj"] != R_NilValue) {
        (*x)->st_reference()->proj() = Rcpp::as<Rcpp::CharacterVector>(Rcpp::as<Rcpp::List>(view["space"])["proj"])[0];
      }
      if (Rcpp::as<Rcpp::List>(view["time"])["t0"] != R_NilValue) {
        std::string tmp = Rcpp::as<Rcpp::String>(Rcpp::as<Rcpp::List>(view["time"])["t0"]);
        (*x)->st_reference()->t0() = datetime::from_string(tmp);
      }
      if (Rcpp::as<Rcpp::List>(view["time"])["t1"] != R_NilValue) {
        std::string tmp = Rcpp::as<Rcpp::String>(Rcpp::as<Rcpp::List>(view["time"])["t1"]);
        (*x)->st_reference()->t1() = datetime::from_string(tmp);
      }
      
      // dt overwrites nt
      if (Rcpp::as<Rcpp::List>(view["time"])["nt"] != R_NilValue) {
        (*x)->st_reference()->nt(Rcpp::as<Rcpp::List>(view["time"])["nt"]);
      }
      if (Rcpp::as<Rcpp::List>(view["time"])["dt"] != R_NilValue) {
        std::string tmp = Rcpp::as<Rcpp::String>(Rcpp::as<Rcpp::List>(view["time"])["dt"]);
        (*x)->st_reference()->dt() = duration::from_string(tmp);
        (*x)->st_reference()->t0().unit() = (*x)->st_reference()->dt().dt_unit; 
        (*x)->st_reference()->t1().unit() = (*x)->st_reference()->dt().dt_unit; 
      }
      
      
      if (view["aggregation"] != R_NilValue) {
        std::string tmp = Rcpp::as<Rcpp::String>(view["aggregation"]);
        std::dynamic_pointer_cast<cube_view>((*x)->st_reference())->aggregation_method() = aggregation::from_string(tmp);
      }
      if (view["resampling"] != R_NilValue) {
        std::string tmp = Rcpp::as<Rcpp::String>(view["resampling"]);
        std::dynamic_pointer_cast<cube_view>((*x)->st_reference())->resampling_method() = resampling::from_string(tmp);
      }
    }
    (*x)->set_chunk_size(chunk_sizes[0], chunk_sizes[1], chunk_sizes[2]);
    
    
    //Rcpp::Rcout << std::dynamic_pointer_cast<cube_view>((*x)->st_reference())->write_json_string() << std::endl;
    
    Rcpp::XPtr< std::shared_ptr<image_collection_cube> > p(x, true) ;
    
    return p;
    
  }
  catch (std::string s) {
    Rcpp::stop(s);
  }
 
  
}


// [[Rcpp::export]]
SEXP libgdalcubes_create_reduce_cube(SEXP pin, std::string reducer) {
  try {
    Rcpp::XPtr< std::shared_ptr<cube> > aa = Rcpp::as<Rcpp::XPtr<std::shared_ptr<cube>>>(pin);
    
    std::shared_ptr<reduce_cube>* x = new std::shared_ptr<reduce_cube>(reduce_cube::create(*aa, reducer));
    Rcpp::XPtr< std::shared_ptr<reduce_cube> > p(x, true) ;
    
    return p;
    
  }
  catch (std::string s) {
    Rcpp::stop(s);
  }
}



// [[Rcpp::export]]
SEXP libgdalcubes_create_join_bands_cube(SEXP pinA, SEXP pinB) {
  try {
    Rcpp::XPtr< std::shared_ptr<cube> > A = Rcpp::as<Rcpp::XPtr<std::shared_ptr<cube>>>(pinA);
    Rcpp::XPtr< std::shared_ptr<cube> > B = Rcpp::as<Rcpp::XPtr<std::shared_ptr<cube>>>(pinB);
    
    std::shared_ptr<join_bands_cube>* x = new std::shared_ptr<join_bands_cube>(join_bands_cube::create(*A, *B));
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
SEXP libgdalcubes_create_apply_pixel_cube(SEXP pin, std::vector<std::string> expr, std::vector<std::string> names = std::vector<std::string>()) {
  try {
    Rcpp::XPtr< std::shared_ptr<cube> > aa = Rcpp::as<Rcpp::XPtr<std::shared_ptr<cube>>>(pin);
    
    std::shared_ptr<apply_pixel_cube>* x = new std::shared_ptr<apply_pixel_cube>(apply_pixel_cube::create(*aa, expr, names));
    Rcpp::XPtr< std::shared_ptr<apply_pixel_cube> > p(x, true) ;
    
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
void libgdalcubes_eval_cube( SEXP pin, std::string outfile) {
  try {
    Rcpp::XPtr< std::shared_ptr<cube> > aa = Rcpp::as<Rcpp::XPtr< std::shared_ptr<cube> >>(pin);
    (*aa)->write_netcdf_file(outfile);
  }
  catch (std::string s) {
    Rcpp::stop(s);
  }
}


// [[Rcpp::export]]
SEXP libgdalcubes_create_stream_cube(SEXP pin, std::string cmd) {
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
void libgdalcubes_set_threads(IntegerVector n) {
  if (n[0] > 1) {
    config::instance()->set_default_chunk_processor(std::dynamic_pointer_cast<chunk_processor>(std::make_shared<chunk_processor_multithread>(n[0])));
  }
  else {
    config::instance()->set_default_chunk_processor(std::dynamic_pointer_cast<chunk_processor>(std::make_shared<chunk_processor_singlethread>()));
  }
}


