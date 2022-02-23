#ifndef GC_ERROR_H
#define GC_ERROR_H

#include <mutex>
#include <string>
#include <sstream>
#include "gdalcubes/src/error.h"

struct r_stderr_buf {
  static std::mutex _m;
  static std::stringstream _s;
  static void print(std::string s="");
};

using namespace gdalcubes;

struct error_handling_r {
  static std::mutex _m_errhandl;
  static std::stringstream _err_stream;
  static bool _defer;
  static std::string _logfile;
  
  static void defer_output();
  static void do_output(); 
  
  static void debug(error_level type, std::string msg, std::string where, int error_code);
  static void standard(error_level type, std::string msg, std::string where, int error_code);
    
  static void standard_file(error_level type, std::string msg, std::string where, int error_code);
  static void debug_file(error_level type, std::string msg, std::string where, int error_code);
};

#endif