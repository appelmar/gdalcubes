
#include "error.h"

#include <thread>
#include <fstream>

// [[Rcpp::plugins("cpp11")]]
#include <Rcpp.h>

static std::thread::id r_main_thread_id = std::this_thread::get_id();

std::mutex r_stderr_buf::_m;
std::stringstream r_stderr_buf::_s;

void r_stderr_buf::print(std::string s) {
  std::lock_guard<std::mutex> lck(_m);
  _s << s;
  if (!_s.str().empty() && r_main_thread_id == std::this_thread::get_id()) {
    REprintf("%s", _s.str().c_str());
    //R_FlushConsole(); // TODO: needed for stderr?
    _s.str("");
  }
}


std::mutex error_handling_r::_m_errhandl;
std::stringstream error_handling_r::_err_stream;
bool error_handling_r::_defer = false;
std::string error_handling_r::_logfile = "gdalcubes.log";



void error_handling_r::defer_output() {
  _m_errhandl.lock();
  _defer = true;
  _m_errhandl.unlock();
}
  
void error_handling_r::do_output() {
  _m_errhandl.lock();
  // TODO: if _err_stream is extremely large,
  // print to a file first and afterwards show only last few lines and 
  // a warning pointing to the file for full output
  if (_err_stream.rdbuf()->in_avail() > 0) {
    r_stderr_buf::print(_err_stream.str());
    _err_stream.str(""); 
  }
  _defer = false;
  _m_errhandl.unlock();
}

void error_handling_r::debug(error_level type, std::string msg, std::string where, int error_code) {
  _m_errhandl.lock();
  std::string code = (error_code != 0) ? " (" + std::to_string(error_code) + ")" : "";
  std::string where_str = (where.empty()) ? "" : " [in " + where + "]";
  if (type == error_level::ERRLVL_ERROR || type == error_level::ERRLVL_FATAL ) {
    _err_stream << "[ERROR] "  << msg << where_str << std::endl;
  } else if (type == error_level::ERRLVL_WARNING) {
    _err_stream << "[WARNING]  " << msg << where_str << std::endl;
  } else if (type == error_level::ERRLVL_INFO) {
    _err_stream << "[INFO] " << msg << where_str << std::endl;
  } else if (type == error_level::ERRLVL_DEBUG) {
    _err_stream << "[DEBUG] "  << msg << where_str << std::endl;
  }
  if (!_defer) {
    if (_err_stream.rdbuf()->in_avail() > 0) {
      r_stderr_buf::print(_err_stream.str());
      _err_stream.str(""); 
    }
  }
  _m_errhandl.unlock();
}

void error_handling_r::standard(error_level type, std::string msg, std::string where, int error_code) {
  _m_errhandl.lock();
  std::string code = (error_code != 0) ? " (" + std::to_string(error_code) + ")" : "";
  if (type == error_level::ERRLVL_ERROR || type == error_level::ERRLVL_FATAL) {
    _err_stream << "[ERROR] " << msg << std::endl;
  } else if (type == error_level::ERRLVL_WARNING) {
    _err_stream << "[WARNING] " << msg << std::endl;
  } else if (type == error_level::ERRLVL_INFO) {
    _err_stream << "## " << msg << std::endl;
  }
  if (!_defer) {
    if (_err_stream.rdbuf()->in_avail() > 0) {
      r_stderr_buf::print(_err_stream.str());
      _err_stream.str(""); 
    }
  }
  _m_errhandl.unlock();
}

void error_handling_r::standard_file(error_level type, std::string msg, std::string where, int error_code) {
  _m_errhandl.lock();
  std::ofstream os;
  os.open(_logfile, std::ios::out | std::ios::app);
  if (!os.is_open()) {
    _m_errhandl.unlock();
    standard(type, msg,  where, error_code);
    return;
  }
  std::string code = (error_code != 0) ? " (" + std::to_string(error_code) + ")" : "";
  if (type == error_level::ERRLVL_ERROR || type == error_level::ERRLVL_FATAL) {
    os << "[ERROR] " << msg << std::endl;
  } else if (type == error_level::ERRLVL_WARNING) {
    os << "[WARNING] " << msg << std::endl;
  } else if (type == error_level::ERRLVL_INFO) {
    os << "## " << msg << std::endl;
  }
  _m_errhandl.unlock();
}

void error_handling_r::debug_file(error_level type, std::string msg, std::string where, int error_code) {
  _m_errhandl.lock();
  std::ofstream os;
  os.open(_logfile, std::ios::out | std::ios::app);
  if (!os.is_open()) {
    _m_errhandl.unlock();
    debug(type, msg,  where, error_code);
    return;
  }
  std::string code = (error_code != 0) ? " (" + std::to_string(error_code) + ")" : "";
  std::string where_str = (where.empty()) ? "" : " [in " + where + "]";
  if (type == error_level::ERRLVL_ERROR || type == error_level::ERRLVL_FATAL ) {
    os << "[ERROR] "  << msg << where_str << std::endl;
  } else if (type == error_level::ERRLVL_WARNING) {
    os << "[WARNING] " << msg << where_str << std::endl;
  } else if (type == error_level::ERRLVL_INFO) {
    os << "[INFO] " << msg << where_str << std::endl;
  } else if (type == error_level::ERRLVL_DEBUG) {
    os << "[DEBUG] "  << msg << where_str << std::endl;
  }
  _m_errhandl.unlock();
}
