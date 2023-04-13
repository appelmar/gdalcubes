/*
    MIT License

    Copyright (c) 2019 Marius Appel <marius.appel@uni-muenster.de>

    Permission is hereby granted, free of charge, to any person obtaining a copy
    of this software and associated documentation files (the "Software"), to deal
    in the Software without restriction, including without limitation the rights
    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    copies of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in all
    copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
    SOFTWARE.
*/

#ifndef ERROR_H
#define ERROR_H

#include <iostream>
#include <mutex>
#include <string>

#include "utils.h"

namespace gdalcubes {

#define GCBS_FATAL(MSG) logger::fatal(MSG, std::string(__FILE__) + ":" + std::string(__func__) + ":" + std::to_string(__LINE__) + "")
#define GCBS_ERROR(MSG) logger::error(MSG, std::string(__FILE__) + ":" + std::string(__func__) + ":" + std::to_string(__LINE__) + "")
#define GCBS_WARN(MSG) logger::warn(MSG, std::string(__FILE__) + ":" + std::string(__func__) + ":" + std::to_string(__LINE__) + "")
#define GCBS_DEBUG(MSG) logger::debug(MSG, std::string(__FILE__) + ":" + std::string(__func__) + ":" + std::to_string(__LINE__) + "")
#define GCBS_INFO(MSG) logger::info(MSG, std::string(__FILE__) + ":" + std::string(__func__) + ":" + std::to_string(__LINE__) + "")
#define GCBS_TRACE(MSG) logger::trace(MSG, std::string(__FILE__) + ":" + std::string(__func__) + ":" + std::to_string(__LINE__) + "")

enum class error_level {
    ERRLVL_TRACE = 6,
    ERRLVL_DEBUG = 5,
    ERRLVL_INFO = 4,
    ERRLVL_WARNING = 3,
    ERRLVL_ERROR = 2,
    ERRLVL_FATAL = 1
};
/**
 * @brief Function pointer prototype for custom error handlers
 */
typedef void (*error_action)(error_level, std::string, std::string, int);

/**
 * @brief Error reporting and logging
 */
class logger {
   public:
    /**
     * @brief Report or log error messages of different types according
     * @param msg message string
     * @param where optional message source
     * @param error_code optional error code number
     * @see config::set_error_handler for customizing log ouptut
     */
    static void error(std::string msg, std::string where = "", int error_code = 0);

    /**
     * @copydoc error
     */
    static void warn(std::string msg, std::string where = "", int error_code = 0);

    /**
    * @copydoc error
    */
    static void debug(std::string msg, std::string where = "", int error_code = 0);

    /**
    * @copydoc error
    */
    static void fatal(std::string msg, std::string where = "", int error_code = 0);

    /**
    * @copydoc error
    */
    static void info(std::string msg, std::string where = "", int error_code = 0);

    /**
    * @copydoc error
    */
    static void trace(std::string msg, std::string where = "", int error_code = 0);

   private:
    static std::mutex _m;
};

/**
 * @brief A few implementations of static error handler functions
 */
class error_handler {
   public:
    /**
     * @brief Default error handler that prints error and warning messages to stderr and stdout respectively
     * @param type error level
     * @param msg message string
     * @param where location where the error / message comes from
     * @param error_code integer error code
     */
    static void default_error_handler(error_level type, std::string msg, std::string where, int error_code) {
#ifndef R_PACKAGE  // avoid std::cerr and std::cout for R package, this is only to reduce R CMD check warnings
        std::string code = (error_code != 0) ? " (" + std::to_string(error_code) + ")" : "";
        std::string where_str = (where.empty()) ? "" : " [in " + where + "]";
        if (type == error_level::ERRLVL_ERROR || type == error_level::ERRLVL_FATAL) {
            std::cerr << "ERROR" << code << ": " << msg << where_str << std::endl;
        } else if (type == error_level::ERRLVL_WARNING) {
            std::cout << "WARNING" << code << ": " << msg << where_str << std::endl;
        }
#endif
    }

    /**
     * @brief Default error handler for debugging that prints messages up to the debug level to stderr and stdout
     * @param type error level
     * @param msg message string
     * @param where location where the error / message comes from
     * @param error_code integer error code
     */
    static void error_handler_debug(error_level type, std::string msg, std::string where, int error_code) {
#ifndef R_PACKAGE  // avoid std::cerr and std::cout for R package, this is only to reduce R CMD check warnings
        std::string code = (error_code != 0) ? " (" + std::to_string(error_code) + ")" : "";
        std::string where_str = (where.empty()) ? "" : " [in " + where + "]";
        if (type == error_level::ERRLVL_ERROR || type == error_level::ERRLVL_FATAL) {
            std::cerr << "ERROR" << code << ": " << msg << where_str << std::endl;
        } else if (type == error_level::ERRLVL_WARNING) {
            std::cout << "WARNING" << code << ": " << msg << where_str << std::endl;
        } else if (type == error_level::ERRLVL_INFO) {
            std::cout << "INFO" << code << ": " << msg << where_str << std::endl;
        } else if (type == error_level::ERRLVL_DEBUG) {
            std::cout << "DEBUG" << code << ": " << msg << where_str << std::endl;
        }
#endif
    }

    /**
     * @brief Error handler for debugging gdalcubes_server including system time in messages
     * @param type error level
     * @param msg message string
     * @param where location where the error / message comes from
     * @param error_code integer error code
     */
    static void error_handler_debug_server(error_level type, std::string msg, std::string where, int error_code) {
#ifndef R_PACKAGE  // avoid std::cerr and std::cout for R package, this is only to reduce R CMD check warnings
        std::string code = (error_code != 0) ? " (" + std::to_string(error_code) + ")" : "";
        std::string where_str = (where.empty()) ? "" : " [in " + where + "]";
        std::string now = "[" + utils::get_curdatetime() + "]";
        if (type == error_level::ERRLVL_ERROR || type == error_level::ERRLVL_FATAL) {
            std::cerr << "ERROR" << code << ": " << msg << where_str << std::endl;
        } else if (type == error_level::ERRLVL_WARNING) {
            std::cout << "WARNING" << code << ": " << msg << where_str << std::endl;
        } else if (type == error_level::ERRLVL_INFO) {
            std::cout << "INFO " << now << ": " << msg << where_str << std::endl;
        } else if (type == error_level::ERRLVL_DEBUG) {
            std::cout << "DEBUG " << now << ": " << msg << where_str << std::endl;
        }
#endif
    }
};

}  // namespace gdalcubes

#endif  // ERROR_H
