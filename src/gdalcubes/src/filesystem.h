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

#ifndef FILESYSTEM_H
#define FILESYSTEM_H

#include <functional>
#include <string>
#include <cstdint> // 2023-01-12: GCC 13 compatibility

#ifdef _WIN32
#define DIR_SEPARATOR "\\"
#else
#define DIR_SEPARATOR "/"
#endif

namespace gdalcubes {

/**
 * @brief A simple wrapper class around GDAL's CPL and VSI interfaces for simple filesystem operations
 */
class filesystem {
   public:
    static bool exists(std::string p);
    static bool is_directory(std::string p);
    static bool is_regular_file(std::string p);
    static std::string stem(std::string p);
    static std::string filename(std::string p);
    static std::string extension(std::string p);
    static std::string directory(std::string p);
    static std::string get_working_dir();
    static std::string make_absolute(std::string p);
    static std::string parent(std::string p);
    static std::string join(std::string p1, std::string p2);
    static void iterate_directory(std::string p, std::function<void(const std::string&)> f);
    static void iterate_directory_recursive(std::string p, std::function<void(const std::string&)> f);
    static void remove(std::string p);
    static void mkdir(std::string p);
    static void mkdir_recursive(std::string p);
    static bool is_relative(std::string p);
    static bool is_absolute(std::string p);
    static std::string get_tempdir();
    static uint32_t file_size(std::string p);
    static void move(std::string src, std::string dest);
    static void copy(std::string src, std::string dest);
};

}  // namespace gdalcubes

#endif  //FILESYSTEM_H
