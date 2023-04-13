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

#ifndef COLLECTION_FORMAT_H
#define COLLECTION_FORMAT_H

#include <string>

#include "config.h"
#include "external/json11/json11.hpp"

namespace gdalcubes {

/**
 * The image collection format describes rules how image collections can be build from a simple list of files / URLs.
 * It defines
 * - how to identify which files correspond to the same image,
 * - the bands of the image collection and in which files they are stored,
 * - and how to derive the date / time from images.
 * The format is specified in a simple JSON document with regular expressions on filenames / URLs.
 * gdalcubes comes with a set of predefined collection formats for typical Earth observation datasets including
 * Sentinel 2, Landsat 8, and  SRTM.
 */
class collection_format {
   public:
    collection_format() {}
    collection_format(std::string filename) { load_file(filename); }

    bool is_null();

    static std::map<std::string, std::string> list_presets();

    /**
     * Construct a collection format from a JSON file
     * @param filename
     */
    void load_file(std::string filename);

    /**
    * @brief Construct a collection format from a JSON string
    * @param jsonstr JSON string
    */
    void load_string(std::string jsonstr);

    /**
     * Returns the raw json document.
     * @return JSON object from json library (see https://github.com/nlohmann/json)
     */
    inline json11::Json& json() { return _j; }

   private:
    json11::Json _j;
};

}  // namespace gdalcubes

#endif  //COLLECTION_FORMAT_H
