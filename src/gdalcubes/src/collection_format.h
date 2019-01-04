/*
   Copyright 2018 Marius Appel <marius.appel@uni-muenster.de>

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.
*/

#ifndef COLLECTION_FORMAT_H
#define COLLECTION_FORMAT_H

#include <sqlite3.h>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include "config.h"
#include "external/json.hpp"
#include "filesystem.h"

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

    bool is_null() {
        return _j.empty();
    }

    static std::map<std::string, std::string> list_presets() {
        std::map<std::string, std::string> out;

        std::vector<std::string> dirs = config::instance()->get_collection_format_preset_dirs();

        // do not use uint here because of descending iteration
        for (int i = dirs.size() - 1; i >= 0; --i) {
            if (!filesystem::exists(dirs[i])) {
                continue;
            }

            filesystem::iterate_directory(dirs[i], [&out](const std::string& p) {
                if (filesystem::is_regular_file(p) && filesystem::extension(p) == "json") {
                    if (out.find(filesystem::stem(p)) == out.end())
                        out.insert(std::pair<std::string, std::string>(filesystem::stem(p), filesystem::make_absolute(p)));
                }
            });
        }
        return out;
    }

    /**
     * Construct a collection format from a JSON file
     *
     *
     * @param filename
     */
    void load_file(std::string filename) {
        _j.clear();

        if ((!filesystem::exists(filename) && !filesystem::is_absolute(filename) && filesystem::parent(filename).empty()) || filesystem::is_directory(filename)) {
            // simple filename without directories
            GCBS_DEBUG("Couldn't find collection format '" + filename + "', looking for a preset with the same name");
            std::map<std::string, std::string> preset_formats = list_presets();
            if (preset_formats.find(filesystem::stem(filename)) != preset_formats.end()) {
                filename = preset_formats[filesystem::stem(filename)];
                GCBS_DEBUG("Found collection format preset at '" + filename + "'");
            }
        }
        if (!filesystem::exists(filename) || filesystem::is_directory(filename))
            throw std::string("ERROR in collection_format::load_file(): image collection format file does not exist.");

        std::ifstream i(filename);
        i >> _j;
    }

    /**
    * @brief Construct a collection format from a JSON string
    * @param jsonstr JSON string
    */
    void load_string(std::string jsonstr) {
        _j.clear();
        std::istringstream ss(jsonstr);
        ss >> _j;
    }

    /**
     * Validates the formal structure of the JSON format description. Validation is only formal,
     * i.e. formal such that the presence of mandatory properties is checked but not the meaning of their values.
     * @note **NOT YET IMPLEMENTED**
     * @return true, if the JSON format description is valid, false otherwise.
     */
    bool validate();

    /**
     * Apply the collection format to a list of filenames or other GDALataset descriptors and build the basic image collection data structure in
     * an SQLite database. This function will not open any files.
     * @note Reading date / time from GDALInfo metadata is currently not implemented
     * @param descriptors list of GDAL dataset descriptors (typically filenames or URLs but can sometimes be database object descriptors or similar)
     * @param db_filename file name of the created result database, if empty string, keep a temporary database, which might be stored in memory.
     * @param strict if true, the creation will cancel as soon as one descriptor produces an error. If false such filles will be skipped.
     * @return SQLite handle of the resulting database, this must eventually be closed by the calling function using `sqlite3_close(db)`
     */
    //sqlite3* apply(const std::vector<std::string>& descriptors, const std::string& db_filename = "", bool strict=true);

    /**
     * Returns the raw json document.
     * @return JSON object from json library (see https://github.com/nlohmann/json)
     */
    inline nlohmann::json& json() { return _j; }

   private:
    nlohmann::json _j;
};

#endif  //COLLECTION_FORMAT_H
