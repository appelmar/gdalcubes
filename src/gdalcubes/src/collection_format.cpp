/*
    MIT License

    Copyright (c) 2019 Marius Appel <marius.appel@hs-bochum.de>

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

#include "collection_format.h"

#include <fstream>
#include <sstream>

#include "filesystem.h"

namespace gdalcubes {

bool collection_format::is_null() {
    return _j.is_null();
}

std::map<std::string, std::string> collection_format::list_presets() {
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

void collection_format::load_file(std::string filename) {
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
    std::stringstream str;
    str << i.rdbuf();

    std::string err;  // TODO: do something with err
    _j = json11::Json::parse(str.str(), err);
}

void collection_format::load_string(std::string jsonstr) {
    std::string err;  // TODO: do something with err
    _j = json11::Json::parse(jsonstr, err);
}

}  // namespace gdalcubes
