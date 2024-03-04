/*
    MIT License

    Copyright (c) 2021 Marius Appel <marius.appel@hs-bochum.de>

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
#ifndef SIMPLE_CUBE_H
#define SIMPLE_CUBE_H

#include <vector>

#include "cube.h"

namespace gdalcubes {

class simple_cube : public cube {
   public:
    /**
      * @brief Create a data cube from a list of image files, their datetime, and possibly bands
      * @note This static creation method should preferably be used instead of the constructors as
      * the constructors will not set connections between cubes properly.
      * @param files vector of filenames pointing to images
      * @param datetime string vector of datetime values, must have the same number of elements as files
      * @param bands string vector of band names if different bands from an image come from different files,
      * must have the same number of elements as files; by default (empty vector), all files are assumed to contain the same bands
      * @param band_names vector of names for bands, applicalbe only if all images contain the same bands (bands is empty)
      * @param dx target pixel size in x direction; if <= 0 (default), highest original resolution from files is used automatically
      * @param dy target pixel size in y direction; if <= 0 (default), highest original resolution from files is used automatically
      * @note Notice that all images must have identical spatial reference systems and spatial extents. No
      * automatic check is performed whether this assumption is fulfilled.
      * @return a shared pointer to the created data cube instance
      */
    static std::shared_ptr<simple_cube> create(std::vector<std::string> files, std::vector<std::string> datetime_values,
                                               std::vector<std::string> bands = {}, std::vector<std::string> band_names = {},
                                               double dx = -1, double dy = -1) {
        return std::make_shared<simple_cube>(files, datetime_values, bands, band_names, dx, dy);
    }

   public:
    simple_cube(std::vector<std::string> files, std::vector<std::string> datetime_values,
                std::vector<std::string> bands = {}, std::vector<std::string> band_names = {},
                double dx = -1, double dy = -1);

   public:
    ~simple_cube() {}

    std::string to_string() override;

    /**
     * @brief Select bands by names
     * @param bands vector of bands to be considered in the cube, if empty, all bands will be selected
     * TODO: implement for simple cube type
     */
    void select_bands(std::vector<std::string> bands) {
        _band_selection.clear();
        if (bands.empty()) {
            _bands = _orig_bands;
        } else {
            band_collection bands_new;
            for (uint16_t i = 0; i < bands.size(); ++i) {
                if (_orig_bands.has(bands[i])) {
                    bands_new.add(_orig_bands.get(bands[i]));
                    _band_selection.push_back(bands[i]);
                } else {
                    GCBS_WARN("Data cube has no band with name '" + bands[i] + "'; band will be skipped");
                }
            }
            if (bands_new.count() > 0) {
                _bands = bands_new;
            } else {
                _bands = _orig_bands;
            }
        }
    }

    std::shared_ptr<chunk_data> read_chunk(chunkid_t id) override;

    // simple cube allows changing chunk sizes from outside!
    // This is important for e.g. streaming.
    void set_chunk_size(uint32_t t, uint32_t y, uint32_t x) {
        _chunk_size = {t, y, x};
    }

    void set_strict(bool s) {
        _strict = s;
    }

    json11::Json make_constructible_json() override {
        json11::Json::object out;
        out["cube_type"] = "simple_cube";
        out["chunk_size"] = json11::Json::array({(int)_chunk_size[0], (int)_chunk_size[1], (int)_chunk_size[2]});

        json11::Json::array arr_files;
        for (uint16_t i = 0; i < _in_files.size(); ++i) {
            arr_files.push_back(_in_files[i]);
        }
        out["files"] = arr_files;

        json11::Json::array arr_datetime;
        for (uint16_t i = 0; i < _in_datetime.size(); ++i) {
            arr_datetime.push_back(_in_datetime[i]);
        }
        out["datetime"] = arr_datetime;

        if (!_in_bands.empty()) {
            json11::Json::array arr_bands;
            for (uint16_t i = 0; i < _in_bands.size(); ++i) {
                arr_bands.push_back(_in_bands[i]);
            }
            out["bands"] = arr_bands;
        }

        if (!_in_band_names.empty()) {
            json11::Json::array arr_band_names;
            for (uint16_t i = 0; i < _in_band_names.size(); ++i) {
                arr_band_names.push_back(_in_band_names[i]);
            }
            out["band_names"] = arr_band_names;
        }

        out["dx"] = _in_dx;
        out["dy"] = _in_dy;
        out["strict"] = _strict;

        return out;
    }

   private:
    // Input arguments
    std::vector<std::string> _in_files;
    std::vector<std::string> _in_datetime;
    std::vector<std::string> _in_bands;
    std::vector<std::string> _in_band_names;
    double _in_dx;
    double _in_dy;
    bool _strict;

    std::map<datetime, std::map<std::string, std::pair<std::string, uint16_t>>> _index;
    band_collection _orig_bands;
    std::vector<std::string> _band_selection;
};

}  // namespace gdalcubes

#endif  //SIMPLE_CUBE_H
