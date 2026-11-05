/*
    MIT License

    Copyright (c) 2026 Marius Appel <marius.appel@hs-bochum.de>

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
#ifndef ZARR_CUBE_H
#define ZARR_CUBE_H

#include "cube.h"

namespace gdalcubes {

/**
 * @brief A data cube that reads data from a Zarr dataset
 *
 * This cube reads data cubes from a single Zarr dataset that has been
 * created with gdalcubes (using cube::write_zarr()).
 * At the moment, it is **not** a general purpose reader for Zarr files from various sources.
 */
class zarr_cube : public cube {
   public:
    /**
     * @brief Create a data cube from a single Zarr dataset
     * @note This static creation method should preferably be used instead of the constructors as
     * the constructors will not set connections between cubes properly.
     * @param path path to a Zarr dataset
     * @param auto_unpack if data values have been packed, offset and scale are applied automatically if true
     * @return a shared pointer to the created data cube instance
     */
    static std::shared_ptr<zarr_cube> create(std::string path, bool auto_unpack = true) {
        return std::make_shared<zarr_cube>(path, auto_unpack);
    }

   public:
    zarr_cube(std::string path, bool auto_unpack = true);

   public:
    ~zarr_cube() {}

    // std::string to_string() override;

    /**
     * @brief Select bands by names
     * @param bands vector of bands to be considered in the cube, if empty, all bands will be selected
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

    json11::Json make_constructible_json() override {
        json11::Json::object out;
        out["cube_type"] = "zarr_cube";
        out["chunk_size"] = json11::Json::array({(int)_chunk_size[0], (int)_chunk_size[1], (int)_chunk_size[2]});
        out["path"] = _path;
        json11::Json::array b;
        for (uint16_t i = 0; i < _band_selection.size(); ++i) {
            b.push_back(_band_selection[i]);
        }
        if (!b.empty()) out["band_selection"] = b;
        out["auto_unpack"] = _auto_unpack;

        return out;
    }

    // zarr_cube allows changing chunk sizes from outside
    void set_chunk_size(uint32_t t, uint32_t y, uint32_t x) {
        _chunk_size = {t, y, x};
    }

   private:
    bool _auto_unpack;
    std::string _path;
    band_collection _orig_bands;
    std::vector<std::string> _band_selection;
    std::mutex _mutex;
};

}  // namespace gdalcubes

#endif  // ZARR_CUBE_H
