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

#ifndef SELECT_BANDS_H
#define SELECT_BANDS_H

#include "cube.h"
#include "image_collection_cube.h"
#include "ncdf_cube.h"
#include "simple_cube.h"

namespace gdalcubes {

/**
 * @brief A data cube that subsets and / or reorders bands of an input data cube
 * @note Where possible select_bands should be called directly on image_collection_cube objects as input
 * because this allows to reduce GDAL RasterIO calls
 */
class select_bands_cube : public cube {
   public:
    /**
     * @brief Create a data cube that subsets bands of an input data cube
     * @note This static creation method should preferably be used instead of the constructors as
     * the constructors will not set connections between cubes properly.
     * @param in input data cube
     * @param bands selected bands given by name
     * @return a shared pointer to the created data cube instance
     */
    static std::shared_ptr<select_bands_cube> create(std::shared_ptr<cube> in, std::vector<std::string> bands) {
        std::shared_ptr<select_bands_cube> out = std::make_shared<select_bands_cube>(in, bands);
        in->add_child_cube(out);
        out->add_parent_cube(in);
        return out;
    }

    /**
     * @brief Create a data cube that subsets bands of an input data cube
     * @note This static creation method should preferably be used instead of the constructors as
     * the constructors will not set connections between cubes properly.
     * @param in input data cube
     * @param bands selected bands given by index
     * @return a shared pointer to the created data cube instance
     */
    static std::shared_ptr<select_bands_cube> create(std::shared_ptr<cube> in, std::vector<uint16_t> bands) {
        std::shared_ptr<select_bands_cube> out = std::make_shared<select_bands_cube>(in, bands);
        in->add_child_cube(out);
        out->add_parent_cube(in);
        return out;
    }

   public:
    select_bands_cube(std::shared_ptr<cube> in, std::vector<std::string> bands) : cube(in->st_reference()->copy()), _in_cube(in), _band_sel(bands), _defer_to_input_cube(false) {  // it is important to duplicate st reference here, otherwise changes will affect input cube as well
        _chunk_size[0] = _in_cube->chunk_size()[0];
        _chunk_size[1] = _in_cube->chunk_size()[1];
        _chunk_size[2] = _in_cube->chunk_size()[2];

        if (std::dynamic_pointer_cast<image_collection_cube>(in)) {
            _defer_to_input_cube = true;
            std::dynamic_pointer_cast<image_collection_cube>(in)->select_bands(bands);
        } else if (std::dynamic_pointer_cast<ncdf_cube>(in)) {
            _defer_to_input_cube = true;
            std::dynamic_pointer_cast<ncdf_cube>(in)->select_bands(bands);
        } else if (std::dynamic_pointer_cast<simple_cube>(in)) {
            _defer_to_input_cube = true;
            std::dynamic_pointer_cast<simple_cube>(in)->select_bands(bands);
        }

        for (uint16_t ib = 0; ib < _band_sel.size(); ++ib) {
            if (!in->bands().has(_band_sel[ib])) {
                GCBS_ERROR("Input cube has no band '" + _band_sel[ib] + "'");
                throw std::string("ERROR in select_bands_cube::select_bands_cube(): Input cube has no band '" + _band_sel[ib] + "'");
            }
            _bands.add(in->bands().get(_band_sel[ib]));
        }
    }

    select_bands_cube(std::shared_ptr<cube> in, std::vector<uint16_t> bands) : cube(in->st_reference()->copy()), _in_cube(in), _band_sel(), _defer_to_input_cube(false) {  // it is important to duplicate st reference here, otherwise changes will affect input cube as well
        _chunk_size[0] = _in_cube->chunk_size()[0];
        _chunk_size[1] = _in_cube->chunk_size()[1];
        _chunk_size[2] = _in_cube->chunk_size()[2];

        std::vector<std::string> bands_str;
        for (uint16_t ib = 0; ib < bands.size(); ++ib) {
            if (bands[ib] < 0 || bands[ib] >= in->bands().count()) {
                GCBS_ERROR("Input cube has no band '" + std::to_string(bands[ib]) + "'");
                throw std::string("ERROR in select_bands_cube::select_bands_cube(): Input cube has no band '" + std::to_string(bands[ib]) + "'");
            }
            bands_str.push_back(in->bands().get(bands[ib]).name);
        }
        _band_sel = bands_str;

        if (std::dynamic_pointer_cast<image_collection_cube>(in)) {
            _defer_to_input_cube = true;
            std::dynamic_pointer_cast<image_collection_cube>(in)->select_bands(bands);
        } else if (std::dynamic_pointer_cast<ncdf_cube>(in)) {
            _defer_to_input_cube = true;
            std::dynamic_pointer_cast<ncdf_cube>(in)->select_bands(bands_str);
        } else if (std::dynamic_pointer_cast<simple_cube>(in)) {
            _defer_to_input_cube = true;
            std::dynamic_pointer_cast<simple_cube>(in)->select_bands(bands_str);
        }

        for (uint16_t ib = 0; ib < _band_sel.size(); ++ib) {
            if (!in->bands().has(_band_sel[ib])) {
                GCBS_ERROR("Input cube has no band '" + _band_sel[ib] + "'");
                throw std::string("ERROR in select_bands_cube::select_bands_cube(): Input cube has no band '" + _band_sel[ib] + "'");
            }
            _bands.add(in->bands().get(_band_sel[ib]));
        }
    }

   public:
    ~select_bands_cube() {}

    std::shared_ptr<chunk_data> read_chunk(chunkid_t id) override;

    json11::Json make_constructible_json() override {
        json11::Json::object out;
        out["cube_type"] = "select_bands";
        out["bands"] = _band_sel;
        out["in_cube"] = _in_cube->make_constructible_json();
        return out;
    }

   private:
    std::shared_ptr<cube> _in_cube;
    std::vector<std::string> _band_sel;
    bool _defer_to_input_cube;
};

}  // namespace gdalcubes

#endif  //SELECT_BANDS_H
