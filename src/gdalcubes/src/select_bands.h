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

#ifndef SELECT_BANDS_H
#define SELECT_BANDS_H

#include "cube.h"
#include "image_collection_cube.h"

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
    select_bands_cube(std::shared_ptr<cube> in, std::vector<std::string> bands) : cube(std::make_shared<cube_st_reference>(*(in->st_reference()))), _in_cube(in), _band_sel(bands), _input_is_image_collection_cube(false) {  // it is important to duplicate st reference here, otherwise changes will affect input cube as well
        _chunk_size[0] = _in_cube->chunk_size()[0];
        _chunk_size[1] = _in_cube->chunk_size()[1];
        _chunk_size[2] = _in_cube->chunk_size()[2];

        if (std::dynamic_pointer_cast<image_collection_cube>(in)) {
            _input_is_image_collection_cube = true;
            std::dynamic_pointer_cast<image_collection_cube>(in)->select_bands(bands);
        }

        for (uint16_t ib = 0; ib < _band_sel.size(); ++ib) {
            if (!in->bands().has(_band_sel[ib])) {
                GCBS_ERROR("Input cube has no band '" + _band_sel[ib] + "'");
                throw std::string("ERROR in select_bands_cube::select_bands_cube(): Input cube has no band '" + _band_sel[ib] + "'");
            }
            _bands.add(in->bands().get(_band_sel[ib]));
        }
    }

    select_bands_cube(std::shared_ptr<cube> in, std::vector<uint16_t> bands) : cube(std::make_shared<cube_st_reference>(*(in->st_reference()))), _in_cube(in), _band_sel(), _input_is_image_collection_cube(false) {  // it is important to duplicate st reference here, otherwise changes will affect input cube as well
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
            _input_is_image_collection_cube = true;
            std::dynamic_pointer_cast<image_collection_cube>(in)->select_bands(bands);
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

    nlohmann::json make_constructible_json() override {
        nlohmann::json out;
        out["cube_type"] = "select_bands";
        out["bands"] = _band_sel;
        out["in_cube"] = _in_cube->make_constructible_json();
        return out;
    }

   private:
    std::shared_ptr<cube> _in_cube;
    std::vector<std::string> _band_sel;
    bool _input_is_image_collection_cube;

    virtual void set_st_reference(std::shared_ptr<cube_st_reference> stref) override {
        // copy fields from st_reference type
        _st_ref->win() = stref->win();
        _st_ref->proj() = stref->proj();
        _st_ref->ny() = stref->ny();
        _st_ref->nx() = stref->nx();
        _st_ref->t0() = stref->t0();
        _st_ref->t1() = stref->t1();
        _st_ref->dt() = stref->dt();
    }
};

#endif  //SELECT_BANDS_H
