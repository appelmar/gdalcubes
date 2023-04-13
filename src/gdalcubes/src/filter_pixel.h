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

#ifndef FILTER_PIXEL_H
#define FILTER_PIXEL_H

#include <algorithm>
#include <string>

#include "cube.h"

struct te_variable;  // forward declaration for add_default_functions

namespace gdalcubes {

/**
     * @brief A data cube that applies one or more arithmetic expressions on band values per pixel
     *
     * @note This class either works either with exprtk or with tinyexpr, depending on whether USE_EXPRTK is defined or not.
     * Please notice that the functionality of these libraries (i.e. the amount of functions they support) may vary. tinyexpr
     * seems to work only with lower case symbols, expressions and band names are automatically converted to lower case then.
     */
class filter_pixel_cube : public cube {
   public:
    /**
         * @brief Create a data cube that applies arithmetic expressions on pixels of an input data cube
         * @note This static creation method should preferably be used instead of the constructors as
         * the constructors will not set connections between cubes properly.
         * @param in input data cube
         * @param expr vector of string expressions, each expression will result in a new band in the resulting cube where values are derived from the input cube according to the specific expression
         * @param band_names specify names for the bands of the resulting cube, if empty, "band1", "band2", "band3", etc. will be used as names
         * @return a shared pointer to the created data cube instance
         */
    static std::shared_ptr<filter_pixel_cube> create(std::shared_ptr<cube> in, std::string predicate) {
        std::shared_ptr<filter_pixel_cube> out = std::make_shared<filter_pixel_cube>(in, predicate);
        in->add_child_cube(out);
        out->add_parent_cube(in);
        return out;
    }

   public:
    /**
         * @brief Create a data cube that applies arithmetic expressions on pixels of an input data cube
         * @param expr vector of string expressions, each expression will result in a new band in the resulting cube where values are derived from the input cube according to the specific expression
         * @param band_names specify names for the bands of the resulting cube, if empty, "band1", "band2", "band3", etc. will be used as names
         */
    filter_pixel_cube(std::shared_ptr<cube> in, std::string predicate) : cube(in->st_reference()->copy()), _in_cube(in), _pred(predicate) {  // it is important to duplicate st reference here, otherwise changes will affect input cube as well
        _chunk_size[0] = _in_cube->chunk_size()[0];
        _chunk_size[1] = _in_cube->chunk_size()[1];
        _chunk_size[2] = _in_cube->chunk_size()[2];

        for (uint16_t i = 0; i < _in_cube->bands().count(); ++i) {
            _bands.add(_in_cube->bands().get(i));
        }

        // tinyexpr works with lower case symbols only

        std::transform(_pred.begin(), _pred.end(), _pred.begin(), ::tolower);

        // parse expressions, currently this is only for validation,
        // expressions will be parsed again in read_chunk(), costs should
        // be negligible compared to the actual evaluation
        if (!parse_predicate()) {
            GCBS_ERROR("Invalid predicate");
            throw std::string("ERROR in filter_pixel_cube::filter_pixel_cube(): Invalid predicate");
        }
    }

   public:
    ~filter_pixel_cube() {}

    std::shared_ptr<chunk_data> read_chunk(chunkid_t id) override;

    json11::Json make_constructible_json() override {
        json11::Json::object out;
        out["cube_type"] = "filter_pixel";
        out["predicate"] = _pred;
        out["in_cube"] = _in_cube->make_constructible_json();
        return out;
    }

   private:
    std::shared_ptr<cube> _in_cube;
    std::string _pred;

    bool parse_predicate();
};

}  // namespace gdalcubes

#endif  //FILTER_PIXEL_H
