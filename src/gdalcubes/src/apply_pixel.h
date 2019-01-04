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

#ifndef APPLY_PIXEL_H
#define APPLY_PIXEL_H

#include "cube.h"
#include <algorithm>
#include <string>

// If possible, take care to NOT include exprtk.hpp in header

/**
 * @brief A data cube that applies one or more arithmetic expressions on band values per pixel
 *
 * @note This class either works either with exprtk or with tinyexpr, depending on whether USE_EXPRTK is defined or not.
 * Please notice that the functionality of these libraries (i.e. the amount of functions they support) may vary. tinyexpr
 * seems to work only with lower case symbols, expressions and band names are automatically converted to lower case then.
 */
class apply_pixel_cube : public cube {
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
    static std::shared_ptr<apply_pixel_cube> create(std::shared_ptr<cube> in, std::vector<std::string> expr, std::vector<std::string> band_names = {}) {
        std::shared_ptr<apply_pixel_cube> out = std::make_shared<apply_pixel_cube>(in, expr, band_names);
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
    apply_pixel_cube(std::shared_ptr<cube> in, std::vector<std::string> expr, std::vector<std::string> band_names = {}) : cube(std::make_shared<cube_st_reference>(*(in->st_reference()))), _in_cube(in), _expr(expr), _band_names(band_names), _band_usage(), _band_usage_all() {  // it is important to duplicate st reference here, otherwise changes will affect input cube as well
        _chunk_size[0] = _in_cube->chunk_size()[0];
        _chunk_size[1] = _in_cube->chunk_size()[1];
        _chunk_size[2] = _in_cube->chunk_size()[2];

        if (_expr.empty()) {
            GCBS_ERROR("No expression given");
            throw std::string("ERROR in apply_pixel_cube::apply_pixel_cube(): No expression given");
        }
        if (!_band_names.empty()) {
            if (_band_names.size() != _expr.size()) {
                GCBS_ERROR("Numbers of provided expressions and band names differ");
                throw std::string("ERROR in apply_pixel_cube::apply_pixel_cube(): Numbers of provided expressions and band names differ");
            }
        }

        for (uint16_t i = 0; i < _expr.size(); ++i) {
            band b(_band_names.empty() ? ("band" + std::to_string(i + 1)) : _band_names[i]);
            b.unit = "";
            b.no_data_value = "nan";
            b.type = "float64";
            b.offset = 0;
            b.scale = 1;
            _bands.add(b);
        }

#ifndef USE_EXPRTK // tinyexpr works with lower case symbols only
        for (uint16_t i = 0; i < _expr.size(); ++i) {
            std::transform(_expr[i].begin(), _expr[i].end(), _expr[i].begin(), ::tolower);
        }
#endif

        // parse expressions, currently this is only for validation,
        // expressions will be parsed again in read_chunk(), costs should
        // be negligible compared to the actual evaluation
        if (!parse_expressions()) {
            GCBS_ERROR("Invalid expression(s)");
            throw std::string("ERROR in apply_pixel_cube::apply_pixel_cube(): Invalid expression(s)");
        }

        // Find out, which expressions are actually used per expression
        for (uint16_t i = 0; i < _expr.size(); ++i) {
            _band_usage.push_back(std::set<std::string>());
            for (uint16_t ib = 0; ib < _in_cube->bands().count(); ++ib) {
                std::string name = _in_cube->bands().get(ib).name;
                std::string temp_name = name;
#ifndef USE_EXPRTK // tinyexpr works with lower case symbols only
                std::transform(temp_name.begin(), temp_name.end(), temp_name.begin(), ::tolower);
#endif
                if (_expr[i].find(temp_name) != std::string::npos) {
                    _band_usage[i].insert(name);
                    _band_usage_all.insert(name);
                }
            }
        }
    }

   public:
    ~apply_pixel_cube() {}

    std::shared_ptr<chunk_data> read_chunk(chunkid_t id) override;

    nlohmann::json make_constructible_json() override {
        nlohmann::json out;
        out["cube_type"] = "apply_pixel";
        out["expr"] = _expr;
        if (!_band_names.empty())
            out["band_names"] = _band_names;
        out["in_cube"] = _in_cube->make_constructible_json();
        return out;
    }

   private:
    std::shared_ptr<cube> _in_cube;
    std::vector<std::string> _expr;
    std::vector<std::string> _band_names;
    std::vector<std::set<std::string>> _band_usage;  // store which bands are really used per expression
    std::set<std::string> _band_usage_all;

    virtual void set_st_reference(std::shared_ptr<cube_st_reference> stref) override {
        _st_ref->win() = stref->win();
        _st_ref->proj() = stref->proj();
        _st_ref->ny() = stref->ny();
        _st_ref->nx() = stref->nx();
        _st_ref->t0() = stref->t0();
        _st_ref->t1() = stref->t1();
        _st_ref->dt() = stref->dt();
    }

    bool parse_expressions();
};

#endif  //APPLY_PIXEL_H
