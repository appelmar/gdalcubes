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

#ifndef REDUCE_SPACE_H
#define REDUCE_SPACE_H

#include "cube.h"

namespace gdalcubes {

/**
 * @brief A data cube that applies reducer functions over selected bands of a data cube over space
 */
class reduce_space_cube : public cube {
   public:
    /**
     * @brief Create a data cube that applies a reducer function on a given input data cube over space
     * @note This static creation method should preferably be used instead of the constructors as
     * the constructors will not set connections between cubes properly.
     * @param in input data cube
     * @param reducer reducer function
     * @param names names for new bands
     * @return a shared pointer to the created data cube instance
     */
    static std::shared_ptr<reduce_space_cube> create(std::shared_ptr<cube> in, std::vector<std::pair<std::string, std::string>> reducer_bands, std::vector<std::string> names) {
        std::shared_ptr<reduce_space_cube> out = std::make_shared<reduce_space_cube>(in, reducer_bands, names);
        in->add_child_cube(out);
        out->add_parent_cube(in);
        return out;
    }
     
    static std::shared_ptr<reduce_space_cube> create(std::shared_ptr<cube> in, std::vector<std::pair<std::string, std::string>> reducer_bands) {
       std::shared_ptr<reduce_space_cube> out = std::make_shared<reduce_space_cube>(in, reducer_bands, std::vector<std::string>());
       in->add_child_cube(out);
       out->add_parent_cube(in);
       return out;
     }

   public:
    reduce_space_cube(std::shared_ptr<cube> in, std::vector<std::pair<std::string, std::string>> reducer_bands, std::vector<std::string> names) : cube(in->st_reference()->copy()), _in_cube(in), _reducer_bands(reducer_bands), _names(names) {  // it is important to duplicate st reference here, otherwise changes will affect input cube as well
        if (cube_stref::type_string(_st_ref) == "cube_stref_regular") {
            std::shared_ptr<cube_stref_regular> stref = std::dynamic_pointer_cast<cube_stref_regular>(_st_ref);
            stref->set_x_axis(_st_ref->left(), _st_ref->right(), (uint32_t)1);
            stref->set_y_axis(_st_ref->bottom(), _st_ref->top(), (uint32_t)1);
        } else if (cube_stref::type_string(_st_ref) == "cube_stref_labeled_time") {
            std::shared_ptr<cube_stref_labeled_time> stref = std::dynamic_pointer_cast<cube_stref_labeled_time>(_st_ref);
            stref->set_x_axis(_st_ref->left(), _st_ref->right(), (uint32_t)1);
            stref->set_y_axis(_st_ref->bottom(), _st_ref->top(), (uint32_t)1);
        }
        assert(_st_ref->nx() == 1 && _st_ref->ny() == 1);

        _chunk_size[0] = _in_cube->chunk_size()[0];
        _chunk_size[1] = 1;
        _chunk_size[2] = 1;

        // TODO: check for duplicate band, reducer pairs?

        if (!names.empty()) {
          if (names.size() != reducer_bands.size())
            throw std::string("ERROR in reduce_space_cube::reduce_space_cube(): The number of provided names must match the number of expressions");
        }
        for (uint16_t i = 0; i < reducer_bands.size(); ++i) {
            std::string reducerstr = reducer_bands[i].first;
            std::string bandstr = reducer_bands[i].second;
            if (!(reducerstr == "min" ||
                  reducerstr == "max" ||
                  reducerstr == "mean" ||
                  reducerstr == "median" ||
                  reducerstr == "count" ||
                  reducerstr == "var" ||
                  reducerstr == "sd" ||
                  reducerstr == "prod" ||
                  reducerstr == "sum"))
                throw std::string("ERROR in reduce_space_cube::reduce_space_cube(): Unknown reducer '" + reducerstr + "'");

            if (!(in->bands().has(bandstr))) {
                throw std::string("ERROR in reduce_space_cube::reduce_space_cube(): Input data cube has no band '" + bandstr + "'");
            }

            band b = in->bands().get(bandstr);
            if (names.empty()) {
              if (in->size_x() > 1 || in->size_y() > 1) {
                  b.name = b.name + "_" + reducerstr;  // Change name only if input is not yet reduced
              }
            }
            else {
              b.name = names[i];
            }
            _bands.add(b);
        }
    }

   public:
    ~reduce_space_cube() {}

    std::shared_ptr<chunk_data> read_chunk(chunkid_t id) override;

    /**
 * Combines all chunks and produces a single GDAL image
 * @param path path to output image file
 * @param format GDAL format (see https://www.gdal.org/formats_list.html)
 * @param co GDAL create options
     * @param p chunk processor instance, defaults to the current global configuration in config::instance()->get_default_chunk_processor()
 */

    json11::Json make_constructible_json() override {
        json11::Json::object out;
        out["cube_type"] = "reduce_space";
        json11::Json::array rb;
        for (uint16_t i = 0; i < _reducer_bands.size(); ++i) {
            rb.push_back(json11::Json::array({_reducer_bands[i].first, _reducer_bands[i].second}));
        }
        if (!_names.empty()) {
          json11::Json::array rn;
          for (uint16_t i = 0; i < _names.size(); ++i) {
            rn.push_back(_names[i]);
          }
          out["names"] = rn;
        }
        out["reducer_bands"] = rb;
        out["in_cube"] = _in_cube->make_constructible_json();
        return out;
    }

   private:
    std::shared_ptr<cube> _in_cube;
    std::vector<std::pair<std::string, std::string>> _reducer_bands;
    std::vector<std::string> _names;
};

}  // namespace gdalcubes

#endif  // REDUCE_SPACE_H
