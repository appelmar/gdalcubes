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

#ifndef REDUCE_TIME_H
#define REDUCE_TIME_H

#include "cube.h"

namespace gdalcubes {

/**
 * @brief A data cube that applies reducer functions over selected bands of a data cube over time
 * @note This is a reimplementation of reduce_cube. The new implementation allows to apply different reducers to different bands instead of just one reducer to all bands of the input data cube
 */
class reduce_time_cube : public cube {
   public:
    /**
     * @brief Create a data cube that applies a reducer function on a given input data cube over time
     * @note This static creation method should preferably be used instead of the constructors as
     * the constructors will not set connections between cubes properly.
     * @param in input data cube
     * @param reducer reducer function
     * @return a shared pointer to the created data cube instance
     */
    static std::shared_ptr<reduce_time_cube> create(std::shared_ptr<cube> in, std::vector<std::pair<std::string, std::string>> reducer_bands) {
        std::shared_ptr<reduce_time_cube> out = std::make_shared<reduce_time_cube>(in, reducer_bands);
        in->add_child_cube(out);
        out->add_parent_cube(in);
        return out;
    }

   public:
    reduce_time_cube(std::shared_ptr<cube> in, std::vector<std::pair<std::string, std::string>> reducer_bands) : cube(in->st_reference()->copy()), _in_cube(in), _reducer_bands(reducer_bands) {  // it is important to duplicate st reference here, otherwise changes will affect input cube as well

        if (cube_stref::type_string(in->st_reference()) == "cube_stref_regular") {
            std::shared_ptr<cube_stref_regular> stref = std::dynamic_pointer_cast<cube_stref_regular>(_st_ref);
            duration dt = (stref->t1() - stref->t0() + 1);
            stref->set_t_axis(stref->t0(), stref->t1(), dt);
        } else if (cube_stref::type_string(_st_ref) == "cube_stref_labeled_time") {
            std::shared_ptr<cube_stref_labeled_time> stref = std::dynamic_pointer_cast<cube_stref_labeled_time>(_st_ref);
            //stref->dt((stref->t1() - stref->t0()) + 1);
            //stref->t1(stref->t0()) ;  // set nt=1
            stref->set_time_labels({stref->t0()});
        }

        assert(_st_ref->nt() == 1);
        _chunk_size[0] = 1;
        _chunk_size[1] = _in_cube->chunk_size()[1];
        _chunk_size[2] = _in_cube->chunk_size()[2];

        // TODO: check for duplicate band, reducer pairs?

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
                  reducerstr == "sum" ||
                  reducerstr == "which_min" ||
                  reducerstr == "which_max" ||
                  reducerstr == "Q1" ||
                  reducerstr == "Q3"))
                throw std::string("ERROR in reduce_time_cube::reduce_time_cube(): Unknown reducer '" + reducerstr + "'");

            if (!(in->bands().has(bandstr))) {
                throw std::string("ERROR in reduce_time_cube::reduce_time_cube(): Input data cube has no band '" + bandstr + "'");
            }

            band b = in->bands().get(bandstr);
            if (in->size_t() > 1) {
                b.name = b.name + "_" + reducerstr;  // Change name only if input is not yet reduced
            }
            _bands.add(b);
        }
    }

   public:
    ~reduce_time_cube() {}

    std::shared_ptr<chunk_data> read_chunk(chunkid_t id) override;

    /**
 * Combines all chunks and produces a single GDAL image
 * @param path path to output image file
 * @param format GDAL format (see https://www.gdal.org/formats_list.html)
 * @param co GDAL create options
     * @param p chunk processor instance, defaults to the current global configuration in config::instance()->get_default_chunk_processor()
 */
    //void write_gdal_image(std::string path, std::string format = "GTiff", std::vector<std::string> co = std::vector<std::string>(), std::shared_ptr<chunk_processor> p = config::instance()->get_default_chunk_processor());

    json11::Json make_constructible_json() override {
        json11::Json::object out;
        out["cube_type"] = "reduce_time";
        json11::Json::array rb;
        for (uint16_t i = 0; i < _reducer_bands.size(); ++i) {
            rb.push_back(json11::Json::array({_reducer_bands[i].first, _reducer_bands[i].second}));
        }
        out["reducer_bands"] = rb;
        out["in_cube"] = _in_cube->make_constructible_json();
        return out;
    }

   private:
    std::shared_ptr<cube> _in_cube;
    std::vector<std::pair<std::string, std::string>> _reducer_bands;
};
}  // namespace gdalcubes

#endif  // REDUCE_TIME_H
