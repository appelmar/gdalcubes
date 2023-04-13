/*
    MIT License

    Copyright (c) 2022 Marius Appel <marius.appel@uni-muenster.de>

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

#ifndef AGGREGATE_SPACE_H
#define AGGREGATE_SPACE_H

#include "cube.h"

namespace gdalcubes {



/**
 * @brief A data cube that applies reducer functions over selected bands of a data cube over time
 * @note This is a reimplementation of reduce_cube. The new implementation allows to apply different reducers to different bands instead of just one reducer to all bands of the input data cube
 */
class aggregate_space_cube : public cube {
   public:

    /**
     * @brief Create a data cube that aggregates pixel time series to lower temporal resolution
     * @note This static creation method should preferably be used instead of the constructors as
     * the constructors will not set connections between cubes properly.
     * @param in input data cube
     * @param func aggregation function
     * @param fact number of cells that become aggregated to a new cell
     * @return a shared pointer to the created data cube instance
     */
    static std::shared_ptr<aggregate_space_cube> create(std::shared_ptr<cube> in, double dx, double dy, std::string func = "mean") {
        std::shared_ptr<aggregate_space_cube> out = std::make_shared<aggregate_space_cube>(in, dx, dy, func);

        // TODO: labeled time axis?
        in->add_child_cube(out);
        out->add_parent_cube(in);
        return out;
    }

    /**
     * @brief Create a data cube that aggregates pixel time series to lower temporal resolution
     * @note This static creation method should preferably be used instead of the constructors as
     * the constructors will not set connections between cubes properly.
     * @param in input data cube
     * @param func aggregation function
     * @param dt new temporal duration of a data cube pixel
     * @return a shared pointer to the created data cube instance
     */
    static std::shared_ptr<aggregate_space_cube> create(std::shared_ptr<cube> in, uint32_t fact, std::string func = "mean") {
        std::shared_ptr<aggregate_space_cube> out = std::make_shared<aggregate_space_cube>(in, double(fact) * in->st_reference()->dx(),
                                                                                           double(fact) * in->st_reference()->dy(), func);
        in->add_child_cube(out);
        out->add_parent_cube(in);
        return out;
    }

   public:
    aggregate_space_cube(std::shared_ptr<cube> in, double dx, double dy, std::string func = "mean") : cube(),_in_cube(in),
                                                                                               _in_func(func), _in_dx(dx),
                                                                                               _in_dy(dy) {


        if (!(func == "min" ||
            func == "max" ||
            func == "mean" ||
            func == "median" ||
            func == "count" ||
            func == "var" ||
            func == "sd" ||
            func == "prod" ||
            func == "sum")) {
                throw std::string("ERROR in aggregate_space_cube::aggregate_space_cube(): unknown aggregation function '" + func + "'");
            }



        std::shared_ptr<cube_stref> stref = in->st_reference()->copy();

        // TODO: check that _in_dx > old dx (and y accordingly)
        std::dynamic_pointer_cast<cube_stref_regular>(stref)->set_x_axis(in->st_reference()->left(), in->st_reference()->right(), _in_dx);
        std::dynamic_pointer_cast<cube_stref_regular>(stref)->set_y_axis(in->st_reference()->bottom(), in->st_reference()->top(), _in_dy);

        _st_ref = stref;
        _chunk_size[0] = _in_cube->chunk_size()[0];
        _chunk_size[1] =  std::min(_st_ref->ny(), _in_cube->chunk_size()[1]);
        _chunk_size[2] =  std::min(_st_ref->nx(), _in_cube->chunk_size()[2]);

        for (uint16_t i = 0; i < in->bands().count(); ++i) {
            _bands.add(in->bands().get(i));
        }
    }

   public:
    ~aggregate_space_cube() {}

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
        out["cube_type"] = "aggregate_space";
        out["func"] = _in_func;
        out["dx"] = _in_dx;
        out["dy"] = _in_dy;
        out["in_cube"] = _in_cube->make_constructible_json();
        return out;
    }

   private:
    std::shared_ptr<cube> _in_cube;
    std::string _in_func;
    double _in_dx;
    double _in_dy;
};
}  // namespace gdalcubes

#endif  // AGGREGATE_SPACE_H
