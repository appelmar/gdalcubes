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

#ifndef AGGREGATE_TIME_H
#define AGGREGATE_TIME_H

#include "cube.h"

namespace gdalcubes {



/**
 * @brief A data cube that applies reducer functions over selected bands of a data cube over time
 * @note This is a reimplementation of reduce_cube. The new implementation allows to apply different reducers to different bands instead of just one reducer to all bands of the input data cube
 */
class aggregate_time_cube : public cube {
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
    static std::shared_ptr<aggregate_time_cube> create(std::shared_ptr<cube> in, std::string dt, std::string func = "mean") {
        std::shared_ptr<aggregate_time_cube> out = std::make_shared<aggregate_time_cube>(in, dt, func);

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
    static std::shared_ptr<aggregate_time_cube> create(std::shared_ptr<cube> in, uint32_t fact, std::string func = "mean") {
        if (!in->st_reference()->has_regular_time()) {
            GCBS_ERROR("Aggregation of data cubes works only by providing a new datetime duration instead of fact");
            throw std::string("Aggregation of data cubes works only by providing a new datetime duration instead of fact");
        }
        duration dt = in->st_reference()->dt();
        dt.dt_interval = (int32_t)fact * dt.dt_interval;
        std::shared_ptr<aggregate_time_cube> out = std::make_shared<aggregate_time_cube>(in, dt.to_string(), func);
        in->add_child_cube(out);
        out->add_parent_cube(in);
        return out;
    }

   public:
    aggregate_time_cube(std::shared_ptr<cube> in, std::string dt, std::string func = "mean") : cube(),_in_cube(in),
                                                                                               _in_func(func), _in_dt(dt),
                                                                                               _dt() {


        if (!(func == "min" ||
            func == "max" ||
            func == "mean" ||
            func == "median" ||
            func == "count" ||
            func == "var" ||
            func == "sd" ||
            func == "prod" ||
            func == "sum")) {
                throw std::string("ERROR in aggregate_time_cube::aggregate_time_cube(): unknown aggregation function '" + func + "'");
            }

        _dt = duration::from_string(dt);


        std::shared_ptr<cube_stref_regular> stref = std::make_shared<cube_stref_regular>();

        stref->srs(in->st_reference()->srs());
        stref->set_x_axis(in->st_reference()->left(), in->st_reference()->right(), in->st_reference()->dx());
        stref->set_y_axis(in->st_reference()->bottom(), in->st_reference()->top(), in->st_reference()->dy());
        if (cube_stref::type_string(in->st_reference()) == "cube_stref_regular") {
            stref->set_t_axis(std::dynamic_pointer_cast<cube_stref_regular>(in->st_reference())->t0(),
                              std::dynamic_pointer_cast<cube_stref_regular>(in->st_reference())->t1(),
                              _dt);
        }
        else if (cube_stref::type_string(in->st_reference()) == "cube_stref_labeled_time") {
            stref->set_t_axis(std::dynamic_pointer_cast<cube_stref_labeled_time>(in->st_reference())->t0(),
                              std::dynamic_pointer_cast<cube_stref_labeled_time>(in->st_reference())->t1(),
                              _dt);
        }


//        if (stref->nt() > _in_cube->st_reference()->nt()) {
//            GCBS_ERROR("ERROR in aggregate_time_cube::aggregate_time_cube(): resulting cube would have more pixels than the input cube, please change target resolution");
//            throw std::string("ERROR in aggregate_time_cube::aggregate_time_cube(): resulting cube would have more pixels than the input cube, please change target resolution");
//        }

        _st_ref = stref;
        _chunk_size[0] = std::min(_st_ref->nt(), _in_cube->chunk_size()[0]);
        _chunk_size[1] = _in_cube->chunk_size()[1];
        _chunk_size[2] = _in_cube->chunk_size()[2];

        for (uint16_t i = 0; i < in->bands().count(); ++i) {
            _bands.add(in->bands().get(i));
        }

        // check whether cell boundaries of input cube are aligned with
        // cell boundaries of result cube
        if (cube_stref::type_string(_in_cube->st_reference()) == "cube_stref_regular") {
            for (uint32_t it=0; it<in->size_t(); ++it) {
                datetime t_cur = _st_ref->datetime_at_index(it);
                datetime t_next =  _st_ref->datetime_at_index( it + 1);

                t_cur.unit(_in_cube->st_reference()->dt_unit());
                t_next.unit( _in_cube->st_reference()->dt_unit());

                uint32_t first = _in_cube->st_reference()->index_at_datetime(t_cur);
                if (_in_cube->st_reference()->datetime_at_index(first) != t_cur) {
                    GCBS_WARN("Some cells of the aggregated cube temporally overlap with two cells of the input cube "
                        "(their boundaries do not align). Aggregation function will select cells of the input cube based on "
                        "their starting datetime. If this is not desired, please change the temporal resolution and/or "
                        "aggregation size to yield aligned cell boundaries.");
                    break; // only once
                }
            }
        }
    }

   public:
    ~aggregate_time_cube() {}

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
        out["cube_type"] = "aggregate_time";
        out["func"] = _in_func;
        out["dt"] = _in_dt;
        out["in_cube"] = _in_cube->make_constructible_json();
        return out;
    }

   private:
    std::shared_ptr<cube> _in_cube;
    std::string _in_func;
    std::string _in_dt;

    duration _dt;
};
}  // namespace gdalcubes

#endif  // AGGREGATE_TIME_H
