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

#ifndef CROP_H
#define CROP_H

#include "cube.h"

namespace gdalcubes {

/**
 * @brief A data cube that crops space and/or time dimensions of a source data cube
 **/
class crop_cube : public cube {
   public:
    /**
     * @brief Create a data cube crops a source data cube by spatial and/or temporal ranges
     * @note This static creation method should preferably be used instead of the constructors as
     * the constructors will not set connections between cubes properly.
     * @param in input data cube
     * @param ix_min integer index of the first cell in the x dimension (left)
     * @param ix_max integer index of the last cell in the x dimension (right)
     * @param iy_min integer index of the first cell in the y dimension (bottom)
     * @param iy_max integer index of the last cell in the y dimension (top)
     * @param it_min integer index of the first cell in the datetime dimension (first)
     * @param it_max integer index of the last cell in the datetime dimension (last)
     * @return a shared pointer to the created data cube instance
     */
    static std::shared_ptr<crop_cube> create(std::shared_ptr<cube> in, int32_t ix_min, int32_t ix_max,
                                             int32_t iy_min, int32_t iy_max, int32_t it_min, int32_t it_max) {
        std::shared_ptr<crop_cube> out = std::make_shared<crop_cube>(in, ix_min, ix_max, iy_min, iy_max, it_min, it_max);
        in->add_child_cube(out);
        out->add_parent_cube(in);
        return out;
    }

   public:
    /**
     * @brief Create a data cube crops a source data cube by spatial and/or temporal ranges
     * @note This static creation method should preferably be used instead of the constructors as
     * the constructors will not set connections between cubes properly.
     * @param in input data cube
     * @param left
     * @param right
     * @param bottom
     * @param top
     * @param t0
     * @param t1
     * @param snap
     * @return a shared pointer to the created data cube instance
     */
    static std::shared_ptr<crop_cube> create(std::shared_ptr<cube> in, double left, double right,
                                             double bottom, double top, std::string t0, std::string t1,
                                             std::string snap = "near") {
        int32_t ix_min = 0;
        int32_t ix_max = 0;
        int32_t iy_min = 0;
        int32_t iy_max = 0;
        int32_t it_min = 0;
        int32_t it_max = 0;

        double x_min = (left - in->st_reference()->left()) / (double)in->st_reference()->dx();
        double x_max = -1 + (right - in->st_reference()->left()) / (double)in->st_reference()->dx();
        // double y_min = (bottom - in->st_reference()->bottom()) / (double)in->st_reference()->dy();
        // double y_max = -1 + (top - in->st_reference()->bottom()) / (double)in->st_reference()->dy();
        double y_min = (in->st_reference()->top() - top) / (double)in->st_reference()->dy();
        double y_max = -1 + (in->st_reference()->top() - bottom) / (double)in->st_reference()->dy();


        if (snap == "near") {
            ix_min = std::round(x_min);
            ix_max = std::round(x_max);
            iy_min = std::round(y_min);
            iy_max = std::round(y_max);
        }
        else if (snap == "in") {
            ix_min = std::ceil(x_min);
            ix_max = std::floor(x_max);
            iy_min = std::ceil(y_min);
            iy_max = std::floor(y_max);
        }
        else if (snap == "out") {
            ix_min = std::floor(x_min);
            ix_max = std::ceil(x_max);
            iy_min = std::floor(y_min);
            iy_max = std::ceil(y_max);
        }
        else {
            std::string msg = "Invalid argument: snap must be one of \"near\", \"in\", and \"out\"";
            GCBS_ERROR(msg);
            throw std::string(msg);
        }

        datetime dt0 = datetime::from_string(t0);
        datetime dt1 = datetime::from_string(t1);
        if (dt0.unit() != in->st_reference()->dt_unit() || dt1.unit() != in->st_reference()->dt_unit()) {
            std::string msg = "Invalid datetime unit of provided temporal extent, please use the datetime unit of the cube";
            GCBS_ERROR(msg);
            throw std::string(msg);
        }
        if (cube_stref::type_string(in->st_reference()) == "cube_stref_regular") {
            auto stref = std::dynamic_pointer_cast<cube_stref_regular>(in->st_reference());
            duration delta1 = (dt0 - stref->t0());
            duration delta2 = (dt1 - stref->t0());
            double t_min = (double)(delta1.dt_interval) / (double)(stref->dt().dt_interval);
            double t_max = (double)(delta2.dt_interval) / (double)(stref->dt().dt_interval);
            if (snap == "near") {
                it_min = std::round(t_min);
                it_max = std::round(t_max);
            }
            else if (snap == "in") {
                it_min = std::ceil(t_min);
                it_max = std::floor(t_max);
            }
            else if (snap == "out") {
                it_min = std::floor(t_min);
                it_max = std::ceil(t_max);
            }
        }
        else if (cube_stref::type_string(in->st_reference()) == "cube_stref_labeled_time") {
            auto stref = std::dynamic_pointer_cast<cube_stref_labeled_time>(in->st_reference());

            std::vector<datetime> labels = stref->get_time_labels();
            // Assuming labels are sorted

            if (dt0 < labels[0] || dt0 > labels[labels.size() - 1] ||
                dt1 < labels[0] || dt1 > labels[labels.size() - 1]) {
                std::string msg = "Temporal range of crop region reaches beyond data cube bounds";
                GCBS_ERROR(msg);
                throw std::string(msg);
            }

            int32_t i = 0;
            while (i < (int32_t)labels.size()) {
                if (labels[i] < dt0) {
                    it_min = i;
                }
                if (labels[i] < dt1) {
                    it_max = i;
                }
                ++i;
            }

            if (snap == "near") {
                // TODO: do we need to check for out-of-bounds cases?
                if ((dt0 - stref->datetime_at_index(it_min)).dt_interval >=
                    (stref->datetime_at_index(it_min + 1) - dt0).dt_interval) {
                    it_min++;
                }
                if ((dt1 - stref->datetime_at_index(it_max)).dt_interval >=
                    (stref->datetime_at_index(it_max + 1) - dt1).dt_interval) {
                    it_max++;
                }
            }
            else if (snap == "in") {
                it_min++;  // TODO: do we need to check for out-of-bounds cases?
            }
            else if (snap == "out") {
                it_max++;  // TODO: do we need to check for out-of-bounds cases?
            }
        }

        std::shared_ptr<crop_cube> out = std::make_shared<crop_cube>(in, ix_min, ix_max, iy_min, iy_max, it_min, it_max);
        in->add_child_cube(out);
        out->add_parent_cube(in);
        return out;
    }

   public:
    crop_cube(std::shared_ptr<cube> in, int32_t ix_min, int32_t ix_max, int32_t iy_min, int32_t iy_max,
              int32_t it_min, int32_t it_max) : cube(in->st_reference()->copy()),
                                                _in_cube(in),
                                                _x_min(ix_min),
                                                _x_max(ix_max),
                                                _y_min(iy_min),
                                                _y_max(iy_max),
                                                _t_min(it_min),
                                                _t_max(it_max) {
        _chunk_size[0] = _in_cube->chunk_size()[0];
        _chunk_size[1] = _in_cube->chunk_size()[1];
        _chunk_size[2] = _in_cube->chunk_size()[2];



        if (_x_max < 0 || _x_min >= (int32_t)in->size_x() ||
            _y_max < 0 || _y_min >= (int32_t)in->size_y() ||
            _t_max < 0 || _t_min >= (int32_t)in->size_t()) {
            std::string msg = "Crop region is completely outside the data cube extent";
            GCBS_ERROR(msg);
            throw std::string(msg);
        }

        if (_x_min < 0 || _x_min >= (int32_t)in->size_x() ||
            _x_max < 0 || _x_max >= (int32_t)in->size_x() ||
            _y_min < 0 || _y_min >= (int32_t)in->size_y() ||
            _y_max < 0 || _y_max >= (int32_t)in->size_y() ||
            _t_min < 0 || _t_min >= (int32_t)in->size_t() ||
            _t_max < 0 || _t_max >= (int32_t)in->size_t()) {
            std::string msg = "Crop region is partially outside the data cube extent";
            GCBS_DEBUG(msg);
        }

        // assert that min < max
        if (_x_min > _x_max) {
            std::swap(_x_min, _x_max);
        }
        if (_y_min > _y_max) {
            std::swap(_y_min, _y_max);
        }
        if (_t_min > _t_max) {
            std::swap(_t_min, _t_max);
        }

        auto stref = std::dynamic_pointer_cast<cube_stref_regular>(_st_ref);  // Notice that this works for labeled time axis too, because labeled st ref inherits from regular
        stref->set_x_axis(in->st_reference()->left() + _x_min * in->st_reference()->dx(),
                          in->st_reference()->left() + (_x_max + 1) * in->st_reference()->dx(),
                          (uint32_t)(_x_max - _x_min + 1));
        stref->set_y_axis(in->st_reference()->top() - (_y_max + 1) * in->st_reference()->dy(),
                          in->st_reference()->top() - (_y_min) * in->st_reference()->dy(),
                          (uint32_t)(_y_max - _y_min + 1));

        if (cube_stref::type_string(in->st_reference()) == "cube_stref_regular") {
            stref->set_t_axis(_in_cube->st_reference()->datetime_at_index(_t_min),
                              _in_cube->st_reference()->datetime_at_index(_t_max),
                              _in_cube->st_reference()->dt());
        }
        else if (cube_stref::type_string(in->st_reference()) == "cube_stref_labeled_time") {
            std::vector<datetime> time_labels = std::dynamic_pointer_cast<cube_stref_labeled_time>(_st_ref)->get_time_labels();
            std::vector<datetime> time_labels_cropped(&time_labels[_t_min], &time_labels[_t_max + 1]);
            std::dynamic_pointer_cast<cube_stref_labeled_time>(_st_ref)->set_time_labels(time_labels_cropped);
        }

        for (uint16_t i = 0; i < _in_cube->bands().count(); ++i) {
            _bands.add(in->bands().get(i));
        }
    }

   public:
    ~crop_cube() {}
    std::shared_ptr<chunk_data> read_chunk(chunkid_t id) override;
    json11::Json make_constructible_json() override {
        json11::Json::object out;
        out["cube_type"] = "crop";
        out["ix_min"] = _x_min;
        out["ix_max"] = _x_max;
        out["iy_min"] = _y_min;
        out["iy_max"] = _y_max;
        out["it_min"] = _t_min;
        out["it_max"] = _t_max;
        out["in_cube"] = _in_cube->make_constructible_json();
        return out;
    }

   private:
    std::shared_ptr<cube> _in_cube;
    int32_t _x_min;
    int32_t _x_max;
    int32_t _y_min;
    int32_t _y_max;
    int32_t _t_min;
    int32_t _t_max;
};

}  // namespace gdalcubes

#endif  // CROP_H
