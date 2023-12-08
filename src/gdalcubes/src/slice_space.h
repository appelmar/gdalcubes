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

#ifndef SLICE_SPACE_H
#define SLICE_SPACE_H

#include "cube.h"

    namespace gdalcubes {

/**
* @brief A data cube that extracts a single time series (spatial slice) from a source data cube
*/
class slice_space_cube : public cube {
   public:
    /**
 * @brief Create a data cube that extracts a time slice from a cube
 * @note This static creation method should preferably be used instead of the constructors as
 * the constructors will not set connections between cubes properly.
 * @param in input data cube
 * @param ix integer x coordinate of the requested slice
 * @param iy integer y coordinate of the requested slice
 * @return a shared pointer to the created data cube instance
     */
    static std::shared_ptr<slice_space_cube> create(std::shared_ptr<cube> in, int32_t ix, int32_t iy) {
        std::shared_ptr<slice_space_cube> out = std::make_shared<slice_space_cube>(in, ix, iy);
        in->add_child_cube(out);
        out->add_parent_cube(in);
        return out;
    }

    /**
 * @brief Create a data cube that extracts a time slice from a cube
 * @note This static creation method should preferably be used instead of the constructors as
 * the constructors will not set connections between cubes properly.
 * @param in input data cube
 * @param x spatial x coordinate of the requested slice, expected to be in the cube's crs
 * @param y spatial y coordinate of the requested slice, expected to be in the cube's crs
 * @return a shared pointer to the created data cube instance
     */
    static std::shared_ptr<slice_space_cube> create(std::shared_ptr<cube> in, double x, double y) {
        std::shared_ptr<slice_space_cube> out = std::make_shared<slice_space_cube>(in, x, y);
        in->add_child_cube(out);
        out->add_parent_cube(in);
        return out;
    }

   public:
    slice_space_cube(std::shared_ptr<cube> in, int32_t ix, int32_t iy) : cube(in->st_reference()->copy()), _in_cube(in), _x_index(ix), _y_index(iy) {
        _chunk_size[0] = _in_cube->chunk_size()[0];
        _chunk_size[1] = 1;
        _chunk_size[2] = 1;

        if (_x_index < 0 || _x_index >= (int32_t)in->size_x() ||
            _y_index < 0 || _y_index >= (int32_t)in->size_y()) {
            std::string msg = "Cell (x=" + std::to_string(_x_index) + ",y=" + std::to_string(_y_index) + ") is out of data cube bounds";
            GCBS_ERROR(msg);
            throw std::string(msg);
        }


        auto stref = std::dynamic_pointer_cast<cube_stref_regular>(_st_ref); // Notice that this works for labeled time axis too, because labeled st ref inherits from regular

        stref->set_x_axis(in->st_reference()->left() + _x_index * in->st_reference()->dx(),
                          in->st_reference()->left() + (_x_index + 1) * in->st_reference()->dx(),
                          uint32_t(1));

        stref->set_y_axis(in->st_reference()->top() - _y_index * in->st_reference()->dy(),
                          in->st_reference()->top() - (_y_index+1) * in->st_reference()->dy(),
                          uint32_t(1));

        if (_st_ref->nx() != 1) {
            std::string msg = "Data cube slice has invalid geometry: nx is not equal to 1";
            GCBS_ERROR(msg);
            throw std::string(msg);
        }
        if (_st_ref->ny() != 1) {
            std::string msg = "Data cube slice has invalid geometry: ny is not equal to 1";
            GCBS_ERROR(msg);
            throw std::string(msg);
        }
        for (uint16_t i = 0; i < _in_cube->bands().count(); ++i) {
            _bands.add(in->bands().get(i));
        }
    }

    slice_space_cube(std::shared_ptr<cube> in, double x, double y) : cube(in->st_reference()->copy()), _in_cube(in), _x_index(-1), _y_index(-1) {
        _x_index = (x - in->st_reference()->left()) / in->st_reference()->dx();
        _y_index = (in->st_reference()->top() - y) / in->st_reference()->dy();

        _chunk_size[0] = _in_cube->chunk_size()[0];
        _chunk_size[1] = 1;
        _chunk_size[2] = 1;

        if (_x_index < 0 || _x_index >= (int32_t)in->size_x() ||
            _y_index < 0 || _y_index >= (int32_t)in->size_y()) {
            std::string msg = "Cell (x=" + std::to_string(_x_index) + ",y=" + std::to_string(_y_index) + ") is out of data cube bounds";
            GCBS_ERROR(msg);
            throw std::string(msg);
        }


        auto stref = std::dynamic_pointer_cast<cube_stref_regular>(_st_ref); // Notice that this works for labeled time axis too, because labeled st ref inherits from regular

        stref->set_x_axis(in->st_reference()->left() + _x_index * in->st_reference()->dx(),
                          in->st_reference()->left() + (_x_index + 1) * in->st_reference()->dx(),
                          (uint32_t)1);
        stref->set_y_axis(in->st_reference()->bottom() + _y_index * in->st_reference()->dy(),
                          in->st_reference()->bottom() + (_y_index+1) * in->st_reference()->dy(),
                          (uint32_t)1);

        if (_st_ref->nx() != 1) {
            std::string msg = "Data cube slice has invalid geometry: nx is not equal to 1";
            GCBS_ERROR(msg);
            throw std::string(msg);
        }
        if (_st_ref->ny() != 1) {
            std::string msg = "Data cube slice has invalid geometry: ny is not equal to 1";
            GCBS_ERROR(msg);
            throw std::string(msg);
        }
        for (uint16_t i = 0; i < _in_cube->bands().count(); ++i) {
            _bands.add(in->bands().get(i));
        }
    }

   public:
    ~slice_space_cube() {}
    std::shared_ptr<chunk_data> read_chunk(chunkid_t id) override;
    json11::Json make_constructible_json() override {
        json11::Json::object out;
        out["cube_type"] = "slice_space";
        out["ix"] = _x_index;
        out["iy"] = _y_index;
        out["in_cube"] = _in_cube->make_constructible_json();
        return out;
    }

   private:
    std::shared_ptr<cube> _in_cube;
    int32_t _x_index;
    int32_t _y_index;

};

}  // namespace gdalcubes

#endif  //SLICE_SPACE_H
