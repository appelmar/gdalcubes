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

#ifndef SLICE_TIME_H
#define SLICE_TIME_H

#include "cube.h"

    namespace gdalcubes {

/**
 * @brief A data cube that extracts a temporal slice from a source data cube
 */
class slice_time_cube : public cube {
   public:
    /**
 * @brief Create a data cube that extracts a time slice from a cube
 * @note This static creation method should preferably be used instead of the constructors as
 * the constructors will not set connections between cubes properly.
 * @param in input data cube
 * @param t datetime string of the temporal slice
 * @return a shared pointer to the created data cube instance
     */
    static std::shared_ptr<slice_time_cube> create(std::shared_ptr<cube> in, std::string t) {
        std::shared_ptr<slice_time_cube> out = std::make_shared<slice_time_cube>(in, t);
        in->add_child_cube(out);
        out->add_parent_cube(in);
        return out;
    }

    /**
 * @brief Create a data cube that extracts a time slice from a cube
 * @note This static creation method should preferably be used instead of the constructors as
 * the constructors will not set connections between cubes properly.
 * @param in input data cube
 * @param t integer index of the temporal slice
 * @return a shared pointer to the created data cube instance
     */
    static std::shared_ptr<slice_time_cube> create(std::shared_ptr<cube> in, int32_t t) {
        std::shared_ptr<slice_time_cube> out = std::make_shared<slice_time_cube>(in, t);
        in->add_child_cube(out);
        out->add_parent_cube(in);
        return out;
    }

   public:
    slice_time_cube(std::shared_ptr<cube> in, int32_t t) : cube(in->st_reference()->copy()), _in_cube(in), _t_index(t) {
        _chunk_size[0] = 1;
        _chunk_size[1] = _in_cube->chunk_size()[1];
        _chunk_size[2] = _in_cube->chunk_size()[2];

        if (_t_index < 0 || _t_index >= (int32_t)in->size_t()) {
            GCBS_ERROR("Datetime is out of data cube bounds");
            throw std::string("Datetime is out of data cube bounds");
        }

        if (cube_stref::type_string(in->st_reference()) == "cube_stref_regular") {
            auto stref = std::dynamic_pointer_cast<cube_stref_regular>(_st_ref);
            stref->set_t_axis(_in_cube->st_reference()->datetime_at_index(_t_index),
                              _in_cube->st_reference()->datetime_at_index(_t_index),
                              in->st_reference()->dt());
        }
        else if (cube_stref::type_string(in->st_reference()) == "cube_stref_labeled_time") {
            auto stref = std::dynamic_pointer_cast<cube_stref_labeled_time>(_st_ref);
            stref->set_time_labels({_in_cube->st_reference()->datetime_at_index(_t_index)});
        }
        for (uint16_t i = 0; i < _in_cube->bands().count(); ++i) {
            _bands.add(in->bands().get(i));
        }
    }

    slice_time_cube(std::shared_ptr<cube> in, std::string t) : cube(in->st_reference()->copy()), _in_cube(in), _t_index(-1) {
        _chunk_size[0] = 1;
        _chunk_size[1] = _in_cube->chunk_size()[1];
        _chunk_size[2] = _in_cube->chunk_size()[2];

        datetime dt = datetime::from_string(t);

        _t_index = in->st_reference()->index_at_datetime(dt);
        if (_t_index < 0 || _t_index >= (int32_t)in->size_t()) {
            GCBS_ERROR("Datetime is out of data cube bounds");
            throw std::string("Datetime is out of data cube bounds");
        }
        if (cube_stref::type_string(in->st_reference()) == "cube_stref_regular") {
            auto stref = std::dynamic_pointer_cast<cube_stref_regular>(_st_ref);
            stref->set_t_axis(_in_cube->st_reference()->datetime_at_index(_t_index),
                              _in_cube->st_reference()->datetime_at_index(_t_index),
                              in->st_reference()->dt());
        }
        else if (cube_stref::type_string(in->st_reference()) == "cube_stref_labeled_time") {
            auto stref = std::dynamic_pointer_cast<cube_stref_labeled_time>(_st_ref);
            stref->set_time_labels({_in_cube->st_reference()->datetime_at_index(_t_index)});
        }
        for (uint16_t i = 0; i < _in_cube->bands().count(); ++i) {
            _bands.add(in->bands().get(i));
        }
    }

   public:
    ~slice_time_cube() {}
    std::shared_ptr<chunk_data> read_chunk(chunkid_t id) override;
    json11::Json make_constructible_json() override {
        json11::Json::object out;
        out["cube_type"] = "slice_time";
        out["t"] = _t_index;
        out["in_cube"] = _in_cube->make_constructible_json();
        return out;
    }

   private:
    std::shared_ptr<cube> _in_cube;
    int32_t _t_index;

};

}  // namespace gdalcubes

#endif  //SLICE_TIME_H
