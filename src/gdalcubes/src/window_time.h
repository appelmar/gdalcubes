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

#ifndef WINDOW_TIME_H
#define WINDOW_TIME_H

#include "cube.h"

namespace gdalcubes {

// TODO: add custom reducer function
// TODO: add boundary fill options (e.g. NA, wrap, repeat, mirror, ...)
/**
 * @brief A data cube that applies reducer functions over selected bands of a data cube over time
 * @note This is a reimplementation of reduce_cube. The new implementation allows to apply different reducers to different bands instead of just one reducer to all bands of the input data cube
 */
class window_time_cube : public cube {
   public:
    /**
         * @brief Create a data cube that applies a reducer function on a given input data cube over time
         * @note This static creation method should preferably be used instead of the constructors as
         * the constructors will not set connections between cubes properly.
         * @param in input data cube
         * @param reducer reducer function
         * @return a shared pointer to the created data cube instance
         */
    static std::shared_ptr<window_time_cube>
    create(std::shared_ptr<cube> in, std::vector<std::pair<std::string, std::string>> reducer_bands,
           uint16_t win_size_l, uint16_t win_size_r) {
        std::shared_ptr<window_time_cube> out = std::make_shared<window_time_cube>(in, reducer_bands, win_size_l,
                                                                                   win_size_r);
        in->add_child_cube(out);
        out->add_parent_cube(in);
        return out;
    }

    /**
        * @brief Create a data cube that applies a convolution kernel on a given input data cube over time
        * @note This static creation method should preferably be used instead of the constructors as
        * the constructors will not set connections between cubes properly.
        * @param in input data cube
        * @param reducer reducer function
        * @return a shared pointer to the created data cube instance
        */
    static std::shared_ptr<window_time_cube>
    create(std::shared_ptr<cube> in, std::vector<double> kernel, uint16_t win_size_l, uint16_t win_size_r) {
        std::shared_ptr<window_time_cube> out = std::make_shared<window_time_cube>(in, kernel, win_size_l,
                                                                                   win_size_r);
        in->add_child_cube(out);
        out->add_parent_cube(in);
        return out;
    }

   public:
    window_time_cube(std::shared_ptr<cube> in, std::vector<std::pair<std::string, std::string>> reducer_bands,
                     uint16_t win_size_l, uint16_t win_size_r) : cube(in->st_reference()->copy()), _in_cube(in), _reducer_bands(reducer_bands), _win_size_l(win_size_l), _win_size_r(win_size_r), _f(), _band_idx_in(), _kernel() {  // it is important to duplicate st reference here, otherwise changes will affect input cube as well
        _chunk_size[0] = _in_cube->chunk_size()[0];
        _chunk_size[1] = _in_cube->chunk_size()[1];
        _chunk_size[2] = _in_cube->chunk_size()[2];

        if (!_st_ref->has_regular_time()) {
            GCBS_WARN("Cube has irregular time dimension, window sizes may vary over time");
        }

        for (uint16_t i = 0; i < reducer_bands.size(); ++i) {
            std::string reducerstr = reducer_bands[i].first;
            std::string bandstr = reducer_bands[i].second;
            _f.push_back(get_default_reducer_by_name(reducerstr));

            band b = in->bands().get(bandstr);
            b.name = b.name + "_" + reducerstr;
            _bands.add(b);

            _band_idx_in.push_back(in->bands().get_index(bandstr));
        }
    }

    window_time_cube(std::shared_ptr<cube> in, std::vector<double> kernel, uint16_t win_size_l, uint16_t win_size_r)
        : cube(in->st_reference()->copy()), _in_cube(in), _reducer_bands(), _win_size_l(win_size_l), _win_size_r(win_size_r), _f(), _band_idx_in(), _kernel(kernel) {  // it is important to duplicate st reference here, otherwise changes will affect input cube as well
        _chunk_size[0] = _in_cube->chunk_size()[0];
        _chunk_size[1] = _in_cube->chunk_size()[1];
        _chunk_size[2] = _in_cube->chunk_size()[2];

        if (!win_size_l + (uint32_t)1 + win_size_r == kernel.size()) {
            GCBS_ERROR("kernel size does not match window size");
            throw std::string(
                "ERROR in window_time_cube::window_time_cube(): Kernel size does not match window size");
        }

        for (uint16_t i = 0; i < in->bands().count(); ++i) {
            _f.push_back(get_kernel_reducer(kernel));

            band b = in->bands().get(i);
            _bands.add(b);

            _band_idx_in.push_back(i);
        }
    }

   public:
    ~window_time_cube() {}

    std::shared_ptr<chunk_data> read_chunk(chunkid_t id) override;

    json11::Json make_constructible_json() override {
        json11::Json::object out;
        out["cube_type"] = "window_time";
        if (!_kernel.empty()) {
            out["kernel"] = _kernel;
        } else {
            json11::Json::array rb;
            for (uint16_t i = 0; i < _reducer_bands.size(); ++i) {
                rb.push_back(json11::Json::array({_reducer_bands[i].first, _reducer_bands[i].second}));
            }
            out["reducer_bands"] = rb;
        }
        out["win_size_l"] = _win_size_l;
        out["win_size_r"] = _win_size_r;
        out["in_cube"] = _in_cube->make_constructible_json();
        return out;
    }

   private:
    std::shared_ptr<cube> _in_cube;
    std::vector<std::pair<std::string, std::string>> _reducer_bands;
    uint16_t _win_size_l;
    uint16_t _win_size_r;
    std::vector<std::function<double(double *buf, uint16_t n)>> _f;
    std::vector<uint16_t> _band_idx_in;
    std::vector<double> _kernel;

    std::function<double(double *buf, uint16_t n)> get_default_reducer_by_name(std::string name);

    std::function<double(double *buf, uint16_t n)> get_kernel_reducer(std::vector<double> kernel);
};

}  // namespace gdalcubes

#endif  // WINDOW_TIME_H
