/*
    MIT License

    Copyright (c) 2024 Marius Appel <marius.appel@hs-bochum.de>

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

#ifndef WINDOW_SPACE_H
#define WINDOW_SPACE_H

#include "cube.h"

namespace gdalcubes {






// TODO: add custom reducer function
// TODO: add boundary fill options (e.g. NA, wrap, repeat, mirror, ...)
/**
 * @brief A data cube that applies a focal operation / kernel over sliding spatial windows
 */
class window_space_cube : public cube {
   public:


    struct padding {
        enum MODE {NONE, CONSTANT, REPLICATE, REFLECT, REFLECT_PIXEL, WRAP}; // see https://processes.openeo.org/#apply_kernel 
        MODE mode;
        double constant_value;
        padding() : mode(NONE), constant_value(NAN) {}
    };
   
    /**
     * @brief Create a data cube that applies a reducer function on a given input data cube over time
     * @note This static creation method should preferably be used instead of the constructors as
     * the constructors will not set connections between cubes properly.
     * @param in input data cube
     * @param reducer reducer function
     * @return a shared pointer to the created data cube instance
     */
    static std::shared_ptr<window_space_cube>
    create(std::shared_ptr<cube> in, std::vector<std::pair<std::string, std::string>> reducer_bands,
           uint16_t win_size_y, uint16_t win_size_x, bool keep_bands, std::string pad_str, double pad_fill=0.0) {
        std::shared_ptr<window_space_cube> out = std::make_shared<window_space_cube>(in, reducer_bands, win_size_y, win_size_x, keep_bands, pad_str, pad_fill);
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
    static std::shared_ptr<window_space_cube>
    create(std::shared_ptr<cube> in, std::vector<double> kernel, uint16_t win_size_y, uint16_t win_size_x, bool keep_bands, std::string pad_str, double pad_fill=0.0) {
        std::shared_ptr<window_space_cube> out = std::make_shared<window_space_cube>(in, kernel, win_size_y, win_size_x, keep_bands, pad_str, pad_fill);
        in->add_child_cube(out);
        out->add_parent_cube(in);
        return out;
    }

   public:

   window_space_cube(std::shared_ptr<cube> in, std::vector<std::pair<std::string, std::string>> reducer_bands,
                     uint16_t win_size_y, uint16_t win_size_x, bool keep_bands, std::string pad_str, double pad_fill=0.0) : cube(in->st_reference()->copy()), _in_cube(in), _reducer_bands(reducer_bands), _win_size_y(win_size_y), _win_size_x(win_size_x), _band_idx_in(), _kernel(), _keep_bands(keep_bands), _pad_str(pad_str), _pad_fill(pad_fill), _pad() {  // it is important to duplicate st reference here, otherwise changes will affect input cube as well
        _chunk_size[0] = _in_cube->chunk_size()[0];
        _chunk_size[1] = _in_cube->chunk_size()[1];
        _chunk_size[2] = _in_cube->chunk_size()[2];

        if (win_size_y % 2 == 0 || win_size_x % 2 == 0) {
            GCBS_ERROR("Window size must not be even");
            throw std::string(
                "ERROR in window_space_cube::window_space_cube(): Window size must not be even");
        }

        if (_keep_bands) {
            // TODO: check for band name conflicts here
            for (uint16_t i = 0; i < _in_cube->size_bands(); ++i) {
                _bands.add(_in_cube->bands().get(i));
            }
        }
        for (uint16_t i = 0; i < reducer_bands.size(); ++i) {
            std::string reducerstr = reducer_bands[i].first;
            std::string bandstr = reducer_bands[i].second;
     
            band b = in->bands().get(bandstr);
            b.name = b.name + "_" + reducerstr;
            _bands.add(b);

            _band_idx_in.push_back(in->bands().get_index(bandstr));
        }

        if (pad_str == "CONSTANT") {
            _pad.mode = padding::MODE::CONSTANT;
            _pad.constant_value = pad_fill;
        }
        else if (pad_str == "REPLICATE") {
            _pad.mode = padding::MODE::REPLICATE;
        }
        else if (pad_str == "REFLECT") {
            _pad.mode = padding::MODE::REFLECT;
        }
        else if (pad_str == "REFLECT_PIXEL") {
            _pad.mode = padding::MODE::REFLECT_PIXEL;
        }
        else {
            if (!pad_str.empty())
                GCBS_WARN("Unknown padding method defined: falling back to default method (no padding)");
            _pad.mode = padding::MODE::NONE;
        }
    }

    
    window_space_cube(std::shared_ptr<cube> in, std::vector<double> kernel, uint16_t win_size_y, uint16_t win_size_x, bool keep_bands, std::string pad_str, double pad_fill=0.0)
        : cube(in->st_reference()->copy()), _in_cube(in), _reducer_bands(), _win_size_y(win_size_y), _win_size_x(win_size_x), _band_idx_in(), _kernel(kernel), _keep_bands(keep_bands), _pad_str(pad_str), _pad_fill(pad_fill), _pad() {  // it is important to duplicate st reference here, otherwise changes will affect input cube as well
        _chunk_size[0] = _in_cube->chunk_size()[0];
        _chunk_size[1] = _in_cube->chunk_size()[1];
        _chunk_size[2] = _in_cube->chunk_size()[2];

        if (win_size_y % 2 == 0 || win_size_x % 2 == 0) {
            GCBS_ERROR("Window size must not be even");
            throw std::string(
                "ERROR in window_space_cube::window_space_cube(): Window size must not be even");
        }

        if (kernel.size() != (win_size_y * win_size_x)) {
            GCBS_ERROR("kernel size does not match window size");
            throw std::string(
                "ERROR in window_space_cube::window_space_cube(): Kernel size does not match window size");
        }
        // Important: keep_bands is ignored if a kernel is used
        for (uint16_t i = 0; i < in->bands().count(); ++i) {
            band b = in->bands().get(i);
            _bands.add(b);
        }

        if (pad_str == "CONSTANT") {
            _pad.mode = padding::MODE::CONSTANT;
            _pad.constant_value = pad_fill;
        }
        else if (pad_str == "REPLICATE") {
            _pad.mode = padding::MODE::REPLICATE;
        }
        else if (pad_str == "REFLECT") {
            _pad.mode = padding::MODE::REFLECT;
        }
        else if (pad_str == "REFLECT_PIXEL") {
            _pad.mode = padding::MODE::REFLECT_PIXEL;
        }
        else {
            if (!pad_str.empty())
                GCBS_WARN("Unknown padding method defined: falling back to default method (no padding)");
            _pad.mode = padding::MODE::NONE;
        }
    }

   public:
    ~window_space_cube() {}

    std::shared_ptr<chunk_data> read_chunk(chunkid_t id) override;

    json11::Json make_constructible_json() override {
        json11::Json::object out;
        out["cube_type"] = "window_space";
        if (!_kernel.empty()) {
            out["kernel"] = _kernel;
        }
        else {
            json11::Json::array rb;
            for (uint16_t i = 0; i < _reducer_bands.size(); ++i) {
                rb.push_back(json11::Json::array({_reducer_bands[i].first, _reducer_bands[i].second}));
            }
            out["reducer_bands"] = rb;
        }
        out["win_size_y"] = _win_size_y;
        out["win_size_x"] = _win_size_x;
        out["pad_str"] = _pad_str;
        out["pad_fill"] = _pad_fill;
        out["keep_bands"] = _keep_bands;
        out["in_cube"] = _in_cube->make_constructible_json();
        return out;
    }

   private:
    std::shared_ptr<cube> _in_cube;
    std::vector<std::pair<std::string, std::string>> _reducer_bands;
    uint16_t _win_size_y; 
    uint16_t _win_size_x;
    std::vector<uint16_t> _band_idx_in;
    std::vector<double> _kernel;
    bool _keep_bands;
    std::string _pad_str;
    double _pad_fill;
    padding _pad; // TODO
};

}  // namespace gdalcubes

#endif  // WINDOW_SPACE_H
