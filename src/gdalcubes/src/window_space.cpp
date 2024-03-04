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
#include "window_space.h"

namespace gdalcubes {





struct window_reducer_singleband {
    virtual ~window_reducer_singleband() {}
    virtual void init() = 0;
    virtual void update(double &v) = 0;
    virtual double finalize() = 0;
};

struct window_reducer_mean : public window_reducer_singleband {
    void init() override {
        sum = 0.0;
        count = 0;
    }
    void update(double &v) override {
        sum += v;
        count++;
    }
    double finalize() override {
        return sum / double(count);
    }
    double sum;
    uint32_t count;
};


struct window_reducer_sum : public window_reducer_singleband {
    void init() override {
        sum = 0.0;
    }
    void update(double &v) override {
        sum += v;
    }
    double finalize() override {
        return sum;
    }
    double sum;
};

struct window_reducer_prod : public window_reducer_singleband {
    void init() override {
        prod = 0.0;
    }
    void update(double &v) override {
        prod *= v;
    }
    double finalize() override {
        return prod;
    }
    double prod;
};

struct window_reducer_max : public window_reducer_singleband {
    void init() override {
        max = NAN;
    }
    void update(double &v) override {
        if (std::isnan(max) || v > max) {
            max = v;
        }
    }
    double finalize() override {
        return max;
    }
    double max;
};

struct window_reducer_min : public window_reducer_singleband {
    void init() override {
        min = NAN;
    }
    void update(double &v) override {
        if (std::isnan(min) || v < min) {
            min = v;
        }
    }
    double finalize() override {
        return min;
    }
    double min;
};


struct window_reducer_count : public window_reducer_singleband {
    void init() override {
        count = 0;
    }
    void update(double &v) override {
        if (std::isfinite(v)) {
            count++;
        }
    }
    double finalize() override {
        return count;
    }
    uint32_t count;
};

struct window_reducer_median : public window_reducer_singleband {
    void init() override {
        values.clear();
    }
    void update(double &v) override {
        if (std::isfinite(v)) {
            values.push_back(v);
        }
    }
    double finalize() override {
        if (values.size() == 0) {
            return NAN;
        } 
        std::sort(values.begin(), values.end());
        if (values.size() % 2 == 1) {
            return values[values.size() / 2];
        } else {
            return (values[values.size() / 2] + values[values.size() / 2 - 1]) / ((double)2);
        }
    }
    std::vector<double> values;
};



struct window_reducer_var : public window_reducer_singleband {
    void init() override {
       mean = 0.0;
       count = 0;
       state = 0.0;
    }
    void update(double &v) override {
        if (std::isfinite(v)) {
            ++count;
            double delta = v - mean;
            mean += delta / count;
            state += delta * (v - mean);
        }
    }
    double finalize() override {
        return count > 1 ? state / (count - 1) : NAN;
    }

    double mean;
    uint32_t count;
    double state;
};


struct window_reducer_sd : public window_reducer_var {
    double finalize() override {
        return count > 1 ? sqrt(state / (count - 1)) : NAN;
    }
};




std::shared_ptr<chunk_data> window_space_cube::read_chunk(chunkid_t id) {
    GCBS_TRACE("window_space_cube::read_chunk(" + std::to_string(id) + ")");
    std::shared_ptr<chunk_data> out = std::make_shared<chunk_data>();
    if (id >= count_chunks())
        return out;  // chunk is outside of the view, we don't need to read anything.

    coords_nd<uint32_t, 3> size_tyx = chunk_size(id);

    coords_nd<uint32_t, 4> size_btyx = {uint32_t(_bands.count()), size_tyx[0], size_tyx[1], size_tyx[2]};
    out->size(size_btyx);

    // Fill buffers accordingly
    out->buf(std::calloc(size_btyx[0] * size_btyx[1] * size_btyx[2] * size_btyx[3], sizeof(double)));
    double* begin = (double*)out->buf();
    double* end = ((double*)out->buf()) + size_btyx[0] * size_btyx[1] * size_btyx[2] * size_btyx[3];
    std::fill(begin, end, NAN);


    // 1. Extract "subwindow" from input chunks
    auto ccoords = chunk_coords_from_id(id);
    std::array<int32_t, 3> lower = {int32_t(ccoords[0] * _chunk_size[0]), 
                                    int32_t(ccoords[1] * _chunk_size[1] - ((_win_size_y - 1) / 2)), 
                                    int32_t(ccoords[2] * _chunk_size[2] - ((_win_size_x - 1) / 2))};
    std::array<int32_t, 3> upper = {int32_t(ccoords[0] * _chunk_size[0] + size_tyx[0] - 1), 
                                    int32_t(ccoords[1] * _chunk_size[1] + size_tyx[1] - 1 + ((_win_size_y - 1) / 2)), 
                                    int32_t(ccoords[2] * _chunk_size[2] + size_tyx[2] - 1 + ((_win_size_x - 1) / 2))};

    auto cwin = _in_cube->read_window(lower, upper);





    // 2. Apply padding (if needed)
    if (_pad.mode != padding::MODE::NONE) {

        double v = NAN;
        if (_pad.mode == padding::MODE::CONSTANT) {
            v = _pad.constant_value;
        }


        // lower y side
        if (lower[1] < 0) {
            int32_t wy = -lower[1];
            for (uint16_t ib=0; ib < cwin->size()[0]; ++ib) {
                for (uint16_t it=0; it < cwin->size()[1]; ++it) {
                    for (int32_t ix=0; ix < int32_t(cwin->size()[3]); ++ix) {
                        for (int32_t off_y=1; off_y <= wy; ++off_y) {
                            int32_t tx = ix;
                            int32_t ty = wy - off_y; 
                            int32_t sx = tx; 
                            int32_t sy = ty; 
                            if (_pad.mode == padding::MODE::REPLICATE) {
                                sx = tx;
                                sy = ty + off_y; 
                            }
                            else if (_pad.mode == padding::MODE::REFLECT) {
                                sx = tx;
                                sy = ty + off_y + off_y - 1; 
                            }
                            else if (_pad.mode == padding::MODE::REFLECT_PIXEL) {
                                sx = tx;
                                sy = ty + off_y + off_y; 
                            }
                            if (_pad.mode != padding::MODE::CONSTANT) {
                                v = ((double*)(cwin->buf()))[
                                    ib * cwin->size()[1] * cwin->size()[2] * cwin->size()[3] +
                                    it * cwin->size()[2] * cwin->size()[3] +
                                    sy * cwin->size()[3] +
                                    sx];
                            }
                            ((double*)(cwin->buf()))[
                                ib * cwin->size()[1] * cwin->size()[2] * cwin->size()[3] +
                                it * cwin->size()[2] * cwin->size()[3] +
                                ty * cwin->size()[3] +
                                tx ] = v;
                        }
                    }
                }  
            }
        }

        // upper y side
        if (upper[1] >= int32_t(size_y())) {
            int32_t wy = upper[1] - size_y() + 1;
            for (uint16_t ib=0; ib < cwin->size()[0]; ++ib) {
                for (uint16_t it=0; it < cwin->size()[1]; ++it) {
                    for (int32_t ix=0; ix < int32_t(cwin->size()[3]); ++ix) {
                        for (int32_t off_y=1; off_y <= wy; ++off_y) {
                            int32_t tx = ix;
                            int32_t ty = cwin->size()[2] - wy - 1 + off_y;  
                            int32_t sx = tx; 
                            int32_t sy = ty; 

                            if (_pad.mode == padding::MODE::REPLICATE) {
                                sx = tx;
                                sy = ty - off_y; 
                            }
                            else if (_pad.mode == padding::MODE::REFLECT) {
                                sx = tx;
                                sy = ty - off_y - off_y + 1; 
                            }
                            else if (_pad.mode == padding::MODE::REFLECT_PIXEL) {
                                sx = tx;
                                sy = ty - off_y - off_y;  
                            }
                            if (_pad.mode != padding::MODE::CONSTANT) {
                                v = ((double*)(cwin->buf()))[
                                    ib * cwin->size()[1] * cwin->size()[2] * cwin->size()[3] +
                                    it * cwin->size()[2] * cwin->size()[3] +
                                    sy * cwin->size()[3] +
                                    sx];
                            }
                            ((double*)(cwin->buf()))[
                                ib * cwin->size()[1] * cwin->size()[2] * cwin->size()[3] +
                                it * cwin->size()[2] * cwin->size()[3] +
                                ty * cwin->size()[3] +
                                tx ] = v;
                        }
                    }                   
                }
            }
        }
        
        

        // lower x side
        if (lower[2] < 0) {
            int32_t wx = -lower[2];
            for (uint16_t ib=0; ib < cwin->size()[0]; ++ib) {
                for (uint16_t it=0; it < cwin->size()[1]; ++it) {      
                    for (int32_t iy=0; iy < int32_t(cwin->size()[2]); ++iy) {
                        for (int32_t off_x=1; off_x <= wx; ++off_x) {
                            
                            int32_t tx = wx - off_x;
                            int32_t ty = iy; 
                            int32_t sx = tx; 
                            int32_t sy = ty; 
                            if (_pad.mode == padding::MODE::REPLICATE) {
                                sy = ty;
                                sx = tx + off_x;
                            }
                            else if (_pad.mode == padding::MODE::REFLECT) {
                                sx = tx + off_x + off_x - 1;
                                sy = ty; 
                            }
                            else if (_pad.mode == padding::MODE::REFLECT_PIXEL) {
                                sx = tx + off_x + off_x;
                                sy = ty;
                            }
                            if (_pad.mode != padding::MODE::CONSTANT) {
                                v = ((double*)(cwin->buf()))[
                                    ib * cwin->size()[1] * cwin->size()[2] * cwin->size()[3] +
                                    it * cwin->size()[2] * cwin->size()[3] +
                                    sy * cwin->size()[3] +
                                    sx];
                            }
                            ((double*)(cwin->buf()))[
                                ib * cwin->size()[1] * cwin->size()[2] * cwin->size()[3] +
                                it * cwin->size()[2] * cwin->size()[3] +
                                ty * cwin->size()[3] +
                                tx ] = v;
                        }
                    }
                }  
            }
        }

        // upper x side
        if (upper[2] >= int32_t(size_x())) {
            int32_t wx = upper[2] - size_x() + 1;
            for (uint16_t ib=0; ib < cwin->size()[0]; ++ib) {
                for (uint16_t it=0; it < cwin->size()[1]; ++it) {
                    for (int32_t iy=0; iy < int32_t(cwin->size()[2]); ++iy) {
                        for (int32_t off_x=1; off_x <= wx; ++off_x) {
                            int32_t tx = cwin->size()[3] - wx - 1 + off_x; 
                            int32_t ty = iy;  
                            int32_t sx = tx; 
                            int32_t sy = ty; 

                            if (_pad.mode == padding::MODE::REPLICATE) {
                                sx = tx - off_x;
                                sy = ty;
                            }
                            else if (_pad.mode == padding::MODE::REFLECT) {
                                sx = tx - off_x - off_x + 1;
                                sy = ty;
                            }
                            else if (_pad.mode == padding::MODE::REFLECT_PIXEL) {
                                sx = tx - off_x - off_x;
                                sy = ty;  
                            }
                            if (_pad.mode != padding::MODE::CONSTANT) {
                                v = ((double*)(cwin->buf()))[
                                    ib * cwin->size()[1] * cwin->size()[2] * cwin->size()[3] +
                                    it * cwin->size()[2] * cwin->size()[3] +
                                    sy * cwin->size()[3] +
                                    sx];
                            }
                            ((double*)(cwin->buf()))[
                                ib * cwin->size()[1] * cwin->size()[2] * cwin->size()[3] +
                                it * cwin->size()[2] * cwin->size()[3] +
                                ty * cwin->size()[3] +
                                tx ] = v;
                        }
                    }                      
                }
            }
        }

        // lower y, lower x side
        if (lower[1] < 0 && lower[2] < 0) {
            int32_t wx = -lower[2];
            int32_t wy = -lower[1];
            for (uint16_t ib=0; ib < cwin->size()[0]; ++ib) {
                for (uint16_t it=0; it < cwin->size()[1]; ++it) {
                    for (int32_t off_y=1; off_y <= wy; ++off_y) {
                        for (int32_t off_x=1; off_x <= wx; ++off_x) {
                            int32_t tx = wx - off_x;
                            int32_t ty = wy - off_y; 
                            int32_t sx = tx; 
                            int32_t sy = ty; 
                            if (_pad.mode == padding::MODE::REPLICATE) {
                                sx = tx + off_y;
                                sy = ty + off_y; 
                            }
                            else if (_pad.mode == padding::MODE::REFLECT) {
                                sx = tx + off_x + off_x - 1;
                                sy = ty + off_y + off_y - 1; 
                            }
                            else if (_pad.mode == padding::MODE::REFLECT_PIXEL) {
                                sx = tx + off_x + off_x;
                                sy = ty + off_y + off_y; 
                            }
                            if (_pad.mode != padding::MODE::CONSTANT) {
                                v = ((double*)(cwin->buf()))[
                                    ib * cwin->size()[1] * cwin->size()[2] * cwin->size()[3] +
                                    it * cwin->size()[2] * cwin->size()[3] +
                                    sy * cwin->size()[3] +
                                    sx];
                            }
                            ((double*)(cwin->buf()))[
                                ib * cwin->size()[1] * cwin->size()[2] * cwin->size()[3] +
                                it * cwin->size()[2] * cwin->size()[3] +
                                ty * cwin->size()[3] +
                                tx ] = v;
                        }
                    }                      
                }
            }    
        }
        // upper y, lower x side
        if (upper[1] >= int32_t(size_y()) && lower[2] < 0) {
            int32_t wx = -lower[2];
            int32_t wy = upper[1] - size_y() + 1;
            for (uint16_t ib=0; ib < cwin->size()[0]; ++ib) {
                for (uint16_t it=0; it < cwin->size()[1]; ++it) {
                    for (int32_t off_y=1; off_y <= wy; ++off_y) {
                        for (int32_t off_x=1; off_x <= wx; ++off_x) {
                            int32_t tx = wx - off_x;
                            int32_t ty = cwin->size()[2] - wy - 1 + off_y;
                            int32_t sx = tx; 
                            int32_t sy = ty; 
                            if (_pad.mode == padding::MODE::REPLICATE) {
                                sx = tx + off_y;
                                sy = ty - off_y; 
                            }
                            else if (_pad.mode == padding::MODE::REFLECT) {
                                sx = tx + off_x + off_x - 1;
                                sy = ty - off_y - off_y + 1; 
                            }
                            else if (_pad.mode == padding::MODE::REFLECT_PIXEL) {
                                sx = tx + off_x + off_x;
                                sy = ty - off_y - off_y;  
                            }
                            if (_pad.mode != padding::MODE::CONSTANT) {
                                v = ((double*)(cwin->buf()))[
                                    ib * cwin->size()[1] * cwin->size()[2] * cwin->size()[3] +
                                    it * cwin->size()[2] * cwin->size()[3] +
                                    sy * cwin->size()[3] +
                                    sx];
                            }
                            ((double*)(cwin->buf()))[
                                ib * cwin->size()[1] * cwin->size()[2] * cwin->size()[3] +
                                it * cwin->size()[2] * cwin->size()[3] +
                                ty * cwin->size()[3] +
                                tx ] = v;
                        }
                    }                      
                }
            }    
        }
        // upper y, upper x side
        if (upper[1] >= int32_t(size_y()) && upper[2] >= int32_t(size_x())) {
            int32_t wx = upper[2] - size_x() + 1;
            int32_t wy = upper[1] - size_y() + 1;
            for (uint16_t ib=0; ib < cwin->size()[0]; ++ib) {
                for (uint16_t it=0; it < cwin->size()[1]; ++it) {
                    for (int32_t off_y=1; off_y <= wy; ++off_y) {
                        for (int32_t off_x=1; off_x <= wx; ++off_x) {
                            int32_t tx = cwin->size()[3] - wx - 1 + off_x;
                            int32_t ty = cwin->size()[2] - wy - 1 + off_y;
                            int32_t sx = tx; 
                            int32_t sy = ty; 
                            if (_pad.mode == padding::MODE::REPLICATE) {
                                sx = tx - off_x;
                                sy = ty - off_y; 
                            }
                            else if (_pad.mode == padding::MODE::REFLECT) {
                                sx = tx - off_x - off_x + 1;
                                sy = ty - off_y - off_y + 1; 
                            }
                            else if (_pad.mode == padding::MODE::REFLECT_PIXEL) {
                                sx = tx - off_x - off_x;
                                sy = ty - off_y - off_y;  
                            }
                            if (_pad.mode != padding::MODE::CONSTANT) {
                                v = ((double*)(cwin->buf()))[
                                    ib * cwin->size()[1] * cwin->size()[2] * cwin->size()[3] +
                                    it * cwin->size()[2] * cwin->size()[3] +
                                    sy * cwin->size()[3] +
                                    sx];
                            }
                            ((double*)(cwin->buf()))[
                                ib * cwin->size()[1] * cwin->size()[2] * cwin->size()[3] +
                                it * cwin->size()[2] * cwin->size()[3] +
                                ty * cwin->size()[3] +
                                tx ] = v;
                        }
                    }                      
                }
            }    
        }
        // lower y, upper x side
        if (upper[1] >= int32_t(size_y()) && upper[2] >= int32_t(size_x())) {
            int32_t wx = upper[2] - size_x() + 1;
            int32_t wy = -lower[1];
            for (uint16_t ib=0; ib < cwin->size()[0]; ++ib) {
                for (uint16_t it=0; it < cwin->size()[1]; ++it) {
                    for (int32_t off_y=1; off_y <= wy; ++off_y) {
                        for (int32_t off_x=1; off_x <= wx; ++off_x) {
                            int32_t tx = cwin->size()[3] - wx - 1 + off_x;
                            int32_t ty = wy - off_y; 
                            int32_t sx = tx; 
                            int32_t sy = ty; 
                            if (_pad.mode == padding::MODE::REPLICATE) {
                                sx = tx - off_x;
                                sy = ty + off_y;  
                            }
                            else if (_pad.mode == padding::MODE::REFLECT) {
                                sx = tx - off_x - off_x + 1;
                                sy = ty + off_y + off_y - 1;
                            }
                            else if (_pad.mode == padding::MODE::REFLECT_PIXEL) {
                                sx = tx - off_x - off_x;
                                sy = ty + off_y + off_y;  
                            }
                            if (_pad.mode != padding::MODE::CONSTANT) {
                                v = ((double*)(cwin->buf()))[
                                    ib * cwin->size()[1] * cwin->size()[2] * cwin->size()[3] +
                                    it * cwin->size()[2] * cwin->size()[3] +
                                    sy * cwin->size()[3] +
                                    sx];
                            }
                            ((double*)(cwin->buf()))[
                                ib * cwin->size()[1] * cwin->size()[2] * cwin->size()[3] +
                                it * cwin->size()[2] * cwin->size()[3] +
                                ty * cwin->size()[3] +
                                tx ] = v;
                        }
                    }                      
                }
            }    
        }
    } // end padding















    // 3. Apply kernel or aggregation
    
    // CASE 1: Convolution using a provided kernel
    if (!_kernel.empty()) {
        
        // Iterate over target buffer and apply kernel.
        for (uint32_t ib = 0; ib < size_btyx[0]; ++ib) {
            for (uint32_t it = 0; it < size_btyx[1]; ++it) {
                for (uint32_t iy = 0; iy < size_btyx[2]; ++iy) {
                    for (uint32_t ix = 0; ix < size_btyx[3]; ++ix) {

                        // apply kernel
                        double sum = 0.0; // TODO: complete.cases only?
                        for (int16_t ky=0; ky < _win_size_y; ++ky) {
                            for (int16_t kx=0; kx < _win_size_x; ++kx) {
                                double v = ((double*)(cwin->buf()))[
                                        ib * cwin->size()[1] * cwin->size()[2] * cwin->size()[3] +
                                        it * cwin->size()[2] * cwin->size()[3] +
                                        (iy+ky) * cwin->size()[3] +
                                        (ix + kx)
                                    ];
                                if (std::isfinite(v)) { //
                                    sum += _kernel[ky * _win_size_x + kx] * v;
                                }
                                else {
                                    sum = NAN;
                                    break;
                                }                                  
                            }
                        }

                        ((double*)(out->buf()))[
                            ib * size_tyx[0] * size_tyx[1] * size_tyx[2] + 
                            it * size_tyx[1] * size_tyx[2] + 
                            iy * size_tyx[2] + 
                            ix] = sum;
                            
                    }
                }   
            }
        }
    }





    // CASE 2: Aggregation of selected bands using built-in aggregation functions (see top of file)
    else {
        std::vector<window_reducer_singleband *> reducers;
        for (uint16_t i = 0; i < _reducer_bands.size(); ++i) {
            window_reducer_singleband *r = nullptr;
            if (_reducer_bands[i].first == "min") {
                r = new window_reducer_min();
            } else if (_reducer_bands[i].first == "max") {
                r = new window_reducer_max();
            } else if (_reducer_bands[i].first == "mean") {
                r = new window_reducer_mean();
            } else if (_reducer_bands[i].first == "median") {
                r = new window_reducer_median();
            } else if (_reducer_bands[i].first == "sum") {
                r = new window_reducer_sum();
            } else if (_reducer_bands[i].first == "count") {
                r = new window_reducer_count();
            } else if (_reducer_bands[i].first == "prod") {
                r = new window_reducer_prod();
            } else if (_reducer_bands[i].first == "var") {
                r = new window_reducer_var();
            } else if (_reducer_bands[i].first == "sd") {
                r = new window_reducer_sd();
            } else
                throw std::string("ERROR in window_space_cube::read_chunk(): Unknown reducer given");

            reducers.push_back(r);
        }

        uint32_t ib=0;
        if (_keep_bands) {
            int32_t offst_y = (_win_size_y-1) / 2;
            int32_t offst_x = (_win_size_x-1) / 2;
            for (; ib < _in_cube->size_bands(); ++ib) {
                for (uint32_t it = 0; it < size_btyx[1]; ++it) {
                    for (uint32_t iy = 0; iy < size_btyx[2]; ++iy) {
                        for (uint32_t ix = 0; ix < size_btyx[3]; ++ix) {

                            double v = ((double*)(cwin->buf()))[
                                        ib * cwin->size()[1] * cwin->size()[2] * cwin->size()[3] +
                                        it * cwin->size()[2] * cwin->size()[3] +
                                        (iy + offst_y) * cwin->size()[3] +
                                        (ix + offst_x)
                                    ];

                            ((double*)(out->buf()))[
                                ib * size_tyx[0] * size_tyx[1] * size_tyx[2] + 
                                it * size_tyx[1] * size_tyx[2] + 
                                iy * size_tyx[2] + 
                                ix] = v;
                        }
                    }
                }
            }
        }
        int16_t offst_b = _keep_bands ? _in_cube->size_bands() : 0;
        for (; ib < size_btyx[0]; ++ib) {
            for (uint32_t it = 0; it < size_btyx[1]; ++it) {
                for (uint32_t iy = 0; iy < size_btyx[2]; ++iy) {
                    for (uint32_t ix = 0; ix < size_btyx[3]; ++ix) {
                        // apply reducer
                        reducers[ib-offst_b]->init();
                        for (int16_t ky=0; ky < _win_size_y; ++ky) {
                            for (int16_t kx=0; kx < _win_size_x; ++kx) {
                                double v = ((double*)(cwin->buf()))[
                                        _band_idx_in[ib-offst_b] * cwin->size()[1] * cwin->size()[2] * cwin->size()[3] +
                                        it * cwin->size()[2] * cwin->size()[3] +
                                        (iy+ky) * cwin->size()[3] +
                                        (ix + kx)
                                    ];
                                reducers[ib-offst_b]->update(v);    
                            }
                        }

                        ((double*)(out->buf()))[
                            ib * size_tyx[0] * size_tyx[1] * size_tyx[2] + 
                            it * size_tyx[1] * size_tyx[2] + 
                            iy * size_tyx[2] + 
                            ix] = reducers[ib-offst_b]->finalize();
                            
                    }
                }   
            }
        }
        for (uint16_t i = 0; i < reducers.size(); ++i) {
            if (reducers[i]) delete reducers[i];
        }
    }

    // check if chunk is completely NAN and if yes, return empty chunk
    if (out->all_nan()) {
        out = std::make_shared<chunk_data>();
    }
    return out;
}


}  // namespace gdalcubes
