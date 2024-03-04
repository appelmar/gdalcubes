/*
    MIT License

    Copyright (c) 2022 Marius Appel <marius.appel@hs-bochum.de>

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


#include "aggregate_space.h"



struct aggregator_space_slice_singleband {
    virtual ~aggregator_space_slice_singleband() {}

    virtual void init(double *out, uint32_t size_t, uint32_t size_y, uint32_t size_x)= 0;
    virtual void combine(double* out, double* v, uint32_t it, uint32_t iy, uint32_t ix,
                         uint32_t size_t, uint32_t size_y, uint32_t size_x) = 0;
    virtual void finalize(double *out, uint32_t size_t, uint32_t size_y, uint32_t size_x)= 0;
};




/**
 * @brief Implementation of reducer to calculate mean values over time
 */
struct mean_aggregtor_space_slice_singleband : public aggregator_space_slice_singleband {
    void init(double *out, uint32_t size_t, uint32_t size_y, uint32_t size_x) override {

        _count = (uint32_t *)std::calloc(size_t * size_x * size_y, sizeof(uint32_t));
        for (uint32_t ixyt = 0; ixyt < size_t * size_x *size_y; ++ixyt) {
            _count[ixyt] = 0;
            out[ixyt] = 0;
        }
    }

    void combine(double* out, double* v, uint32_t it, uint32_t iy, uint32_t ix,
                         uint32_t size_t, uint32_t size_y, uint32_t size_x) override {
        if (!std::isnan(*v)) {
            out[it * size_y * size_x + iy * size_x + ix] += *v;
            ++_count[it * size_y * size_x + iy * size_x + ix];
        }
    }

    void finalize(double *out, uint32_t size_t, uint32_t size_y, uint32_t size_x) override {
        // divide by count;
        for (uint32_t ixyt = 0; ixyt < size_t * size_x * size_y; ++ixyt) {
            out[ixyt] =  _count[ixyt] > 0 ? out[ixyt] / _count[ixyt] : NAN;
        }
        std::free(_count);
    }

   private:
    uint32_t *_count;
};


struct min_aggregtor_space_slice_singleband : public aggregator_space_slice_singleband {
    void init(double *out, uint32_t size_t, uint32_t size_y, uint32_t size_x) override {
        for (uint32_t ixyt = 0; ixyt < size_t * size_x *size_y; ++ixyt) {
            out[ixyt] = NAN;
        }
    }
    void combine(double* out, double* v, uint32_t it, uint32_t iy, uint32_t ix,
                         uint32_t size_t, uint32_t size_y, uint32_t size_x) override {
        if (!std::isnan(*v)) {
            uint32_t ixyt = it * size_y * size_x + iy * size_x + ix;
            if (std::isnan(out[ixyt])) {
                out[ixyt] = *v;
            }
            else {
                out[ixyt] = std::min(out[ixyt], *v);
            }
        }
    }

    void finalize(double *out, uint32_t size_t, uint32_t size_y, uint32_t size_x) override {}
};

struct max_aggregtor_space_slice_singleband : public aggregator_space_slice_singleband {
    void init(double *out, uint32_t size_t, uint32_t size_y, uint32_t size_x) override {
        for (uint32_t ixyt = 0; ixyt < size_t * size_x *size_y; ++ixyt) {
            out[ixyt] = NAN;
        }
    }
    void combine(double* out, double* v, uint32_t it, uint32_t iy, uint32_t ix,
                         uint32_t size_t, uint32_t size_y, uint32_t size_x) override {
        if (!std::isnan(*v)) {
            uint32_t ixyt = it * size_y * size_x + iy * size_x + ix;
            if (std::isnan(out[ixyt])) {
                out[ixyt] = *v;
            }
            else {
                out[ixyt] = std::max(out[ixyt], *v);
            }
        }
    }

    void finalize(double *out, uint32_t size_t, uint32_t size_y, uint32_t size_x) override {}
};


struct sum_aggregtor_space_slice_singleband : public aggregator_space_slice_singleband {
    void init(double *out, uint32_t size_t, uint32_t size_y, uint32_t size_x) override {
        for (uint32_t ixyt = 0; ixyt < size_t * size_x *size_y; ++ixyt) {
            out[ixyt] = NAN;
        }
    }
    void combine(double* out, double* v, uint32_t it, uint32_t iy, uint32_t ix,
                         uint32_t size_t, uint32_t size_y, uint32_t size_x) override {

        if (!std::isnan(*v)) {
            uint32_t ixyt = it * size_y * size_x + iy * size_x + ix;
            if (std::isnan(out[ixyt])) {
                out[ixyt] = *v;
            }
            else {
                out[ixyt] += *v;
            }
        }
    }

    void finalize(double *out, uint32_t size_t, uint32_t size_y, uint32_t size_x) override {}
};



struct prod_aggregtor_space_slice_singleband : public aggregator_space_slice_singleband {
    void init(double *out, uint32_t size_t, uint32_t size_y, uint32_t size_x) override {
        for (uint32_t ixyt = 0; ixyt < size_t * size_x *size_y; ++ixyt) {
            out[ixyt] = NAN;
        }
    }
    void combine(double* out, double* v, uint32_t it, uint32_t iy, uint32_t ix,
                         uint32_t size_t, uint32_t size_y, uint32_t size_x) override {
        if (!std::isnan(*v)) {
            uint32_t ixyt = it * size_y * size_x + iy * size_x + ix;
            if (std::isnan(out[ixyt])) {
                out[ixyt] = *v;
            }
            else {
                out[ixyt] *= *v;
            }
        }
    }

    void finalize(double *out, uint32_t size_t, uint32_t size_y, uint32_t size_x) override {}
};



struct count_aggregtor_space_slice_singleband : public aggregator_space_slice_singleband {
    void init(double *out, uint32_t size_t, uint32_t size_y, uint32_t size_x) override {
        for (uint32_t ixyt = 0; ixyt < size_t * size_x *size_y; ++ixyt) {
            out[ixyt] = 0;
        }
    }
    void combine(double* out, double* v, uint32_t it, uint32_t iy, uint32_t ix,
                         uint32_t size_t, uint32_t size_y, uint32_t size_x) override {
        if (!std::isnan(*v)) {
            uint32_t ixyt = it * size_y * size_x + iy * size_x + ix;
            out[ixyt] += 1.0;
        }
    }

    void finalize(double *out, uint32_t size_t, uint32_t size_y, uint32_t size_x)override {}
};




struct median_aggregtor_space_slice_singleband : public aggregator_space_slice_singleband {
    void init(double *out, uint32_t size_t, uint32_t size_y, uint32_t size_x) override {
        _m_buckets.resize(size_t * size_x *size_y, std::vector<double>());
        for (uint32_t ixyt = 0; ixyt < size_t *size_x *size_y; ++ixyt) {
            out[ixyt] = NAN;
        }
    }

    void combine(double* out, double* v, uint32_t it, uint32_t iy, uint32_t ix,
                         uint32_t size_t, uint32_t size_y, uint32_t size_x) override {


        if (!std::isnan(*v)) {
            uint32_t ixyt = it * size_y * size_x + iy * size_x + ix;
            _m_buckets[ixyt].push_back(*v);
        }
    }

    void finalize(double *out, uint32_t size_t, uint32_t size_y, uint32_t size_x) override {
        for (uint32_t ixyt = 0; ixyt < size_t * size_x * size_y; ++ixyt) {
            std::vector<double> &list = _m_buckets[ixyt];
            std::sort(list.begin(), list.end());
            if (list.size() == 0) {
                out[ixyt] = NAN;
            } else if (list.size() % 2 == 1) {
                out[ixyt] = list[list.size() / 2];
            } else {
                out[ixyt] = (list[list.size() / 2] + list[list.size() / 2 - 1]) / ((double)2);
            }
        }
    }
   private:
    std::vector<std::vector<double>> _m_buckets;
};


/**
 * @brief Implementation of reducer to calculate variance values over time using Welford's Online algorithm
 */
struct var_aggregtor_space_slice_singleband : public aggregator_space_slice_singleband {
    void init(double *out, uint32_t size_t, uint32_t size_y, uint32_t size_x) override {
        _count = (uint32_t *)std::calloc(size_t * size_x * size_y, sizeof(uint32_t));
        _mean = (double *)std::calloc(size_t * size_x * size_y, sizeof(double));

        for (uint32_t ixyt = 0; ixyt < size_t * size_x * size_y; ++ixyt) {
            _count[ixyt] = 0;
            _mean[ixyt] = 0;
            out[ixyt] = 0;
        }
    }
    void combine(double* out, double* v, uint32_t it, uint32_t iy, uint32_t ix,
                         uint32_t size_t, uint32_t size_y, uint32_t size_x) override {

        if (!std::isnan(*v)) {
            uint32_t ixyt = it * size_y * size_x + iy * size_x + ix;
            double &mean = _mean[ixyt];
            uint32_t &count = _count[ixyt];
            ++count;
            double delta = *v - mean;
            mean += delta / count;
            out[ixyt] += delta * (*v - mean);
        }
    }

    void finalize(double *out, uint32_t size_t, uint32_t size_y, uint32_t size_x) override {
        for (uint32_t ixyt = 0; ixyt < size_t * size_x * size_y; ++ixyt) {
            out[ixyt] = _count[ixyt] > 1 ? out[ixyt] / (_count[ixyt] - 1) : NAN;
        }
        std::free(_count);
        std::free(_mean);
    }

   protected:
    uint32_t *_count;
    double *_mean;
};


struct sd_aggregtor_space_slice_singleband : public var_aggregtor_space_slice_singleband {
    void finalize(double *out, uint32_t size_t, uint32_t size_y, uint32_t size_x) override {
        for (uint32_t ixyt = 0; ixyt < size_t * size_x * size_y; ++ixyt) {
            out[ixyt] = _count[ixyt] > 1 ? std::sqrt(out[ixyt] / (_count[ixyt] - 1)) : NAN;
        }
        std::free(_count);
        std::free(_mean);
    }
};


namespace gdalcubes {

std::shared_ptr<chunk_data> aggregate_space_cube::read_chunk(chunkid_t id) {
    GCBS_TRACE("aggregate_space_cube::read_chunk(" + std::to_string(id) + ")");
    std::shared_ptr<chunk_data> out = std::make_shared<chunk_data>();
    if (id >= count_chunks())
        return out;  // chunk is outside of the view, we don't need to read anything.

    coords_nd<uint32_t, 3> size_tyx = chunk_size(id);
    coords_nd<uint32_t, 4> size_btyx = {uint32_t(_bands.count()), size_tyx[0], size_tyx[1], size_tyx[2]};
    out->size(size_btyx);
    auto ccoords = chunk_coords_from_id(id);
    auto cbounds = bounds_from_chunk(id);



    // 1. find chunks from in cube that intersect with curent chunk
    // 2. iterate over the chunks
    // 3. iterate over pixels of input chunks and update aggregator

    double ix_from = (cbounds.s.left - _in_cube->st_reference()->left()) / _in_cube->st_reference()->dx();
    double ix_to = (- 1 + (cbounds.s.right - _in_cube->st_reference()->left()) / _in_cube->st_reference()->dx());
    double iy_from = (_in_cube->st_reference()->top() - cbounds.s.top) / _in_cube->st_reference()->dy();
    double iy_to =  -1 + (_in_cube->st_reference()->top() - cbounds.s.bottom) / _in_cube->st_reference()->dy();

    // some check to be safe
    //if (ix_from < 0) ix_from = 0;
    if (ix_from >= _in_cube->size_x())  return out;

    //if (ix_to >= _in_cube->size_x()) ix_to = _in_cube->size_x() - 1;
    if (ix_to < 0)  return out;

   // if (iy_from < 0) iy_from = 0;
    if (iy_from >= _in_cube->size_y())  return out;

    //if (iy_to >= _in_cube->size_y()) iy_to = _in_cube->size_y() - 1;
    if (iy_to < 0)  return out;

    if (ix_to < ix_from) return out;
    if (iy_to < iy_from) return out;

    // NOTE: Take care because ix_to and iy_to can be > than input cube has pixels and  ix_from and iy_from can be < 0


    std::vector<aggregator_space_slice_singleband*> agg;
    for (uint16_t ib=0; ib<_bands.count(); ++ib) {
        if (_in_func == "min") {
            agg.push_back(new min_aggregtor_space_slice_singleband());
        }
        else if (_in_func == "max") {
            agg.push_back(new max_aggregtor_space_slice_singleband());
        }
        else if (_in_func == "mean") {
            agg.push_back(new mean_aggregtor_space_slice_singleband());
        }
        else if (_in_func == "median") {
            agg.push_back(new median_aggregtor_space_slice_singleband());
        }
        else if (_in_func == "count") {
            agg.push_back(new count_aggregtor_space_slice_singleband());
        }
        else if (_in_func == "var") {
            agg.push_back(new var_aggregtor_space_slice_singleband());
        }
        else if (_in_func == "sd") {
            agg.push_back(new sd_aggregtor_space_slice_singleband());
        }
        else if (_in_func == "prod") {
            agg.push_back(new prod_aggregtor_space_slice_singleband());
        }
        else if (_in_func == "sum") {
            agg.push_back(new sum_aggregtor_space_slice_singleband());
        }
    }

    // find out which chunks need to be read
    chunk_coordinate_tyx in_ccords_from = ccoords;
    in_ccords_from[1] = std::max(int32_t(std::floor(iy_from) / _in_cube->chunk_size()[1]), 0);
    in_ccords_from[2] = std::max(int32_t(std::floor(ix_from) / _in_cube->chunk_size()[2]), 0);

    chunk_coordinate_tyx in_ccords_to = ccoords;
    in_ccords_to[1] = std::min(int32_t(std::floor(iy_to) / _in_cube->chunk_size()[1]), int32_t(_in_cube->count_chunks_y() - 1));
    in_ccords_to[2] = std::min(int32_t(std::floor(ix_to) / _in_cube->chunk_size()[2]), int32_t(_in_cube->count_chunks_x() - 1));

    bool chunk_initialized = false;

    double in_cube_left = _in_cube->st_reference()->left();
    double in_cube_top = _in_cube->st_reference()->top();
    double in_cube_dx = _in_cube->st_reference()->dx();
    double in_cube_dy = _in_cube->st_reference()->dy();
    int32_t in_cube_chunksize_x = _in_cube->chunk_size()[2];
    int32_t in_cube_chunksize_y = _in_cube->chunk_size()[1];

    double out_cube_left = _st_ref->left();
    double out_cube_top = _st_ref->top();
    double out_cube_dx = _st_ref->dx();
    double out_cube_dy = _st_ref->dy();



    // iterate over these chunks
    for (uint32_t ch_y = in_ccords_from[1]; ch_y <= in_ccords_to[1]; ++ch_y) {
        for (uint32_t ch_x = in_ccords_from[2]; ch_x <= in_ccords_to[2]; ++ch_x) {
            chunkid_t input_chunk_id = _in_cube->chunk_id_from_coords({ccoords[0], ch_y, ch_x});
            std::shared_ptr<chunk_data> in_chunk = _in_cube->read_chunk(input_chunk_id);

            // propagate chunk status
            if (in_chunk->status() == chunk_data::chunk_status::ERROR) {
                out->set_status(chunk_data::chunk_status::ERROR);
            }
            else if (in_chunk->status() == chunk_data::chunk_status::INCOMPLETE && out->status() != chunk_data::chunk_status::ERROR) {
                out->set_status(chunk_data::chunk_status::INCOMPLETE);
            }


            if (in_chunk->empty()) {
                continue;
            }

            auto in_chunk_size = in_chunk->size();

            // initialize chunk buffer and aggregators
            if (!chunk_initialized) {
                out->buf(std::calloc(size_btyx[0] * size_btyx[1] * size_btyx[2] * size_btyx[3], sizeof(double)));
                //double *begin = (double *)out->buf();
                //double *end = ((double *)out->buf()) + size_btyx[0] * size_btyx[1] * size_btyx[2] * size_btyx[3];
                //std::fill(begin, end, NAN);
                for (uint16_t ib=0; ib<_bands.count(); ++ib) {
                    agg[ib]->init(&(((double*)(out->buf()))[ib * size_btyx[1] * size_btyx[2] * size_btyx[3]]), size_btyx[1], size_btyx[2], size_btyx[3]);
                }
                chunk_initialized = true;
            }

            // Find out which part of the chunk is needed
            int32_t in_abs_low_x = ch_x  * _in_cube->chunk_size()[2];
            int32_t in_abs_high_x = in_abs_low_x + in_chunk->size()[3] - 1;
            int32_t in_abs_low_y = ch_y  *_in_cube->chunk_size()[1];
            int32_t in_abs_high_y = in_abs_low_y + in_chunk->size()[2]- 1;


            // TODO: consider pixel center????
            int32_t start_x = std::max(in_abs_low_x, int32_t(ix_from));
            int32_t end_x = std::min(in_abs_high_x, int32_t(ix_to));
            int32_t start_y = std::max(in_abs_low_y, int32_t(iy_from));
            int32_t end_y = std::min(in_abs_high_y, int32_t(iy_to));

            // iterate over rectangular subset of input chunk
            // and for each pixel find out to which pixel it will contribute

            double* in_buf = (double*)in_chunk->buf();
            double* out_buf = (double*)out->buf();


            for (int32_t ib = 0; ib<(int32_t)size_btyx[0]; ++ib) {
                for (int32_t it=0; it < (int32_t)size_btyx[1]; ++it) {
                    for (int32_t iy=start_y; iy <= end_y; ++iy) {
                        for (int32_t ix=start_x; ix <= end_x; ++ix) {
                            // find
                            // coordinates of center point from current pixel of the input chunk
                            double in_center_x = (in_cube_left + (ix + 0.5) * in_cube_dx);
                            double in_center_y = (in_cube_top - (iy + 0.5) * in_cube_dy);

                            // corresponding cube coordinates for result cube
                            int32_t out_x_global = int32_t(std::floor((in_center_x - out_cube_left) / out_cube_dx));
                            int32_t out_y_global = int32_t(std::floor((out_cube_top - in_center_y) / out_cube_dy));

                            // local coordinates within chunk
                            int32_t out_x_in_chunk = out_x_global % chunk_size()[2];
                            int32_t out_y_in_chunk = out_y_global % chunk_size()[1];

                            int32_t in_x_in_chunk = ix % in_cube_chunksize_x;
                            int32_t in_y_in_chunk = iy % in_cube_chunksize_y;

                            double* v = &in_buf[ib * in_chunk_size[1] * in_chunk_size[2] * in_chunk_size[3] + it * in_chunk_size[2] * in_chunk_size[3] + in_y_in_chunk *  in_chunk_size[3] + in_x_in_chunk];
                            agg[ib]->combine(&(out_buf[ib * size_btyx[1] * size_btyx[2] * size_btyx[3]]),
                                             v, it, out_y_in_chunk, out_x_in_chunk, size_btyx[1], size_btyx[2], size_btyx[3]);
                        }
                    }
                }
            }


        }
    }

    for (uint16_t ib=0; ib<_bands.count(); ++ib) {
        agg[ib]->finalize(((double*)out->buf()) + ib * size_tyx[0] * size_tyx[1] * size_tyx[2],  size_tyx[0],  size_tyx[1],  size_tyx[2]);
    }
    for (uint16_t ib=0; ib<_bands.count(); ++ib) {
        if (agg[ib] != nullptr) delete agg[ib];
    }
    return out;
}
}  // namespace gdalcubes
