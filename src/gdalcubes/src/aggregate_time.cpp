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
#include "aggregate_time.h"



struct aggregator_time_slice_singleband {
    virtual ~aggregator_time_slice_singleband() {}

    virtual void init(double *out, uint32_t size_x, uint32_t size_y) = 0;
    virtual void combine(double* out, double* in, uint32_t size_x, uint32_t size_y) = 0;
    virtual void finalize(double *out, uint32_t size_x, uint32_t size_y) = 0;
};




/**
 * @brief Implementation of reducer to calculate mean values over time
 */
struct mean_aggregtor_time_slice_singleband : public aggregator_time_slice_singleband {
    void init(double *out, uint32_t size_x, uint32_t size_y) override {

        _count = (uint32_t *)std::calloc(size_x * size_y, sizeof(uint32_t));
        for (uint32_t ixy = 0; ixy < size_x *size_y; ++ixy) {
            _count[ixy] = 0;
            out[ixy] = 0;
        }
    }

    void combine(double* out, double* in, uint32_t size_x, uint32_t size_y) override {
        for (uint32_t ixy = 0; ixy < size_x * size_y; ++ixy) {
            double v = in[ixy];
            if (!std::isnan(v)) {
                out[ixy] += v;
                ++_count[ixy];
            }
        }
    }

    void finalize(double *out, uint32_t size_x, uint32_t size_y) override {
        // divide by count;
        for (uint32_t ixy = 0; ixy < size_x * size_y; ++ixy) {
            out[ixy] =  _count[ixy] > 0 ? out[ixy] / _count[ixy] : NAN;
        }
        std::free(_count);
    }

   private:
    uint32_t *_count;
};


struct min_aggregtor_time_slice_singleband : public aggregator_time_slice_singleband {
    void init(double *out, uint32_t size_x, uint32_t size_y) override {
        for (uint32_t ixy = 0; ixy < size_x *size_y; ++ixy) {
            out[ixy] = NAN;
        }
    }
    void combine(double* out, double* in, uint32_t size_x, uint32_t size_y) override {
        for (uint32_t ixy = 0; ixy < size_x * size_y; ++ixy) {
            double v = in[ixy];
            if (!std::isnan(v)) {
                if (std::isnan(out[ixy])) {
                    out[ixy] = v;
                }
                else {
                    out[ixy] = std::min(out[ixy], v);
                }
            }
        }
    }

    void finalize(double *out, uint32_t size_x, uint32_t size_y) override {}
};

struct max_aggregtor_time_slice_singleband : public aggregator_time_slice_singleband {
    void init(double *out, uint32_t size_x, uint32_t size_y) override {
        for (uint32_t ixy = 0; ixy < size_x *size_y; ++ixy) {
            out[ixy] = NAN;
        }
    }
    void combine(double* out, double* in, uint32_t size_x, uint32_t size_y) override {
        for (uint32_t ixy = 0; ixy < size_x * size_y; ++ixy) {
            double v = in[ixy];
            if (!std::isnan(v)) {
                if (std::isnan(out[ixy])) {
                    out[ixy] = v;
                }
                else {
                    out[ixy] = std::max(out[ixy], v);
                }
            }
        }
    }

    void finalize(double *out, uint32_t size_x, uint32_t size_y) override {}
};


struct sum_aggregtor_time_slice_singleband : public aggregator_time_slice_singleband {
    void init(double *out, uint32_t size_x, uint32_t size_y) override {
        for (uint32_t ixy = 0; ixy < size_x *size_y; ++ixy) {
            out[ixy] = NAN;
        }
    }
    void combine(double* out, double* in, uint32_t size_x, uint32_t size_y) override {
        for (uint32_t ixy = 0; ixy < size_x * size_y; ++ixy) {
            double v = in[ixy];
            if (!std::isnan(v)) {
                if (std::isnan(out[ixy])) {
                    out[ixy] = v;
                }
                else {
                    out[ixy] += v;
                }
            }
        }
    }

    void finalize(double *out, uint32_t size_x, uint32_t size_y) override {}
};



struct prod_aggregtor_time_slice_singleband : public aggregator_time_slice_singleband {
    void init(double *out, uint32_t size_x, uint32_t size_y) override {
        for (uint32_t ixy = 0; ixy < size_x *size_y; ++ixy) {
            out[ixy] = NAN;
        }
    }
    void combine(double* out, double* in, uint32_t size_x, uint32_t size_y) override {
        for (uint32_t ixy = 0; ixy < size_x * size_y; ++ixy) {
            double v = in[ixy];
            if (!std::isnan(v)) {
                if (std::isnan(out[ixy])) {
                    out[ixy] = v;
                }
                else {
                    out[ixy] *= v;
                }
            }
        }
    }

    void finalize(double *out, uint32_t size_x, uint32_t size_y) override {}
};



struct count_aggregtor_time_slice_singleband : public aggregator_time_slice_singleband {
    void init(double *out, uint32_t size_x, uint32_t size_y) override {
        for (uint32_t ixy = 0; ixy < size_x *size_y; ++ixy) {
            out[ixy] = 0;
        }
    }
    void combine(double* out, double* in, uint32_t size_x, uint32_t size_y) override {
        for (uint32_t ixy = 0; ixy < size_x * size_y; ++ixy) {
            double v = in[ixy];
            if (!std::isnan(v)) {
                out[ixy] += 1.0;
            }
        }
    }

    void finalize(double *out, uint32_t size_x, uint32_t size_y) override {}
};




struct median_aggregtor_time_slice_singleband : public aggregator_time_slice_singleband {
    void init(double *out, uint32_t size_x, uint32_t size_y) override {
        _m_buckets.resize(size_x *size_y, std::vector<double>());
        for (uint32_t ixy = 0; ixy < size_x *size_y; ++ixy) {
            out[ixy] = NAN;
        }
    }

    void combine(double* out, double* in, uint32_t size_x, uint32_t size_y) override {
        for (uint32_t ixy = 0; ixy < size_x * size_y; ++ixy) {
            double v = in[ixy];
            if (!std::isnan(v)) {
                _m_buckets[ixy].push_back(v);
            }
        }
    }

    void finalize(double *out, uint32_t size_x, uint32_t size_y) override {
        for (uint32_t ixy = 0; ixy < size_x * size_y; ++ixy) {
            std::vector<double> &list = _m_buckets[ixy];
            std::sort(list.begin(), list.end());
            if (list.size() == 0) {
                out[ixy] = NAN;
            } else if (list.size() % 2 == 1) {
                out[ixy] = list[list.size() / 2];
            } else {
                out[ixy] = (list[list.size() / 2] + list[list.size() / 2 - 1]) / ((double)2);
            }
        }
    }
   private:
    std::vector<std::vector<double>> _m_buckets;
};


/**
 * @brief Implementation of reducer to calculate variance values over time using Welford's Online algorithm
 */
struct var_aggregtor_time_slice_singleband : public aggregator_time_slice_singleband {
    void init(double *out, uint32_t size_x, uint32_t size_y) override {
        _count = (uint32_t *)std::calloc(size_x * size_y, sizeof(uint32_t));
        _mean = (double *)std::calloc(size_x * size_y, sizeof(double));

        for (uint32_t ixy = 0; ixy < size_x *size_y; ++ixy) {
            _count[ixy] = 0;
            _mean[ixy] = 0;
            out[ixy] = 0;
        }
    }
    void combine(double* out, double* in, uint32_t size_x, uint32_t size_y) override {
        for (uint32_t ixy = 0; ixy < size_x * size_y; ++ixy) {
            double v = in[ixy];
            if (!std::isnan(v)) {
                double &mean = _mean[ixy];
                uint32_t &count = _count[ixy];
                ++count;
                double delta = v - mean;
                mean += delta / count;
                out[ixy] += delta * (v - mean);
            }
        }
    }

    void finalize(double *out, uint32_t size_x, uint32_t size_y) override {
        for (uint32_t ixy = 0; ixy < size_x *size_y; ++ixy) {
            out[ixy] = _count[ixy] > 1 ? out[ixy] / (_count[ixy] - 1) : NAN;
        }
        std::free(_count);
        std::free(_mean);
    }

   protected:
    uint32_t *_count;
    double *_mean;
};


struct sd_aggregtor_time_slice_singleband : public var_aggregtor_time_slice_singleband {
    void finalize(double *out, uint32_t size_x, uint32_t size_y) override {
        for (uint32_t ixy = 0; ixy < size_x *size_y; ++ixy) {
            out[ixy] = _count[ixy] > 1 ? std::sqrt(out[ixy] / (_count[ixy] - 1)) : NAN;
        }
        std::free(_count);
        std::free(_mean);
    }
};


namespace gdalcubes {

std::shared_ptr<chunk_data> aggregate_time_cube::read_chunk(chunkid_t id) {
    GCBS_TRACE("aggregate_time_cube::read_chunk(" + std::to_string(id) + ")");
    std::shared_ptr<chunk_data> out = std::make_shared<chunk_data>();
    if (id >= count_chunks())
        return out;  // chunk is outside of the view, we don't need to read anything.


    coords_nd<uint32_t, 3> size_tyx = chunk_size(id);
    coords_nd<uint32_t, 4> size_btyx = {uint32_t(_bands.count()), size_tyx[0], size_tyx[1], size_tyx[2]};
    out->size(size_btyx);

    // Fill buffers accordingly
    out->buf(std::calloc(size_btyx[0] * size_btyx[1] * size_btyx[2] * size_btyx[3], sizeof(double)));
    double *begin = (double *)out->buf();
    double *end = ((double *)out->buf()) + size_btyx[0] * size_btyx[1] * size_btyx[2] * size_btyx[3];
    std::fill(begin, end, NAN);

    auto climits = chunk_limits(id);
    auto ccoords = chunk_coords_from_id(id);

    std::map<chunkid_t, std::shared_ptr<chunk_data>> chunk_cache;

    std::vector<aggregator_time_slice_singleband*> agg;
    for (uint16_t ib=0; ib<_bands.count(); ++ib) {
        if (_in_func == "min") {
            agg.push_back(new min_aggregtor_time_slice_singleband());
        }
        else if (_in_func == "max") {
            agg.push_back(new max_aggregtor_time_slice_singleband());
        }
        else if (_in_func == "mean") {
            agg.push_back(new mean_aggregtor_time_slice_singleband());
        }
        else if (_in_func == "median") {
            agg.push_back(new median_aggregtor_time_slice_singleband());
        }
        else if (_in_func == "count") {
            agg.push_back(new count_aggregtor_time_slice_singleband());
        }
        else if (_in_func == "var") {
            agg.push_back(new var_aggregtor_time_slice_singleband());
        }
        else if (_in_func == "sd") {
            agg.push_back(new sd_aggregtor_time_slice_singleband());
        }
        else if (_in_func == "prod") {
            agg.push_back(new prod_aggregtor_time_slice_singleband());
        }
        else if (_in_func == "sum") {
            agg.push_back(new sum_aggregtor_time_slice_singleband());
        }
    }



    for (uint32_t it=0; it<size_tyx[0]; ++it) {

        datetime t_cur = _st_ref->datetime_at_index(climits.low[0] + it);
        datetime t_next =  _st_ref->datetime_at_index(climits.low[0] + it + 1);

        uint32_t first = 0;
        uint32_t last = 0;

        t_cur.unit(_in_cube->st_reference()->dt_unit());
        t_next.unit( _in_cube->st_reference()->dt_unit());

        if (cube_stref::type_string(_in_cube->st_reference()) == "cube_stref_regular") {
            first = _in_cube->st_reference()->index_at_datetime(t_cur);
            if (_in_cube->st_reference()->datetime_at_index(first) < t_cur) {
                ++first;
            }
            last = _in_cube->st_reference()->index_at_datetime(t_next);
            if (_in_cube->st_reference()->datetime_at_index(last) >= t_next) {
                --last;
            }
        }
        else if (cube_stref::type_string(_in_cube->st_reference()) == "cube_stref_labeled_time") {
            auto p = std::dynamic_pointer_cast<cube_stref_labeled_time>(_in_cube->st_reference());
            while (p->datetime_at_index(first) < t_cur) {
                ++first;
            }
            if (p->datetime_at_index(first) >= t_next) {
                continue; // labeled time axis of input cube as a gap larger than new time duration of cells
            }
            last = first;
            while (p->datetime_at_index(last) < t_next) {
                ++last;
            }
            --last;
        }
        if (last < first) {
            // TODO: exception or empty time slice if labeled time axis?!
            GCBS_DEBUG("Aggregation state points to invalid input time points, ignoring time slice");
            continue;
        }

        for (uint16_t ib=0; ib<_bands.count(); ++ib) {
            agg[ib]->init(((double*)out->buf()) + ib * size_tyx[0] * size_tyx[1] * size_tyx[2] + it * size_tyx[1] * size_tyx[2],
                              size_tyx[1], size_tyx[2]);
        }
        for (uint32_t i=first; i<=last; ++i) {

            if (i >= _in_cube->st_reference()->nt()) {
                break;
            }

            // Now, find out which chunk from input cube needs to be read
            chunk_coordinate_tyx in_ccords = ccoords;
            in_ccords[0] = i / _in_cube->chunk_size()[0];
            chunkid_t cur_in_chunk = _in_cube->chunk_id_from_coords(in_ccords);

            if (chunk_cache.find(cur_in_chunk) == chunk_cache.end()) {
                chunk_cache[cur_in_chunk] = _in_cube->read_chunk(cur_in_chunk);
                // TODO: remove old chunk from "cache" and use a single pointer instead of map?!
            }

            std::shared_ptr<chunk_data> in_chunk = chunk_cache[cur_in_chunk];

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

            // feed aggregator function for each band
            for (uint16_t ib=0; ib<_bands.count(); ++ib) {
                agg[ib]->combine(((double*)out->buf()) + ib * size_tyx[0] * size_tyx[1] * size_tyx[2] + it * size_tyx[1] * size_tyx[2],
                             ((double*)in_chunk->buf()) +  ib * in_chunk->size()[1] * in_chunk->size()[2] * in_chunk->size()[3] + (i % _in_cube->chunk_size()[0]) * in_chunk->size()[2] * in_chunk->size()[3],
                             size_tyx[1], size_tyx[2]);
            }
        }
        for (uint16_t ib=0; ib<_bands.count(); ++ib) {
            agg[ib]->finalize(((double*)out->buf()) + ib * size_tyx[0] * size_tyx[1] * size_tyx[2] + it * size_tyx[1] * size_tyx[2],
                              size_tyx[1], size_tyx[2]);
        }
    }

    for (uint16_t ib=0; ib<_bands.count(); ++ib) {
        if (agg[ib] != nullptr) delete agg[ib];
    }
    return out;
}
}  // namespace gdalcubes
