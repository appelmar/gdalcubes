/*
    MIT License

    Copyright (c) 2019 Marius Appel <marius.appel@hs-bochum.de>

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

#include "reduce_space.h"

namespace gdalcubes {

struct reducer_singleband_s {
    virtual ~reducer_singleband_s() {}

    /**
     * @brief Initialization function for reducers that is automatically called before reading data from the input cube
     * @param a chunk data where reduction results are written to later
     * @param band_idx_in over which band of the chunk data (zero-based index) shall the reducer be applied?
     * @param band_idx_out to which band of the result chunk (zero-based index) shall the reducer write?
     */
    virtual void init(std::shared_ptr<chunk_data> a, uint16_t band_idx_in, uint16_t band_idx_out, std::shared_ptr<cube> in_cube) = 0;

    /**
     * @brief Combines a chunk of data from the input cube with the current state of the result chunk according to the specific reducer
     * @param a output chunk of the reduction
     * @param b one input chunk of the input cube, which is aligned with the output chunk in space
     */
    virtual void combine(std::shared_ptr<chunk_data> a, std::shared_ptr<chunk_data> b, chunkid_t chunk_id) = 0;

    /**
     * @brief Finallizes the reduction, i.e., frees additional buffers and postprocesses the result (e.g. dividing by n for mean reducer)
     * @param a result chunk
     */
    virtual void finalize(std::shared_ptr<chunk_data> a) = 0;
};

/**
 * @brief Implementation of reducer to calculate sum values over time
 */
struct sum_reducer_singleband_s : public reducer_singleband_s {
    void init(std::shared_ptr<chunk_data> a, uint16_t band_idx_in, uint16_t band_idx_out, std::shared_ptr<cube> in_cube) override {
        _band_idx_in = band_idx_in;
        _band_idx_out = band_idx_out;
        for (uint32_t it = 0; it < a->size()[1]; ++it) {
            ((double *)a->buf())[_band_idx_out * a->size()[1] + it] = 0;
        }
    }

    void combine(std::shared_ptr<chunk_data> a, std::shared_ptr<chunk_data> b, chunkid_t chunk_id) override {
        for (uint32_t it = 0; it < b->size()[1]; ++it) {
            for (uint32_t ixy = 0; ixy < b->size()[2] * b->size()[3]; ++ixy) {
                double v = ((double *)b->buf())[_band_idx_in * b->size()[1] * b->size()[2] * b->size()[3] + it * b->size()[2] * b->size()[3] + ixy];
                if (!std::isnan(v)) {
                    ((double *)a->buf())[_band_idx_out * a->size()[1] * a->size()[2] * a->size()[3] + it] += v;
                }
            }
        }
    }
    void finalize(std::shared_ptr<chunk_data> a) override {}

   private:
    uint16_t _band_idx_in;
    uint16_t _band_idx_out;
};

/**
 * @brief Implementation of reducer to calculate product values over time
 */
struct prod_reducer_singleband_s : public reducer_singleband_s {
    void init(std::shared_ptr<chunk_data> a, uint16_t band_idx_in, uint16_t band_idx_out, std::shared_ptr<cube> in_cube) override {
        _band_idx_in = band_idx_in;
        _band_idx_out = band_idx_out;
        for (uint32_t it = 0; it < a->size()[1]; ++it) {
            ((double *)a->buf())[_band_idx_out * a->size()[1] + it] = 1;
        }
    }

    void combine(std::shared_ptr<chunk_data> a, std::shared_ptr<chunk_data> b, chunkid_t chunk_id) override {
        for (uint32_t it = 0; it < b->size()[1]; ++it) {
            for (uint32_t ixy = 0; ixy < b->size()[2] * b->size()[3]; ++ixy) {
                double v = ((double *)b->buf())[_band_idx_in * b->size()[1] * b->size()[2] * b->size()[3] + it * b->size()[2] * b->size()[3] + ixy];
                if (!std::isnan(v)) {
                    ((double *)a->buf())[_band_idx_out * a->size()[1] * a->size()[2] * a->size()[3] + it] *= v;
                }
            }
        }
    }
    void finalize(std::shared_ptr<chunk_data> a) override {}

   private:
    uint16_t _band_idx_in;
    uint16_t _band_idx_out;
};

/**
 * @brief Implementation of reducer to calculate mean values over time
 */
struct mean_reducer_singleband_s : public reducer_singleband_s {
    void init(std::shared_ptr<chunk_data> a, uint16_t band_idx_in, uint16_t band_idx_out, std::shared_ptr<cube> in_cube) override {
        _band_idx_in = band_idx_in;
        _band_idx_out = band_idx_out;
        _count = (uint32_t *)std::calloc(a->size()[1], sizeof(uint32_t));
        for (uint32_t it = 0; it < a->size()[1]; ++it) {
            _count[it] = 0;
            ((double *)a->buf())[_band_idx_out * a->size()[1] + it] = 0;
        }
    }

    void combine(std::shared_ptr<chunk_data> a, std::shared_ptr<chunk_data> b, chunkid_t chunk_id) override {
        for (uint32_t it = 0; it < b->size()[1]; ++it) {
            for (uint32_t ixy = 0; ixy < b->size()[2] * b->size()[3]; ++ixy) {
                double v = ((double *)b->buf())[_band_idx_in * b->size()[1] * b->size()[2] * b->size()[3] + it * b->size()[2] * b->size()[3] + ixy];
                if (!std::isnan(v)) {
                    ((double *)a->buf())[_band_idx_out * a->size()[1] * a->size()[2] * a->size()[3] + it] += v;
                    ++_count[it];
                }
            }
        }
    }

    void finalize(std::shared_ptr<chunk_data> a) override {
        // divide by count;
        for (uint32_t it = 0; it < a->size()[1]; ++it) {
            ((double *)a->buf())[_band_idx_out * a->size()[1] + it] = _count[it] > 0 ? ((double *)a->buf())[_band_idx_out * a->size()[1] + it] / _count[it] : NAN;
        }
        std::free(_count);
    }

   private:
    uint32_t *_count;
    uint16_t _band_idx_in;
    uint16_t _band_idx_out;
};

/**
 * @brief Implementation of reducer to calculate minimum values over time
 */
struct min_reducer_singleband_s : public reducer_singleband_s {
    void init(std::shared_ptr<chunk_data> a, uint16_t band_idx_in, uint16_t band_idx_out, std::shared_ptr<cube> in_cube) override {
        _band_idx_in = band_idx_in;
        _band_idx_out = band_idx_out;
        for (uint32_t it = 0; it < a->size()[1]; ++it) {
            ((double *)a->buf())[_band_idx_out * a->size()[1] + it] = NAN;
        }
    }

    void combine(std::shared_ptr<chunk_data> a, std::shared_ptr<chunk_data> b, chunkid_t chunk_id) override {
        for (uint32_t it = 0; it < b->size()[1]; ++it) {
            double *w = &(((double *)a->buf())[_band_idx_out * a->size()[1] * a->size()[2] * a->size()[3] + it]);
            for (uint32_t ixy = 0; ixy < b->size()[2] * b->size()[3]; ++ixy) {
                double v = ((double *)b->buf())[_band_idx_in * b->size()[1] * b->size()[2] * b->size()[3] + it * b->size()[2] * b->size()[3] + ixy];
                if (!std::isnan(v)) {
                    if (std::isnan(*w))
                        *w = v;
                    else
                        *w = std::min(*w, v);
                }
            }
        }
    }

    void finalize(std::shared_ptr<chunk_data> a) override {}

   private:
    uint16_t _band_idx_in;
    uint16_t _band_idx_out;
};

/**
 * @brief Implementation of reducer to calculate maximum values over time
 */
struct max_reducer_singleband_s : public reducer_singleband_s {
    void init(std::shared_ptr<chunk_data> a, uint16_t band_idx_in, uint16_t band_idx_out, std::shared_ptr<cube> in_cube) override {
        _band_idx_in = band_idx_in;
        _band_idx_out = band_idx_out;
        for (uint32_t it = 0; it < a->size()[1]; ++it) {
            ((double *)a->buf())[_band_idx_out * a->size()[1] + it] = NAN;
        }
    }

    void combine(std::shared_ptr<chunk_data> a, std::shared_ptr<chunk_data> b, chunkid_t chunk_id) override {
        for (uint32_t it = 0; it < b->size()[1]; ++it) {
            double *w = &(((double *)a->buf())[_band_idx_out * a->size()[1] * a->size()[2] * a->size()[3] + it]);
            for (uint32_t ixy = 0; ixy < b->size()[2] * b->size()[3]; ++ixy) {
                double v = ((double *)b->buf())[_band_idx_in * b->size()[1] * b->size()[2] * b->size()[3] + it * b->size()[2] * b->size()[3] + ixy];
                if (!std::isnan(v)) {
                    if (std::isnan(*w))
                        *w = v;
                    else
                        *w = std::max(*w, v);
                }
            }
        }
    }

    void finalize(std::shared_ptr<chunk_data> a) override {}

   private:
    uint16_t _band_idx_in;
    uint16_t _band_idx_out;
};

struct count_reducer_singleband_s : public reducer_singleband_s {
    void init(std::shared_ptr<chunk_data> a, uint16_t band_idx_in, uint16_t band_idx_out, std::shared_ptr<cube> in_cube) override {
        _band_idx_in = band_idx_in;
        _band_idx_out = band_idx_out;
        for (uint32_t it = 0; it < a->size()[1]; ++it) {
            ((double *)a->buf())[_band_idx_out * a->size()[1] + it] = 0;
        }
    }

    void combine(std::shared_ptr<chunk_data> a, std::shared_ptr<chunk_data> b, chunkid_t chunk_id) override {
        for (uint32_t it = 0; it < b->size()[1]; ++it) {
            for (uint32_t ixy = 0; ixy < b->size()[2] * b->size()[3]; ++ixy) {
                double v = ((double *)b->buf())[_band_idx_in * b->size()[1] * b->size()[2] * b->size()[3] + it * b->size()[2] * b->size()[3] + ixy];
                if (!std::isnan(v)) {
                    double *w = &(((double *)a->buf())[_band_idx_out * a->size()[1] * a->size()[2] * a->size()[3] + it]);
                    *w += 1;
                }
            }
        }
    }

    void finalize(std::shared_ptr<chunk_data> a) override {}

   private:
    uint16_t _band_idx_in;
    uint16_t _band_idx_out;
};

/**
 * @brief Implementation of reducer to calculate median values over time
 * @note Calculating the exact median has a strong memory overhead, approximate median reducers might be implemented in the future
 */
struct median_reducer_singleband_s : public reducer_singleband_s {
    void init(std::shared_ptr<chunk_data> a, uint16_t band_idx_in, uint16_t band_idx_out, std::shared_ptr<cube> in_cube) override {
        _band_idx_in = band_idx_in;
        _band_idx_out = band_idx_out;
        _m_buckets.resize(a->size()[1], std::vector<double>());
    }

    void combine(std::shared_ptr<chunk_data> a, std::shared_ptr<chunk_data> b, chunkid_t chunk_id) override {
        for (uint32_t it = 0; it < b->size()[1]; ++it) {
            for (uint32_t ixy = 0; ixy < b->size()[2] * b->size()[3]; ++ixy) {
                double v = ((double *)b->buf())[_band_idx_in * b->size()[1] * b->size()[2] * b->size()[3] + it * b->size()[2] * b->size()[3] + ixy];
                if (!std::isnan(v)) {
                    _m_buckets[it].push_back(v);
                }
            }
        }
    }

    void finalize(std::shared_ptr<chunk_data> a) override {
        for (uint32_t it = 0; it < a->size()[1]; ++it) {
            std::vector<double> &list = _m_buckets[it];
            std::sort(list.begin(), list.end());
            if (list.size() == 0) {
                ((double *)a->buf())[_band_idx_out * a->size()[1] + it] = NAN;
            } else if (list.size() % 2 == 1) {
                ((double *)a->buf())[_band_idx_out * a->size()[1] + it] = list[list.size() / 2];
            } else {
                ((double *)a->buf())[_band_idx_out * a->size()[1] + it] = (list[list.size() / 2] + list[list.size() / 2 - 1]) / ((double)2);
            }
        }
    }

   private:
    std::vector<std::vector<double>> _m_buckets;
    uint16_t _band_idx_in;
    uint16_t _band_idx_out;
};

/**
 * @brief Implementation of reducer to calculate variance values over time using Welford's Online algorithm
 */
struct var_reducer_singleband_s : public reducer_singleband_s {
    void init(std::shared_ptr<chunk_data> a, uint16_t band_idx_in, uint16_t band_idx_out, std::shared_ptr<cube> in_cube) override {
        _band_idx_in = band_idx_in;
        _band_idx_out = band_idx_out;
        _count = (uint32_t *)std::calloc(a->size()[1], sizeof(uint32_t));
        _mean = (double *)std::calloc(a->size()[1], sizeof(double));
        for (uint32_t it = 0; it < a->size()[1]; ++it) {
            _count[it] = 0;
            _mean[it] = 0;
            ((double *)a->buf())[_band_idx_out * a->size()[1] + it] = 0;
        }
    }

    void combine(std::shared_ptr<chunk_data> a, std::shared_ptr<chunk_data> b, chunkid_t chunk_id) override {
        for (uint32_t it = 0; it < b->size()[1]; ++it) {
            for (uint32_t ixy = 0; ixy < b->size()[2] * b->size()[3]; ++ixy) {
                double &v = ((double *)b->buf())[_band_idx_in * b->size()[1] * b->size()[2] * b->size()[3] + it * b->size()[2] * b->size()[3] + ixy];
                if (!std::isnan(v)) {
                    double &mean = _mean[it];
                    uint32_t &count = _count[it];
                    ++count;
                    double delta = v - mean;
                    mean += delta / count;
                    ((double *)a->buf())[_band_idx_out * a->size()[1] * a->size()[2] * a->size()[3] + it] += delta * (v - mean);
                }
            }
        }
    }

    virtual void finalize(std::shared_ptr<chunk_data> a) override {
        // divide by count - 1;
        for (uint32_t it = 0; it < a->size()[1]; ++it) {
            ((double *)a->buf())[_band_idx_out * a->size()[1] + it] = _count[it] > 1 ? ((double *)a->buf())[_band_idx_out * a->size()[1] + it] / (_count[it] - 1) : NAN;
        }
        std::free(_count);
        std::free(_mean);
    }

   protected:
    uint32_t *_count;
    double *_mean;
    uint16_t _band_idx_in;
    uint16_t _band_idx_out;
};

/**
 * @brief Implementation of reducer to calculate standard deviation values over time using Welford's Online algorithm
 */
struct sd_reducer_singleband_s : public var_reducer_singleband_s {
    void finalize(std::shared_ptr<chunk_data> a) override {
        // divide by count - 1;
        for (uint32_t it = 0; it < a->size()[1]; ++it) {
            ((double *)a->buf())[_band_idx_out * a->size()[1] + it] = _count[it] > 1 ? sqrt(((double *)a->buf())[_band_idx_out * a->size()[1] + it] / (_count[it] - 1)) : NAN;
        }
        std::free(_count);
        std::free(_mean);
    }
};

std::shared_ptr<chunk_data> reduce_space_cube::read_chunk(chunkid_t id) {
    GCBS_TRACE("reduce_space_cube::read_chunk(" + std::to_string(id) + ")");
    std::shared_ptr<chunk_data> out = std::make_shared<chunk_data>();
    if (id >= count_chunks())
        return out;  // chunk is outside of the view, we don't need to read anything.

    // If input cube is already "reduced", simply return corresponding input chunk
    if (_in_cube->size_y() == 1 && _in_cube->size_x() == 1) {
        return _in_cube->read_chunk(id);
    }

    coords_nd<uint32_t, 3> size_tyx = chunk_size(id);
    coords_nd<uint32_t, 4> size_btyx = {uint32_t(_reducer_bands.size()), size_tyx[0], 1, 1};
    out->size(size_btyx);

    std::vector<reducer_singleband_s *> reducers;
    for (uint16_t i = 0; i < _reducer_bands.size(); ++i) {
        reducer_singleband_s *r = nullptr;
        if (_reducer_bands[i].first == "min") {
            r = new min_reducer_singleband_s();
        } else if (_reducer_bands[i].first == "max") {
            r = new max_reducer_singleband_s();
        } else if (_reducer_bands[i].first == "mean") {
            r = new mean_reducer_singleband_s();
        } else if (_reducer_bands[i].first == "median") {
            r = new median_reducer_singleband_s();
        } else if (_reducer_bands[i].first == "sum") {
            r = new sum_reducer_singleband_s();
        } else if (_reducer_bands[i].first == "count") {
            r = new count_reducer_singleband_s();
        } else if (_reducer_bands[i].first == "prod") {
            r = new prod_reducer_singleband_s();
        } else if (_reducer_bands[i].first == "var") {
            r = new var_reducer_singleband_s();
        } else if (_reducer_bands[i].first == "sd") {
            r = new sd_reducer_singleband_s();
        } else
            throw std::string("ERROR in reduce_time_cube::read_chunk(): Unknown reducer given");

        reducers.push_back(r);
    }

    // iterate over all chunks that must be read from the input cube to compute this chunk
    bool empty = true;
    bool initialized = false; // lazy initialization after the first non-empty chunk
    for (chunkid_t i = id * _in_cube->count_chunks_x() * _in_cube->count_chunks_y(); i < (id + 1) * _in_cube->count_chunks_x() * _in_cube->count_chunks_y(); ++i) {
        std::shared_ptr<chunk_data> x = _in_cube->read_chunk(i);

        // propagate chunk status
        if (x->status() == chunk_data::chunk_status::ERROR) {
            out->set_status(chunk_data::chunk_status::ERROR);
        }
        else if (x->status() == chunk_data::chunk_status::INCOMPLETE && out->status() != chunk_data::chunk_status::ERROR) {
            out->set_status(chunk_data::chunk_status::INCOMPLETE);
        }

        if (!x->empty()) {
            if (!initialized) {
                // Fill buffers with NAN
                out->buf(std::calloc(size_btyx[0] * size_btyx[1] * size_btyx[2] * size_btyx[3], sizeof(double)));
                double *begin = (double *)out->buf();
                double *end = ((double *)out->buf()) + size_btyx[0] * size_btyx[1] * size_btyx[2] * size_btyx[3];
                std::fill(begin, end, NAN);
                for (uint16_t ib = 0; ib < _reducer_bands.size(); ++ib) {
                    uint16_t band_idx_in = _in_cube->bands().get_index(_reducer_bands[ib].second);
                    reducers[ib]->init(out, band_idx_in, ib, _in_cube);
                }
                initialized = true;
            }
            for (uint16_t ib = 0; ib < _reducer_bands.size(); ++ib) {
                reducers[ib]->combine(out, x, i);
            }
            empty = false;
        }
    }
    if (empty) {
        auto s = out->status();
        out = std::make_shared<chunk_data>();
        out->set_status(s);
    }
    else {
        for (uint16_t i = 0; i < _reducer_bands.size(); ++i) {
            reducers[i]->finalize(out);
        }
    }
    for (uint16_t i = 0; i < reducers.size(); ++i) {
        if (reducers[i]) delete reducers[i];
    }
    return out;
}

}  // namespace gdalcubes
