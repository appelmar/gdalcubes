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
#include "window_time.h"

namespace gdalcubes {

std::function<double(double* buf, uint16_t n)> window_time_cube::get_default_reducer_by_name(std::string name) {
    if (name == "mean") {
        return std::function<double(double* buf, uint16_t n)>([](double* buf, uint16_t n) {
            double sum = 0.0;
            uint32_t count = 0;
            for (uint16_t i = 0; i < n; ++i) {
                if (!std::isnan(buf[i])) {
                    sum += buf[i];
                    ++count;
                }
            }
            return sum / (double)count;
        });
    } else if (name == "sum") {
        return std::function<double(double* buf, uint16_t n)>([](double* buf, uint16_t n) {
            double sum = 0.0;
            for (uint16_t i = 0; i < n; ++i) {
                if (!std::isnan(buf[i]))
                    sum += buf[i];
            }
            return sum;
        });
    } else if (name == "count") {
        return std::function<double(double* buf, uint16_t n)>([](double* buf, uint16_t n) {
            double count = 0.0;
            for (uint16_t i = 0; i < n; ++i) {
                if (!std::isnan(buf[i])) ++count;
            }
            return count;
        });
    } else if (name == "prod") {
        return std::function<double(double* buf, uint16_t n)>([](double* buf, uint16_t n) {
            double prod = 1.0;
            for (uint16_t i = 0; i < n; ++i) {
                if (!std::isnan(buf[i])) {
                    prod *= buf[i];
                }
            }
            return prod;
        });
    } else if (name == "min") {
        return std::function<double(double* buf, uint16_t n)>([](double* buf, uint16_t n) {
            double min = std::numeric_limits<double>::max();
            for (uint16_t i = 0; i < n; ++i) {
                if (!std::isnan(buf[i])) {
                    min = std::min(min, buf[i]);
                }
            }
            return min;
        });
    } else if (name == "max") {
        return std::function<double(double* buf, uint16_t n)>([](double* buf, uint16_t n) {
            double max = std::numeric_limits<double>::min();
            for (uint16_t i = 0; i < n; ++i) {
                if (!std::isnan(buf[i])) {
                    max = std::max(max, buf[i]);
                }
            }
            return max;
        });
    } else if (name == "median") {
        return std::function<double(double* buf, uint16_t n)>([](double* buf, uint16_t n) {
            std::vector<double> val;
            val.assign(buf, buf+n); // TODO: handle NAN?
            std::sort(val.begin(), val.end());
             if (val.size() % 2 == 1) {
                return val[val.size() / 2];
        } else {
            return (val[val.size() / 2] + val[val.size() / 2 - 1]) / ((double)2);
        } });
    } else {
        throw std::string("ERROR in window_time_cube::get_default_reducer_by_name(): Unknown reducer '" + name + "'");
    }
}

std::function<double(double* buf, uint16_t n)> window_time_cube::get_kernel_reducer(std::vector<double> kernel) {
    if (kernel.size() != (uint32_t)(_win_size_l + 1 + _win_size_r)) {
        throw std::string("ERROR in window_time_cube::get_kernel_reducer(): Size of kernel does not match size of window");
    }

    return std::function<double(double* buf, uint16_t n)>([kernel](double* buf, uint16_t n) -> double {
        double v = 0.0;
        for (uint16_t i = 0; i < n; ++i) {
            if (std::isnan(buf[i])) {
                return NAN;
            }
            v += buf[i] * kernel[i];
        }
        return v;
    });
}

std::shared_ptr<chunk_data> window_time_cube::read_chunk(chunkid_t id) {
    GCBS_TRACE("window_time_cube::read_chunk(" + std::to_string(id) + ")");
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

    uint32_t chunk_count_l = (uint32_t)std::ceil((double)_win_size_l / (double)(_in_cube->chunk_size()[0]));
    uint32_t chunk_count_r = (uint32_t)std::ceil((double)_win_size_r / (double)(_in_cube->chunk_size()[0]));

    std::shared_ptr<chunk_data> this_chunk = _in_cube->read_chunk(id);
    std::vector<std::shared_ptr<chunk_data>> l_chunks;
    std::vector<std::shared_ptr<chunk_data>> r_chunks;

    // Read needed chunks depending on window and chunk sizes
    for (uint16_t i = 1; i <= chunk_count_l; ++i) {
        // read l chunks
        int32_t tid = id - i * (_in_cube->count_chunks_x() * _in_cube->count_chunks_y());
        if (tid < 0) break;
        l_chunks.push_back(_in_cube->read_chunk(tid));
    }
    for (uint16_t i = 1; i <= chunk_count_r; ++i) {
        // read l chunks
        int32_t tid = id + i * (_in_cube->count_chunks_x() * _in_cube->count_chunks_y());
        if (tid >= (int32_t)_in_cube->count_chunks()) break;
        r_chunks.push_back(_in_cube->read_chunk(tid));
    }

    // buffer for a single time series including data from adjacent chunks for all used input bands
    uint32_t cur_ts_length = _win_size_l + size_tyx[0] + _win_size_r;
    double* cur_ts = (double*)std::calloc(cur_ts_length * _bands.count(), sizeof(double));
    std::fill(cur_ts, cur_ts + cur_ts_length * _bands.count(), NAN);

    for (uint32_t ixy = 0; ixy < size_tyx[1] * size_tyx[2]; ++ixy) {
        // fill values from l chunks
        int32_t tsidx = _win_size_l - 1;
        for (uint16_t i = 0; i < chunk_count_l; ++i) {
            if (i >= l_chunks.size()) {
                // fill NA
                for (uint32_t ic = 0; ic < _in_cube->chunk_size()[0]; ++ic) {
                    if (tsidx >= 0 && tsidx < (int32_t)cur_ts_length) {
                        for (uint16_t ib = 0; ib < _bands.count(); ++ib) {
                            cur_ts[ib * cur_ts_length + tsidx] = NAN;
                        }
                        tsidx--;
                    }
                }
            } else {
                for (uint32_t ic = 0; ic < l_chunks[i]->size()[1]; ++ic) {
                    if (tsidx >= 0 && tsidx < (int32_t)cur_ts_length) {  // read only up to window size even if current chunk is larger
                        for (uint16_t ib = 0; ib < _bands.count(); ++ib) {
                            if (l_chunks[i]->empty()) {
                                cur_ts[ib * cur_ts_length + tsidx] = NAN;
                            }
                            else {
                                cur_ts[ib * cur_ts_length + tsidx] = ((double*)(l_chunks[i]->buf()))[_band_idx_in[ib] * (l_chunks[i]->size()[1] * l_chunks[i]->size()[2] * l_chunks[i]->size()[3]) + (l_chunks[i]->size()[1] - 1 - ic) * (l_chunks[i]->size()[2] * l_chunks[i]->size()[3]) + ixy];
                            }
                        }
                        tsidx--;
                    }
                }
            }
        }

        // fill values from current chunk
        tsidx = _win_size_l;
        for (uint32_t ic = 0; ic < this_chunk->size()[1]; ++ic) {
            if (tsidx >= 0 && tsidx < (int32_t)cur_ts_length) {  // read only up to window size even if current chunk is larger
                for (uint16_t ib = 0; ib < _bands.count(); ++ib) {
                    if (this_chunk->empty()) {
                        cur_ts[ib * cur_ts_length + tsidx] = NAN;
                    }
                    else {
                        cur_ts[ib * cur_ts_length + tsidx] = ((double*)(this_chunk->buf()))[_band_idx_in[ib] * (this_chunk->size()[1] * this_chunk->size()[2] * this_chunk->size()[3]) +
                                                                                            ic * (this_chunk->size()[2] * this_chunk->size()[3]) + ixy];
                    }
                    tsidx++;
                }
            }
        }

        // fill values from r chunks
        for (uint16_t i = 0; i < chunk_count_r; ++i) {
            if (i >= r_chunks.size()) {
                // fill NA
                for (uint32_t ic = 0; ic < _in_cube->chunk_size()[0]; ++ic) {
                    if (tsidx >= 0 && tsidx < (int32_t)cur_ts_length) {
                        for (uint16_t ib = 0; ib < _bands.count(); ++ib) {
                            cur_ts[ib * cur_ts_length + tsidx] = NAN;
                        }
                        tsidx++;
                    }
                }

            } else {
                for (uint32_t ic = 0; ic < r_chunks[i]->size()[1]; ++ic) {
                    if (tsidx >= 0 && tsidx < (int32_t)cur_ts_length) {  // read only up to window size even if current chunk is larger
                        for (uint16_t ib = 0; ib < _bands.count(); ++ib) {
                            if (r_chunks[i]->empty()) {
                                cur_ts[ib * cur_ts_length + tsidx] = NAN;
                            }
                            else {
                                cur_ts[ib * cur_ts_length + tsidx] = ((double*)(r_chunks[i]->buf()))[_band_idx_in[ib] * (r_chunks[i]->size()[1] * r_chunks[i]->size()[2] * r_chunks[i]->size()[3]) + ic * (r_chunks[i]->size()[2] * r_chunks[i]->size()[3]) + ixy];
                            }
                        }
                        tsidx++;
                    }
                }
            }
        }

        // compute new value
        for (uint32_t ic = 0; ic < size_tyx[0]; ++ic) {
            // call function over all windows
            for (uint16_t ib = 0; ib < _bands.count(); ++ib) {
                ((double*)out->buf())[ib * (size_tyx[0] * size_tyx[1] * size_tyx[2]) + ic * size_tyx[1] * size_tyx[2] + ixy] = _f[ib](&(cur_ts[ic]), _win_size_l + 1 + _win_size_r);
            }
        }
    }
    std::free(cur_ts);

    // check if chunk is completely NAN and if yes, return empty chunk
    if (out->all_nan()) {
        out = std::make_shared<chunk_data>();
    }
    return out;
}

}  // namespace gdalcubes
