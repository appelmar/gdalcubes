#include "fill_time.h"

#include <unordered_map>

namespace gdalcubes {

std::shared_ptr<chunk_data> fill_time_cube::read_chunk(chunkid_t id) {
    GCBS_TRACE("fill_time_cube::read_chunk(" + std::to_string(id) + ")");
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

    //std::shared_ptr<chunk_data> this_chunk = _in_cube->read_chunk(id);
    //std::vector<std::shared_ptr<chunk_data>> l_chunks;
    //std::vector<std::shared_ptr<chunk_data>> r_chunks;

    std::unordered_map<chunkid_t, std::shared_ptr<chunk_data>> in_chunks;
    in_chunks.insert(std::pair<chunkid_t, std::shared_ptr<chunk_data>>(id, _in_cube->read_chunk(id)));

    if (in_chunks[id]->empty()) {  // if input chunk is empty, fill with NANs
        in_chunks[id]->size(size_btyx);
        in_chunks[id]->buf(std::calloc(size_btyx[0] * size_btyx[1] * size_btyx[2] * size_btyx[3], sizeof(double)));
        std::fill((double*)(in_chunks[id]->buf()), ((double*)(in_chunks[id]->buf())) + size_btyx[0] * size_btyx[1] * size_btyx[2] * size_btyx[3], NAN);
    }

    // iterate over all pixel time series
    for (uint32_t ixy = 0; ixy < size_btyx[2] * size_btyx[3]; ++ixy) {
        // and all bands...
        for (uint32_t ib = 0; ib < size_btyx[0]; ++ib) {
            int32_t cur_t = 0;
            int32_t next_t = 0;  // points to next non NAN cell in current  time series
            int32_t prev_t = 0;  // points to next non NAN cell in current  time series

            int32_t next_chunk = id;
            int32_t prev_chunk = id;

            // 1. INITIALIZE POINTERS
            // if value of current pixel and current band at cur_t is nan:
            // - forwards iterate next_t until non-NAN found (potentially exceeding chunk) or end of time series
            // - backwars iterate prev_t until non_NAN fount (potentially exceeding chunk) or end of time series

            if (std::isnan(((double*)in_chunks[id]->buf())[ib * size_btyx[1] * size_btyx[2] * size_btyx[3] + cur_t * size_btyx[2] * size_btyx[3] + ixy])) {
                // backwards iterate prev_t
                bool found = false;
                prev_chunk -= _in_cube->count_chunks_x() * _in_cube->count_chunks_y();
                while (prev_chunk >= 0 && !found) {
                    // load chunk (only if needed)
                    if (in_chunks.find(prev_chunk) == in_chunks.end()) {
                        in_chunks.insert(std::pair<chunkid_t, std::shared_ptr<chunk_data>>(prev_chunk, _in_cube->read_chunk(prev_chunk)));
                    }
                    if (!in_chunks[prev_chunk]->empty()) {
                        prev_t = _in_cube->chunk_size()[0] - 1;
                        chunk_size_tyx cs = _in_cube->chunk_size(prev_chunk);
                        while (prev_t >= 0 && !found) {
                            // if not nan break both loops
                            if (!std::isnan(((double*)in_chunks[prev_chunk]->buf())[ib * cs[0] * cs[1] * cs[2] +
                                                                                    prev_t * cs[1] * cs[2] + ixy])) {
                                found = true;
                            } else {
                                --prev_t;
                            }
                        }
                    }
                    if (!found)
                        prev_chunk -= _in_cube->count_chunks_x() * _in_cube->count_chunks_y();
                }
                if (!found) {
                    // reset, there is no non NAN value before
                    prev_chunk = id;
                    prev_t = 0;
                }

                // forwards iterate next_t
                found = false;
                ++next_t;
                while (next_chunk < (int32_t)_in_cube->count_chunks() && !found) {
                    // load chunk (only if needed)
                    if (in_chunks.find(next_chunk) == in_chunks.end()) {
                        in_chunks.insert(std::pair<chunkid_t, std::shared_ptr<chunk_data>>(next_chunk, _in_cube->read_chunk(next_chunk)));
                    }
                    if (!in_chunks[next_chunk]->empty()) {
                        chunk_size_tyx cs = _in_cube->chunk_size(next_chunk);
                        while (next_t < (int32_t)cs[0] && !found) {
                            // if not nan break both loops
                            if (!std::isnan(((double*)in_chunks[next_chunk]->buf())[ib * cs[0] * cs[1] * cs[2] +
                                                                                    next_t * cs[1] * cs[2] + ixy])) {
                                found = true;
                            } else {
                                ++next_t;
                            }
                        }
                    }
                    if (!found) {
                        next_t = 0;
                        next_chunk += _in_cube->count_chunks_x() * _in_cube->count_chunks_y();
                    }
                }
                if (!found) {
                    // reset, there is no non NAN value after
                    next_chunk = id;
                    next_t = 0;
                }
            }

            // next_t and prev_t should now point to the next data point (or NAN if end of time series

            // 2. FILL VALUES
            // iterate cur_t from 0 to end of chunk (size_btyx[1]-1)
            // if value of current pixel is NAN, fill based on cur_t, next_t, prev_t and _method
            // else, set prev_t = cur_t, forward iterate next_t
            // ++ cur_t

            while (cur_t < (int32_t)size_btyx[1]) {
                if (std::isnan(((double*)in_chunks[id]->buf())[ib * size_btyx[1] * size_btyx[2] * size_btyx[3] + cur_t * size_btyx[2] * size_btyx[3] + ixy])) {
                    // fill based on cur_t, next_t, prev_t and _method
                    chunk_size_tyx cs = _in_cube->chunk_size(prev_chunk);
                    double v0 = ((double*)in_chunks[prev_chunk]->buf())[ib * cs[0] * cs[1] * cs[2] + prev_t * cs[1] * cs[2] + ixy];
                    cs = _in_cube->chunk_size(next_chunk);
                    double v1 = ((double*)in_chunks[next_chunk]->buf())[ib * cs[0] * cs[1] * cs[2] + next_t * cs[1] * cs[2] + ixy];

                    int32_t prev_dist = ((id - prev_chunk) / (_in_cube->count_chunks_x() * _in_cube->count_chunks_y())) * _in_cube->chunk_size()[0] + cur_t - prev_t;
                    int32_t next_dist = ((next_chunk - id) / (_in_cube->count_chunks_x() * _in_cube->count_chunks_y())) * _in_cube->chunk_size()[0] - cur_t + next_t;

                    double* res = &((double*)out->buf())[ib * size_btyx[1] * size_btyx[2] * size_btyx[3] + cur_t * size_btyx[2] * size_btyx[3] + ixy];

                    if (_method == "linear") {  // linear interpolation
                        if (std::isnan(v0) && std::isnan(v1)) {
                            *res = NAN;
                        } else if (std::isnan(v0)) {
                            *res = v1;
                        } else if (std::isnan(v1)) {
                            *res = v0;
                        } else {
                            *res = v0 * ((double)next_dist / ((double)prev_dist + (double)next_dist)) + v1 * ((double)prev_dist / ((double)prev_dist + (double)next_dist));
                        }
                    } else if (_method == "locf") {
                        *res = v0;
                    } else if (_method == "nocb") {
                        *res = v1;
                    } else {  // nearest
                        if (std::isnan(v0) && std::isnan(v1)) {
                            *res = NAN;
                        } else if (std::isnan(v0)) {
                            *res = v1;
                        } else if (std::isnan(v1)) {
                            *res = v0;
                        } else {
                            *res = (prev_dist <= next_dist) ? v0 : v1;
                        }
                    }

                } else {
                    ((double*)out->buf())[ib * size_btyx[1] * size_btyx[2] * size_btyx[3] + cur_t * size_btyx[2] * size_btyx[3] + ixy] = ((double*)in_chunks[id]->buf())[ib * size_btyx[1] * size_btyx[2] * size_btyx[3] + cur_t * size_btyx[2] * size_btyx[3] + ixy];

                    prev_t = cur_t;
                    prev_chunk = id;

                    // forwards iterate next_t (if needed)
                    bool found = false;
                    ++next_t;
                    while (next_chunk < (int32_t)_in_cube->count_chunks() && !found) {
                        // load chunk (only if needed)
                        if (in_chunks.find(next_chunk) == in_chunks.end()) {
                            in_chunks.insert(std::pair<chunkid_t, std::shared_ptr<chunk_data>>(next_chunk, _in_cube->read_chunk(next_chunk)));
                        }
                        if (!in_chunks[next_chunk]->empty()) {
                            chunk_size_tyx cs = _in_cube->chunk_size(next_chunk);

                            while (next_t < (int32_t)cs[0] && !found) {
                                // if not nan break both loops
                                if (!std::isnan(((double*)in_chunks[next_chunk]->buf())[ib * cs[0] * cs[1] * cs[2] +
                                                                                        next_t * cs[1] * cs[2] +
                                                                                        ixy])) {
                                    found = true;
                                } else {
                                    ++next_t;
                                }
                            }
                        }
                        if (!found) {
                            next_t = 0;
                            next_chunk += _in_cube->count_chunks_x() * _in_cube->count_chunks_y();
                        }
                    }
                    if (!found) {
                        // reset, there is no non NAN value after
                        next_chunk = id;
                        next_t = cur_t + 1;  // This should still be OK in case the temporal chunk size == 1, because overall while loop iterates only once
                    }
                }
                ++cur_t;
            }
        }
    }

    return out;
}

}  // namespace gdalcubes
