#include "crop.h"

namespace gdalcubes {

std::shared_ptr<chunk_data> crop_cube::read_chunk(chunkid_t id) {
    GCBS_TRACE("crop_cube::read_chunk(" + std::to_string(id) + ")");
    std::shared_ptr<chunk_data> out = std::make_shared<chunk_data>();
    if (id >= count_chunks())
        return out;  // chunk is outside of the view, we don't need to read anything.

    coords_nd<uint32_t, 3> size_tyx = chunk_size(id);
    coords_nd<uint32_t, 4> size_btyx = {uint32_t(_bands.count()), size_tyx[0], size_tyx[1], size_tyx[2]};

    auto output_chunk_coords = chunk_coords_from_id(id);
    bool chunk_is_initialized = false;

    //  integer coordinates of current chunk with regard to the input cube
    int32_t abs_low_x = output_chunk_coords[2]  * _chunk_size[2] + _x_min;
    int32_t abs_high_x = abs_low_x + size_tyx[2] - 1; // TODO: check and test for different chunk sizes
    if (abs_low_x >= (int32_t)_in_cube->size_x() ||
        abs_high_x < 0) {
        return out; // completely outside input cube
    }
    if (abs_low_x < 0) abs_low_x = 0;
    if (abs_high_x >= (int32_t)_in_cube->size_x()) abs_high_x = _in_cube->size_x() - 1;

    int32_t abs_low_y = output_chunk_coords[1]  * _chunk_size[1] + _y_min;
    int32_t abs_high_y = abs_low_y + size_tyx[1] - 1; // TODO: check and test for different chunk sizes
    if (abs_low_y >= (int32_t)_in_cube->size_y() ||
        abs_high_y < 0) {
        return out; // completely outside input cube
    }
    if (abs_low_y < 0) abs_low_y = 0;
    if (abs_high_y >= (int32_t)_in_cube->size_y()) abs_high_y = _in_cube->size_y() - 1;

    int32_t abs_low_t = output_chunk_coords[0]  * _chunk_size[0] + _t_min;
    int32_t abs_high_t = abs_low_t + size_tyx[0] - 1; // TODO: check and test for different chunk sizes
    if (abs_low_t >= (int32_t)_in_cube->size_t() ||
        abs_high_t < 0) {
        return out; // completely outside input cube
    }
    if (abs_low_t < 0) abs_low_t = 0;
    if (abs_high_t >= (int32_t)_in_cube->size_t()) abs_high_t = _in_cube->size_t() - 1;

    chunk_coordinate_tyx input_chunk_coords_low = {abs_low_t /  _in_cube->chunk_size()[0],
                                                   abs_low_y /  _in_cube->chunk_size()[1],
                                                   abs_low_x /  _in_cube->chunk_size()[2]};

    chunk_coordinate_tyx input_chunk_coords_high = {abs_high_t /  _in_cube->chunk_size()[0],
                                                    abs_high_y /  _in_cube->chunk_size()[1],
                                                    abs_high_x /  _in_cube->chunk_size()[2]};



    for (uint16_t ch_t = input_chunk_coords_low[0]; ch_t <= input_chunk_coords_high[0]; ++ch_t) {
        for (uint16_t ch_y = input_chunk_coords_low[1]; ch_y <= input_chunk_coords_high[1]; ++ch_y) {
            for (uint16_t ch_x = input_chunk_coords_low[2]; ch_x <= input_chunk_coords_high[2]; ++ch_x) {
                chunkid_t input_chunk_id = _in_cube->chunk_id_from_coords({ch_t, ch_y, ch_x});
                std::shared_ptr<chunk_data> in_chunk = _in_cube->read_chunk(input_chunk_id);

                if (in_chunk->empty()) {
                    continue;
                }
                if (!chunk_is_initialized) {
                    out->size(size_btyx);
                    // Fill buffers accordingly
                    out->buf(std::calloc(size_btyx[0] * size_btyx[1] * size_btyx[2] * size_btyx[3], sizeof(double)));
                    double* begin = (double*)out->buf();
                    double* end = ((double*)out->buf()) + size_btyx[0] * size_btyx[1] * size_btyx[2] * size_btyx[3];
                    std::fill(begin, end, NAN);
                    chunk_is_initialized = true;
                }

                // Find out which part of the chunk is needed
                int32_t in_abs_low_x = ch_x  * _in_cube->chunk_size()[2];
                int32_t in_abs_high_x = in_abs_low_x + in_chunk->size()[3] - 1; // TODO: check and test for different chunk sizes
                int32_t in_abs_low_y = ch_y  *_in_cube->chunk_size()[1];
                int32_t in_abs_high_y = in_abs_low_y + in_chunk->size()[2]- 1; // TODO: check and test for different chunk sizes
                int32_t in_abs_low_t = ch_t  * _in_cube->chunk_size()[0];
                int32_t in_abs_high_t = in_abs_low_t + in_chunk->size()[1] - 1; // TODO: check and test for different chunk sizes

                int32_t start_x = std::max(in_abs_low_x, abs_low_x);
                int32_t end_x = std::min(in_abs_high_x, abs_high_x);
                int32_t start_y = std::max(in_abs_low_y, abs_low_y);
                int32_t end_y = std::min(in_abs_high_y, abs_high_y);
                int32_t start_t = std::max(in_abs_low_t, abs_low_t);
                int32_t end_t = std::min(in_abs_high_t, abs_high_t);

                for (int32_t ib = 0; ib<(int32_t)size_btyx[0]; ++ib) {
                    for (int32_t it=start_t; it <= end_t; ++it) {
                        for (int32_t iy=start_y; iy <= end_y; ++iy) {
                            for (int32_t ix=start_x; ix <= end_x; ++ix) {
                                int32_t out_t = (it - _t_min) % (int32_t)chunk_size()[0];
                                int32_t out_y = (iy - _y_min) % (int32_t)chunk_size()[1];
                                int32_t out_x = (ix - _x_min) % (int32_t)chunk_size()[2];
                                int32_t in_t =  it  % (int32_t)_in_cube->chunk_size()[0];
                                int32_t in_y =  iy  % (int32_t)_in_cube->chunk_size()[1];
                                int32_t in_x =  ix  % (int32_t)_in_cube->chunk_size()[2];

                                ((double*)out->buf())[ib * size_btyx[1] * size_btyx[2] * size_btyx[3] + out_t * size_btyx[2] * size_btyx[3] + out_y * size_btyx[3] + out_x] =((double*)in_chunk->buf())[ib * in_chunk->size()[1] * in_chunk->size()[2] * in_chunk->size()[3] + in_t * in_chunk->size()[2] * in_chunk->size()[3] + in_y *  in_chunk->size()[3] + in_x];
                            }
                        }
                    }
                }
            }
        }
    }

    return out;
}

}  // namespace gdalcubes
