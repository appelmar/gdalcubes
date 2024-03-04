#include "slice_time.h"
#include <cstring>

namespace gdalcubes {

std::shared_ptr<chunk_data> slice_time_cube::read_chunk(chunkid_t id) {
    GCBS_TRACE("slice_time_cube::read_chunk(" + std::to_string(id) + ")");
    std::shared_ptr<chunk_data> out = std::make_shared<chunk_data>();
    if (id >= count_chunks())
        return out;  // chunk is outside of the view, we don't need to read anything.

    coords_nd<uint32_t, 3> size_tyx = chunk_size(id);
    coords_nd<uint32_t, 4> size_btyx = {uint32_t(_bands.count()), size_tyx[0], size_tyx[1], size_tyx[2]};

    auto output_chunk_coords = chunk_coords_from_id(id);
    auto input_chunk_coords = output_chunk_coords;
    input_chunk_coords[0] = _t_index / _in_cube->chunk_size()[0];
    chunkid_t input_chunk_id = _in_cube->chunk_id_from_coords(input_chunk_coords);

    std::shared_ptr<chunk_data> in_chunk = _in_cube->read_chunk(input_chunk_id);
    out->set_status(in_chunk->status());  // propagate chunk status
    if (in_chunk && !in_chunk->empty()) {
        out->size(size_btyx);
        // Fill buffers accordingly
        out->buf(std::calloc(size_btyx[0] * size_btyx[1] * size_btyx[2] * size_btyx[3], sizeof(double)));
        double* begin = (double*)out->buf();
        double* end = ((double*)out->buf()) + size_btyx[0] * size_btyx[1] * size_btyx[2] * size_btyx[3];
        std::fill(begin, end, NAN);

        for (uint16_t ib = 0; ib < size_btyx[0]; ++ib) {
            std::memcpy(&(((double*)out->buf())[ib * size_btyx[1] * size_btyx[2] * size_btyx[3]]),
                        &(((double*)in_chunk->buf())[ib * in_chunk->size()[1] * in_chunk->size()[2] * in_chunk->size()[3] + (_t_index % _in_cube->chunk_size()[0]) * in_chunk->size()[2] * in_chunk->size()[3]]),
                        size_btyx[2] * size_btyx[3] * sizeof(double));
        }
    }
    return out;
}

}  // namespace gdalcubes
