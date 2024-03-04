#include "slice_space.h"

namespace gdalcubes {

std::shared_ptr<chunk_data> slice_space_cube::read_chunk(chunkid_t id) {
    GCBS_TRACE("slice_space_cube::read_chunk(" + std::to_string(id) + ")");
    std::shared_ptr<chunk_data> out = std::make_shared<chunk_data>();
    if (id >= count_chunks())
        return out;  // chunk is outside of the view, we don't need to read anything.

    coords_nd<uint32_t, 3> size_tyx = chunk_size(id);
    coords_nd<uint32_t, 4> size_btyx = {uint32_t(_bands.count()), size_tyx[0], size_tyx[1], size_tyx[2]};

    auto output_chunk_coords = chunk_coords_from_id(id);
    auto input_chunk_coords = output_chunk_coords;
    input_chunk_coords[1] = _y_index / _in_cube->chunk_size()[1];
    input_chunk_coords[2] = _x_index / _in_cube->chunk_size()[2];
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
            for (uint32_t it = 0; it < size_btyx[1]; ++it) {
                // It is:  size_btyx[1] == size_btyx[2] == 1
                ((double*)out->buf())[ib * size_btyx[1] + it] =  ((double*)in_chunk->buf())[ib * in_chunk->size()[1] * in_chunk->size()[2] * in_chunk->size()[3] + it * in_chunk->size()[2] * in_chunk->size()[3] + (_y_index % _in_cube->chunk_size()[1]) * in_chunk->size()[3] + (_x_index % _in_cube->chunk_size()[2])];
            }
        }
    }
    return out;
}

}  // namespace gdalcubes
