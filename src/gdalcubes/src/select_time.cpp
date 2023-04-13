#include "select_time.h"
#include <cstring>

namespace gdalcubes {

std::shared_ptr<chunk_data> select_time_cube::read_chunk(chunkid_t id) {
    GCBS_TRACE("select_time_cube::read_chunk(" + std::to_string(id) + ")");
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

    std::shared_ptr<chunk_data> in_chunk = nullptr;
    chunkid_t cur_input_chunk_id = 0;

    auto output_chunk_coords = chunk_coords_from_id(id);

    // TODO: what if input chunk already has irregular time dimension
    for (uint32_t it = 0; it < size_tyx[0]; ++it) {
        datetime t = std::dynamic_pointer_cast<cube_stref_labeled_time>(_st_ref)->datetime_at_index(output_chunk_coords[0] * _chunk_size[0] + it);
        uint32_t iin = _in_cube->st_reference()->index_at_datetime(t);
        if (iin >= 0 && iin < _in_cube->size_t()) {
            // COPY values
            auto input_chunk_coords = output_chunk_coords;
            input_chunk_coords[0] = iin / _in_cube->chunk_size()[0];
            chunkid_t input_chunk_id = _in_cube->chunk_id_from_coords(input_chunk_coords);
            if (!in_chunk) {
                in_chunk = _in_cube->read_chunk(input_chunk_id);
                cur_input_chunk_id = input_chunk_id;
            } else {
                if (cur_input_chunk_id != input_chunk_id) {
                    in_chunk = _in_cube->read_chunk(input_chunk_id);
                    cur_input_chunk_id = input_chunk_id;
                }
            }

            // for all bands
            if (!in_chunk->empty()) {
                for (uint16_t ib = 0; ib < size_btyx[0]; ++ib) {
                    //                assert(size_btyx[0] == in_chunk->size()[0]);
                    //                assert(size_btyx[2] == in_chunk->size()[2]);
                    //                assert(size_btyx[3] == in_chunk->size()[3]);
                    std::memcpy(&(((double*)out->buf())[ib * size_btyx[1] * size_btyx[2] * size_btyx[3] + it * size_btyx[2] * size_btyx[3]]),
                                &(((double*)in_chunk->buf())[ib * in_chunk->size()[1] * in_chunk->size()[2] * in_chunk->size()[3] + (iin % _in_cube->chunk_size()[0]) * in_chunk->size()[2] * in_chunk->size()[3]]),
                                size_btyx[2] * size_btyx[3] * sizeof(double));
                }
            }
        } else {
            GCBS_WARN("Cube does not contain date/time " + t.to_string());
        }
    }
    return out;
}

}  // namespace gdalcubes
