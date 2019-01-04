/*
   Copyright 2018 Marius Appel <marius.appel@uni-muenster.de>

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.
*/

#include "join_bands.h"

std::shared_ptr<chunk_data> join_bands_cube::read_chunk(chunkid_t id) {
    GCBS_DEBUG("join_bands_cube::read_chunk(" + std::to_string(id) + ")");
    std::shared_ptr<chunk_data> out = std::make_shared<chunk_data>();
    if (id < 0 || id >= count_chunks())
        return out;  // chunk is outside of the view, we don't need to read anything.

    coords_nd<uint32_t, 3> size_tyx = chunk_size(id);
    coords_nd<uint32_t, 4> size_btyx = {_bands.count(), size_tyx[0], size_tyx[1], size_tyx[2]};
    out->size(size_btyx);

    std::shared_ptr<chunk_data> dat_A = _in_A->read_chunk(id);
    std::shared_ptr<chunk_data> dat_B = _in_B->read_chunk(id);

    // Fill buffers accordingly
    out->buf(std::calloc(size_btyx[0] * size_btyx[1] * size_btyx[2] * size_btyx[3], sizeof(double)));
    double *begin = (double *)out->buf();
    double *end = ((double *)out->buf()) + size_btyx[0] * size_btyx[1] * size_btyx[2] * size_btyx[3];
    std::fill(begin, end, NAN);

    memcpy(((double *)out->buf()), ((double *)dat_A->buf()), dat_A->size()[0] * dat_A->size()[1] * dat_A->size()[2] * dat_A->size()[3] * sizeof(double));
    memcpy(((double *)out->buf()) + dat_A->size()[0] * dat_A->size()[1] * dat_A->size()[2] * dat_A->size()[3], ((double *)dat_B->buf()), dat_B->size()[0] * dat_B->size()[1] * dat_B->size()[2] * dat_B->size()[3] * sizeof(double));

    return out;
}