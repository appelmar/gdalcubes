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

#include "select_bands.h"

std::shared_ptr<chunk_data> select_bands_cube::read_chunk(chunkid_t id) {
    GCBS_DEBUG("select_bands::read_chunk(" + std::to_string(id) + ")");
    if (id < 0 || id >= count_chunks())
        return std::shared_ptr<chunk_data>();  // chunk is outside of the view, we don't need to read anything.

    // if input cube is image_collection_cube, delegate (since in->select_bands has been called in the cosntructor)
    if (_input_is_image_collection_cube) {
        return _in_cube->read_chunk(id);
    }

    std::shared_ptr<chunk_data> in = _in_cube->read_chunk(id);

    // Fill buffers accordingly
    std::shared_ptr<chunk_data> out = std::make_shared<chunk_data>();
    out->size({_bands.count(), in->size()[1], in->size()[2], in->size()[3]});
    out->buf(std::calloc(_bands.count() * in->size()[1] * in->size()[2] * in->size()[3], sizeof(double)));

    // We do not need to fill with NAN because we can be sure that it is completeley filled from the input cube
    //double *begin = (double *)out->buf();
    //double *end = ((double *)out->buf()) + size_btyx[0] * size_btyx[1] * size_btyx[2] * size_btyx[3];
    // std::fill(begin, end, NAN);

    for (uint16_t i = 0; i < _bands.count(); ++i) {
        uint16_t orig_idx = _in_cube->bands().get_index(_bands.get(i).name);
        memcpy(((double*)out->buf()) + i * in->size()[1] * in->size()[2] * in->size()[3], ((double*)in->buf()) + orig_idx * in->size()[1] * in->size()[2] * in->size()[3], in->size()[1] * in->size()[2] * in->size()[3] * sizeof(double));
    }

    return out;
}