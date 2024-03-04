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

#include "select_bands.h"

namespace gdalcubes {

std::shared_ptr<chunk_data> select_bands_cube::read_chunk(chunkid_t id) {
    GCBS_TRACE("select_bands::read_chunk(" + std::to_string(id) + ")");
    if (id >= count_chunks())
        return  std::make_shared<chunk_data>();  // chunk is outside of the view, we don't need to read anything.

    // if input cube is image_collection_cube, delegate (since in->select_bands has been called in the cosntructor)
    if (_defer_to_input_cube) {
        return _in_cube->read_chunk(id);
    }

    std::shared_ptr<chunk_data> in = _in_cube->read_chunk(id);

    // propagate chunk status

    std::shared_ptr<chunk_data> out = std::make_shared<chunk_data>();
    out->set_status(in->status());  // propagate chunk status
    if (in->empty()) {
        return out;
    }

    // Fill buffers accordingly

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

}  // namespace gdalcubes
