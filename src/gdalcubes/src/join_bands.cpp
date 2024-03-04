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

#include "join_bands.h"
#include <cstring>

namespace gdalcubes {

std::shared_ptr<chunk_data> join_bands_cube::read_chunk(chunkid_t id) {
    GCBS_TRACE("join_bands_cube::read_chunk(" + std::to_string(id) + ")");
    std::shared_ptr<chunk_data> out = std::make_shared<chunk_data>();
    if (id >= count_chunks())
        return out;  // chunk is outside of the view, we don't need to read anything.

    coords_nd<uint32_t, 3> size_tyx = chunk_size(id);
    coords_nd<uint32_t, 4> size_btyx = {_bands.count(), size_tyx[0], size_tyx[1], size_tyx[2]};
    out->size(size_btyx);

    //
    //    std::shared_ptr<chunk_data> dat_A = _in_A->read_chunk(id);
    //    std::shared_ptr<chunk_data> dat_B = _in_B->read_chunk(id);
    //    if (dat_A->empty() && dat_B->empty()) {
    //        return out;
    //    }

    // Fill buffers accordingly
    out->buf(std::calloc(size_btyx[0] * size_btyx[1] * size_btyx[2] * size_btyx[3], sizeof(double)));
    double *begin = (double *)out->buf();
    double *end = ((double *)out->buf()) + size_btyx[0] * size_btyx[1] * size_btyx[2] * size_btyx[3];
    std::fill(begin, end, NAN);

    uint32_t offset = 0;
    bool allempty = true;
    for (uint16_t i = 0; i < _in.size(); ++i) {
        std::shared_ptr<chunk_data> dat = _in[i]->read_chunk(id);
        // propagate chunk status
        if (dat->status() == chunk_data::chunk_status::ERROR) {
            out->set_status(chunk_data::chunk_status::ERROR);
        }
        else if (dat->status() == chunk_data::chunk_status::INCOMPLETE && out->status() != chunk_data::chunk_status::ERROR) {
            out->set_status(chunk_data::chunk_status::INCOMPLETE);
        }

        if (!dat->empty()) {
            allempty = false;
            std::memcpy(((double *)out->buf()) + offset, ((double *)dat->buf()), dat->size()[0] * dat->size()[1] * dat->size()[2] * dat->size()[3] * sizeof(double));
        }
        offset += _in[i]->size_bands() * size_tyx[0] * size_tyx[1] * size_tyx[2];
    }
    if (allempty) {
        auto s = out->status();
        // propagated empty chunk if all input chunks are emtpy
        out = std::make_shared<chunk_data>();
        out->set_status(s);
    }
    return out;
}

}  // namespace gdalcubes
