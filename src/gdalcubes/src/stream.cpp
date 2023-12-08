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

#include "stream.h"

#include <stdlib.h>
#include <fstream>
#include <cstring>

#include "external/tiny-process-library/process.hpp"

namespace gdalcubes {

std::shared_ptr<chunk_data> stream_cube::read_chunk(chunkid_t id) {
    GCBS_TRACE("stream_cube::read_chunk(" + std::to_string(id) + ")");
    std::shared_ptr<chunk_data> out = std::make_shared<chunk_data>();
    if (id >= count_chunks()) {
        // chunk is outside of the cube, we don't need to read anything.
        GCBS_WARN("Chunk id " + std::to_string(id) + " is out of range");
        return out;
    }
    out = stream_chunk_file(_in_cube->read_chunk(id), id);
    if (out->empty()) {
        GCBS_DEBUG("Streaming returned empty chunk " + std::to_string(id));
    }
    return out;
}


std::shared_ptr<chunk_data> stream_cube::stream_chunk_file(std::shared_ptr<chunk_data> data, chunkid_t id) {
    std::shared_ptr<chunk_data> out = std::make_shared<chunk_data>();

    int size[] = {(int)data->size()[0], (int)data->size()[1], (int)data->size()[2], (int)data->size()[3]};
    if (size[0] * size[1] * size[2] * size[3] == 0) {
        return out;
    }

    // generate in and out filename
    std::string f_in = filesystem::join(config::instance()->get_streaming_dir(), utils::generate_unique_filename(12, ".stream_" + std::to_string(id) + "_", "_in"));
    std::string f_out = filesystem::join(config::instance()->get_streaming_dir(), utils::generate_unique_filename(12, ".stream_" + std::to_string(id) + "_", "_out"));

    std::string errstr;  // capture error string

    // write input data
    std::ofstream f_in_stream(f_in, std::ios::out | std::ios::binary | std::ios::trunc);
    if (!f_in_stream.is_open()) {
        GCBS_ERROR("Cannot write streaming input data to file '" + f_in + "'");
        throw std::string("ERROR in stream_cube::stream_chunk_file(): cannot write streaming input data to file '" + f_in + "'");
    }

    std::string proj = _in_cube->st_reference()->srs();
    f_in_stream.write((char *)(size), sizeof(int) * 4);
    for (uint16_t i = 0; i < _in_cube->bands().count(); ++i) {
        int str_size = _in_cube->bands().get(i).name.size();
        f_in_stream.write((char *)(&str_size), sizeof(int));
        f_in_stream.write(_in_cube->bands().get(i).name.c_str(), sizeof(char) * str_size);
    }

    double *dims = (double *)std::calloc(size[1] + size[2] + size[3], sizeof(double));
    int i = 0;
    for (int it = 0; it < size[1]; ++it) {
        dims[i] = (_in_cube->st_reference()->datetime_at_index(it + _in_cube->chunk_size()[0] * _in_cube->chunk_limits(id).low[0])).to_double();
        ++i;
    }
    bounds_st cextent = this->bounds_from_chunk(id);  // implemented in derived classes
    for (int iy = 0; iy < size[2]; ++iy) {
        dims[i] = cextent.s.top - (iy + 0.5) * st_reference()->dy();  // cell center
        ++i;
    }
    for (int ix = 0; ix < size[3]; ++ix) {
        dims[i] = cextent.s.left + (ix + 0.5) * st_reference()->dx();
        ++i;
    }

    f_in_stream.write((char *)(dims), sizeof(double) * (size[1] + size[2] + size[3]));
    std::free(dims);

    int str_size = proj.size();
    f_in_stream.write((char *)(&str_size), sizeof(int));
    f_in_stream.write(proj.c_str(), sizeof(char) * str_size);
    f_in_stream.write(((char *)(data->buf())), sizeof(double) * data->size()[0] * data->size()[1] * data->size()[2] * data->size()[3]);
    f_in_stream.close();

    /* setenv / _putenv is not thread-safe, we need to get a mutex until the child process has been started. */
    static std::mutex mtx;
    mtx.lock();
    utils::env::instance().set({
        {"GDALCUBES_STREAMING", "1"},
        {"GDALCUBES_STREAMING_CHUNK_ID", std::to_string(id)},
        {"GDALCUBES_STREAMING_FILE_IN", f_in},
        {"GDALCUBES_STREAMING_FILE_OUT", f_out}});

    // start process
    TinyProcessLib::Process process(
        _cmd, "", [](const char *bytes, std::size_t n) {}, [&errstr](const char *bytes, std::size_t n) {
        errstr = std::string(bytes, n);
        GCBS_DEBUG(errstr); }, false);
    utils::env::instance().unset_all();
    mtx.unlock();
    auto exit_status = process.get_exit_status();
    filesystem::remove(f_in);
    if (exit_status != 0) {
        GCBS_ERROR("Child process failed with exit code " + std::to_string(exit_status));
        GCBS_ERROR("Child process output: " + errstr);
        if (filesystem::exists(f_out)) {
            filesystem::remove(f_out);
        }
        throw std::string("ERROR in stream_cube::read_chunk(): external program returned exit code " + std::to_string(exit_status));
    }

    // read output data
    std::ifstream f_out_stream(f_out, std::ios::in | std::ios::binary);
    if (!f_out_stream.is_open()) {
        GCBS_ERROR("Cannot read streaming output data from file '" + f_out + "'");
        throw std::string("ERROR in stream_cube::stream_chunk_file(): cannot read streaming output data from file '" + f_out + "'");
    }

    f_out_stream.seekg(0, f_out_stream.end);
    int length = f_out_stream.tellg();
    f_out_stream.seekg(0, f_out_stream.beg);
    char *buffer = (char *)std::calloc(length, sizeof(char));
    f_out_stream.read(buffer, length);
    f_out_stream.close();

    chunk_size_btyx out_size = {(uint32_t)(((int *)buffer)[0]), (uint32_t)(((int *)buffer)[1]),
                                (uint32_t)(((int *)buffer)[2]), (uint32_t)(((int *)buffer)[3])};
    out->size(out_size);
    out->buf(std::calloc(out_size[0] * out_size[1] * out_size[2] * out_size[3], sizeof(double)));
    std::memcpy(out->buf(), buffer + (4 * sizeof(int)), length - 4 * sizeof(int));
    std::free(buffer);

    if (filesystem::exists(f_out)) {
        filesystem::remove(f_out);
    }

    return out;
}

}  // namespace gdalcubes
