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

#include "stream.h"
#include "external/tiny-process-library/process.hpp"

std::shared_ptr<chunk_data> stream_cube::read_chunk(chunkid_t id) {
    GCBS_DEBUG("stream_cube::read_chunk(" + std::to_string(id) + ")");
    std::shared_ptr<chunk_data> out = std::make_shared<chunk_data>();
    if (id < 0 || id >= count_chunks()) {
        // chunk is outside of the cube, we don't need to read anything.
        GCBS_WARN("Chunk id " + std::to_string(id) + " is out of range");
        return out;
    }

    out = stream_chunk_stdin(_in_cube->read_chunk(id));
    if (out->empty()) {
        GCBS_DEBUG("Streaming returned empty chunk " + std::to_string(id));
    }
    return out;
}

std::shared_ptr<chunk_data> stream_cube::stream_chunk_stdin(std::shared_ptr<chunk_data> data) {
    std::shared_ptr<chunk_data> out = std::make_shared<chunk_data>();

    int size[] = {(int)data->size()[0], (int)data->size()[1], (int)data->size()[2], (int)data->size()[3]};
    if (size[0] * size[1] * size[2] * size[3] == 0) {
        return out;
    }

    // 1. new process with env "GDALCUBES_STREAMING=1"
    std::string errstr;
    // TODO: must _cmd be splitted by arguments as vector?

    uint32_t databytes_read = 0;
    TinyProcessLib::Process process(_cmd, "", {{"GDALCUBES_STREAMING", "1"}}, [out, &databytes_read](const char *bytes, std::size_t n) {

        if (databytes_read == 0) {
            // Assumption is that at least 4 integers with chunk size area always contained in the first call of this function
            if (n >= 4 * sizeof(int)) {
                chunk_size_btyx out_size = {(uint32_t)(((int *)bytes)[0]), (uint32_t)(((int *)bytes)[1]),
                                            (uint32_t)(((int *)bytes)[2]), (uint32_t)(((int *)bytes)[3])};
                out->size(out_size);
                // Fill buffers accordingly
                out->buf(std::calloc(out_size[0] * out_size[1] * out_size[2] * out_size[3], sizeof(double)));
                memcpy(out->buf(), bytes + (4 * sizeof(int)), n -  4 * sizeof(int));
                databytes_read += n - (4 * sizeof(int));
            } else {
                GCBS_WARN("Cannot read streaming result, returning empty chunk");
                databytes_read = std::numeric_limits<uint32_t>::max(); // prevent further calls to write anything into the chunk buffer
            }
        } else {
            // buffer overflow check
            if ((char*)(out->buf()) + databytes_read + n <= (char*)(out->buf()) + out->size()[0] * out->size()[1] * out->size()[2] * out->size()[3] * sizeof(double)) {
                memcpy((char*)(out->buf()) + databytes_read, bytes, n);
                databytes_read += n;
            }
        } }, [&errstr, this](const char *bytes, std::size_t n) {
    errstr = std::string(bytes, n);
    if (_log_output == "stdout") {
        std::cout << errstr << std::endl;
    } else if (_log_output == "stderr") {
        std::cerr << errstr << std::endl;
    } else if (!_log_output.empty()) {
        std::ofstream flog(_log_output, std::ios_base::out | std::ios_base::app);
        if (flog.fail()) {
            GCBS_WARN("Failed to open file '" + _log_output + "' for writing streaming output");
        } else {
            flog << errstr;
            flog.close();
        }
    } }, true);

    // Write to stdin
    std::string proj = _in_cube->st_reference()->proj();
    process.write((char *)(size), sizeof(int) * 4);
    for (uint16_t i = 0; i < _in_cube->bands().count(); ++i) {
        int str_size = _in_cube->bands().get(i).name.size();
        process.write((char *)(&str_size), sizeof(int));
        process.write(_in_cube->bands().get(i).name.c_str(), sizeof(char) * str_size);
    }
    double *dims = (double *)std::calloc(size[1] + size[2] + size[3], sizeof(double));
    for (int i = 0; i < size[1]; ++i) {
        dims[i] = (_in_cube->st_reference()->t0() + _in_cube->st_reference()->dt() * i).to_double();
    }
    for (int i = size[1]; i < size[1] + size[2]; ++i) {
        dims[i] = _in_cube->st_reference()->win().bottom + i * _in_cube->st_reference()->dy();
    }
    for (int i = size[1] + size[2]; i < size[1] + size[2] + size[3]; ++i) {
        dims[i] = _in_cube->st_reference()->win().left + i * _in_cube->st_reference()->dx();
    }
    process.write((char *)(dims), sizeof(double) * (size[1] + size[2] + size[3]));
    std::free(dims);

    int str_size = proj.size();
    process.write((char *)(&str_size), sizeof(int));
    process.write(proj.c_str(), sizeof(char) * str_size);
    process.write(((char *)(data->buf())), sizeof(double) * data->size()[0] * data->size()[1] * data->size()[2] * data->size()[3]);

    process.close_stdin();  // needed?

    auto exit_status = process.get_exit_status();
    if (exit_status != 0) {
        GCBS_ERROR("Child process failed with exit code " + std::to_string(exit_status));
        GCBS_ERROR("Child process output: " + errstr);
        throw std::string("ERROR in stream_cube::read_chunk(): external program returned exit code " + std::to_string(exit_status));
    }
    return out;
}