#include "stream_reduce_space.h"

#include <fstream>
#include <cstring>

#include "external/tiny-process-library/process.hpp"

namespace gdalcubes {

std::shared_ptr<chunk_data> stream_reduce_space_cube::read_chunk(chunkid_t id) {
    GCBS_TRACE("stream_reduce_space_cube::read_chunk(" + std::to_string(id) + ")");
    std::shared_ptr<chunk_data> out = std::make_shared<chunk_data>();
    if (id >= count_chunks())
        return out;  // chunk is outside of the view, we don't need to read anything.

    coords_nd<uint32_t, 3> size_tyx = chunk_size(id);
    coords_nd<uint32_t, 4> size_btyx = {_nbands, size_tyx[0], 1, 1};
    out->size(size_btyx);



    // 1. read everything to input buffer (for first version, can be memory-intensive, same as rechunk_merge_time)
    std::shared_ptr<chunk_data> inbuf = std::make_shared<chunk_data>();
    coords_nd<uint32_t, 4> in_size_btyx = {uint32_t(_in_cube->size_bands()), size_tyx[0], _in_cube->size_y(),
                                           _in_cube->size_x()};
    inbuf->size(in_size_btyx);

    uint32_t nchunks_in_space = _in_cube->count_chunks_x() * _in_cube->count_chunks_x();
    bool empty = true;
    bool initialized = false;
    for (chunkid_t i = 0; i < nchunks_in_space; ++i) {
        uint32_t in_chunk_id = id * nchunks_in_space + i;
        auto in_chunk_coords = _in_cube->chunk_coords_from_id(in_chunk_id);
        std::shared_ptr<chunk_data> x = _in_cube->read_chunk(in_chunk_id);

        if (!x->empty()) {
            if (!initialized) {
                // Fill buffers with NAN
                out->buf(std::calloc(size_btyx[0] * size_btyx[1] * size_btyx[2] * size_btyx[3], sizeof(double)));
                double *begin = (double *)out->buf();
                double *end = ((double *)out->buf()) + size_btyx[0] * size_btyx[1] * size_btyx[2] * size_btyx[3];
                std::fill(begin, end, NAN);

                inbuf->buf(std::calloc(in_size_btyx[0] * in_size_btyx[1] * in_size_btyx[2] * in_size_btyx[3], sizeof(double)));
                double *inbegin = (double *)inbuf->buf();
                double *inend =
                    ((double *)inbuf->buf()) + in_size_btyx[0] * in_size_btyx[1] * in_size_btyx[2] * in_size_btyx[3];
                std::fill(inbegin, inend, NAN);

                initialized = true;
            }

            for (uint16_t ib = 0; ib < x->size()[0]; ++ib) {
                for (uint32_t it = 0; it < x->size()[1]; ++it) {
                    for (uint32_t iy = 0; iy < x->size()[2]; ++iy) {
                        for (uint32_t ix = 0; ix < x->size()[3]; ++ix) {
                            ((double *)inbuf->buf())[ib * x->size()[1] * _in_cube->size_y() * _in_cube->size_x() +
                                                     (it * _in_cube->size_y() * _in_cube->size_x()) +
                                                     (iy + in_chunk_coords[1] * _in_cube->chunk_size()[1]) * _in_cube->size_x() +
                                                     (ix + in_chunk_coords[2] * _in_cube->chunk_size()[2])] =
                                ((double *)x->buf())[ib * x->size()[1] * x->size()[2] * x->size()[3] +
                                                     it * x->size()[2] * x->size()[3] +
                                                     iy * x->size()[3] + ix];
                        }
                    }
                }
            }
            empty = false;
        }
    }
    // check if inbuf is completely empty and if yes, avoid streaming at all and return empty chunk
    if (empty) {
        return std::make_shared<chunk_data>();
    }

    // generate in and out filename
    std::string f_in = filesystem::join(config::instance()->get_streaming_dir(),
                                        utils::generate_unique_filename(12, ".stream_" + std::to_string(id) + "_", "_in"));
    std::string f_out = filesystem::join(config::instance()->get_streaming_dir(),
                                         utils::generate_unique_filename(12, ".stream_" + std::to_string(id) + "_", "_out"));

    std::string errstr;  // capture error string

    // write input data
    std::ofstream f_in_stream(f_in, std::ios::out | std::ios::binary | std::ios::trunc);
    if (!f_in_stream.is_open()) {
        GCBS_ERROR("Cannot write streaming input data to file '" + f_in + "'");
        throw std::string(
            "ERROR in stream_reduce_time_cube::stream_chunk_file(): cannot write streaming input data to file '" + f_in +
            "'");
    }

    int size[] = {(int)in_size_btyx[0], (int)in_size_btyx[1], (int)in_size_btyx[2], (int)in_size_btyx[3]};
    if (size[0] * size[1] * size[2] * size[3] == 0) {
        return out;
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
        dims[i] = (_in_cube->st_reference()->datetime_at_index(it + id * chunk_size()[0])).to_double();
        ++i;
    }

    for (int iy = 0; iy < size[2]; ++iy) {
        dims[i] = _in_cube->st_reference()->top() - (iy + 0.5) * st_reference()->dy();  // cell center
        ++i;
    }
    for (int ix = 0; ix < size[3]; ++ix) {
        dims[i] = _in_cube->st_reference()->left() + (ix + 0.5) * st_reference()->dx();  // cell center
        ++i;
    }
    f_in_stream.write((char *)(dims), sizeof(double) * (size[1] + size[2] + size[3]));
    std::free(dims);

    int str_size = proj.size();
    f_in_stream.write((char *)(&str_size), sizeof(int));
    f_in_stream.write(proj.c_str(), sizeof(char) * str_size);
    f_in_stream.write(((char *)(inbuf->buf())),
                      sizeof(double) * inbuf->size()[0] * inbuf->size()[1] * inbuf->size()[2] * inbuf->size()[3]);
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
    TinyProcessLib::Config pconf;
    pconf.show_window = TinyProcessLib::Config::ShowWindow::hide;
    TinyProcessLib::Process process(
        _cmd, "", [](const char *bytes, std::size_t n) {},
        [&errstr](const char *bytes, std::size_t n) {
            std::string s(bytes, n);
            errstr = errstr + s;
        },
        false, pconf);
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
        throw std::string("ERROR in stream_reduce_space_cube::read_chunk(): external program returned exit code " +
                          std::to_string(exit_status));
    }
    GCBS_DEBUG(errstr);

    // read output data
    std::ifstream f_out_stream(f_out, std::ios::in | std::ios::binary);
    if (!f_out_stream.is_open()) {
        GCBS_ERROR("Cannot read streaming output data from file '" + f_out + "'");
        throw std::string(
            "ERROR in stream_reduce_space_cube::read_chunk(): cannot read streaming output data from file '" +
            f_out + "'");
    }

    f_out_stream.seekg(0, f_out_stream.end);
    int length = f_out_stream.tellg();
    f_out_stream.seekg(0, f_out_stream.beg);
    char *buffer = (char *)std::calloc(length, sizeof(char));
    f_out_stream.read(buffer, length);
    f_out_stream.close();

    // Copy results to chunk buffer, at most the size of the output
    std::memcpy(out->buf(), buffer + (4 * sizeof(int)), std::min(length - 4 * sizeof(int), size_btyx[0] * size_btyx[1] * size_btyx[2] * size_btyx[3] * sizeof(double)));
    std::free(buffer);

    if (filesystem::exists(f_out)) {
        filesystem::remove(f_out);
    }

    return out;
}

}  // namespace gdalcubes