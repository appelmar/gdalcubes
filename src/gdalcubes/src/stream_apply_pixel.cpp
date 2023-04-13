#include "stream_apply_pixel.h"

#include <fstream>
#include <cstring>
#include "external/tiny-process-library/process.hpp"

namespace gdalcubes {

std::shared_ptr<chunk_data> stream_apply_pixel_cube::read_chunk(chunkid_t id) {
    GCBS_TRACE("stream_apply_pixel_cube::read_chunk(" + std::to_string(id) + ")");
    std::shared_ptr<chunk_data> out = std::make_shared<chunk_data>();
    if (id >= count_chunks())
        return out;  // chunk is outside of the view, we don't need to read anything.


    std::shared_ptr<chunk_data> inbuf = _in_cube->read_chunk(id);
    // check whether input chunk is empty and if yes, avoid computations
    if (inbuf->empty()) {
        return out;
    }
    coords_nd<uint32_t, 3> size_tyx = chunk_size(id);
    coords_nd<uint32_t, 4> size_btyx = {_bands.count(), size_tyx[0], size_tyx[1], size_tyx[2]};
    out->size(size_btyx);

    // Fill buffers accordingly
    out->buf(std::calloc(size_btyx[0] * size_btyx[1] * size_btyx[2] * size_btyx[3], sizeof(double)));
    double *begin = (double *)out->buf();
    double *end = ((double *)out->buf()) + size_btyx[0] * size_btyx[1] * size_btyx[2] * size_btyx[3];
    std::fill(begin, end, NAN);

    coords_nd<uint32_t, 4> in_size_btyx = {uint32_t(_in_cube->size_bands()), size_tyx[0], size_tyx[1],
                                           size_tyx[2]};

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
            "ERROR in stream_apply_pixel_cube::stream_chunk_file(): cannot write streaming input data to file '" + f_in +
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
        dims[i] = _in_cube->st_reference()->datetime_at_index(it + _in_cube->chunk_size()[0] * _in_cube->chunk_limits(id).low[0]).to_double();
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
    f_in_stream.write(((char *)(inbuf->buf())), sizeof(double) * inbuf->size()[0] * inbuf->size()[1] * inbuf->size()[2] * inbuf->size()[3]);
    f_in_stream.close();

    /* setenv / _putenv is not thread-safe, we need to get a mutex until the child process has been started. */
    static std::mutex mtx;
    utils::env::instance().set({
        {"GDALCUBES_STREAMING", "1"},
        {"GDALCUBES_STREAMING_CHUNK_ID", std::to_string(id)},
        {"GDALCUBES_STREAMING_FILE_IN", f_in},
        {"GDALCUBES_STREAMING_FILE_OUT", f_out}});
    mtx.lock();

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
        throw std::string("ERROR in stream_apply_pixel_cube::read_chunk(): external program returned exit code " +
                          std::to_string(exit_status));
    }
    GCBS_DEBUG(errstr);

    // read output data
    std::ifstream f_out_stream(f_out, std::ios::in | std::ios::binary);
    if (!f_out_stream.is_open()) {
        GCBS_ERROR("Cannot read streaming output data from file '" + f_out + "'");
        throw std::string(
            "ERROR in stream_apply_pixel_cube::read_chunk(): cannot read streaming output data from file '" +
            f_out + "'");
    }

    f_out_stream.seekg(0, f_out_stream.end);
    int length = f_out_stream.tellg();
    f_out_stream.seekg(0, f_out_stream.beg);
    char *buffer = (char *)std::calloc(length, sizeof(char));
    f_out_stream.read(buffer, length);
    f_out_stream.close();

    // Copy results to chunk buffer, at most the size of the output

    uint32_t offset = _keep_bands ? (inbuf->size()[0] * inbuf->size()[1] * inbuf->size()[2] * inbuf->size()[3]) : 0;
    if (_keep_bands) {
        std::memcpy(out->buf(), inbuf->buf(),
                    sizeof(double) * offset);
    }

    std::memcpy(((double *)(out->buf())) + offset, buffer + (4 * sizeof(int)), std::min(length - 4 * sizeof(int), _nbands * size_btyx[1] * size_btyx[2] * size_btyx[3] * sizeof(double)));
    std::free(buffer);

    if (filesystem::exists(f_out)) {
        filesystem::remove(f_out);
    }

    return out;
}

}  // namespace gdalcubes