#include "reduce_time_stream.h"
#include "gdalcubes/src/external/tiny-process-library/process.hpp"


std::shared_ptr<chunk_data> reduce_time_stream_cube::read_chunk(chunkid_t id) {
  GCBS_TRACE("reduce_time_stream_cube::read_chunk(" + std::to_string(id) + ")");
  std::shared_ptr<chunk_data> out = std::make_shared<chunk_data>();
  if (id >= count_chunks())
    return out;  // chunk is outside of the view, we don't need to read anything.
  
  coords_nd<uint32_t, 3> size_tyx = chunk_size(id);
  coords_nd<uint32_t, 4> size_btyx = {_nbands, 1, size_tyx[1], size_tyx[2]};
  out->size(size_btyx);
  
  // Fill buffers accordingly
  out->buf(std::calloc(size_btyx[0] * size_btyx[1] * size_btyx[2] * size_btyx[3], sizeof(double)));
  double *begin = (double *)out->buf();
  double *end = ((double *)out->buf()) + size_btyx[0] * size_btyx[1] * size_btyx[2] * size_btyx[3];
  std::fill(begin, end, NAN);

  
  // 1. read everything to input buffer (for first version, can be memory-intensive, same as rechunk_merge_time)
  std::shared_ptr<chunk_data> inbuf = std::make_shared<chunk_data>();
  coords_nd<uint32_t, 3> in_size_tyx = chunk_size(id);
  coords_nd<uint32_t, 4> in_size_btyx = {uint32_t(_in_cube->size_bands()), _in_cube->size_t(), size_tyx[1], size_tyx[2]};
  inbuf->size(in_size_btyx);
  inbuf->buf(std::calloc(in_size_btyx[0] * in_size_btyx[1] * in_size_btyx[2] * in_size_btyx[3], sizeof(double)));
  double *inbegin = (double *)inbuf->buf();
  double *inend = ((double *)inbuf->buf()) + in_size_btyx[0] * in_size_btyx[1] * in_size_btyx[2] * in_size_btyx[3];
  std::fill(inbegin, inend, NAN);

  uint32_t ichunk=0;
  for (chunkid_t i = id; i < _in_cube->count_chunks(); i += _in_cube->count_chunks_x() * _in_cube->count_chunks_y()) {
    std::shared_ptr<chunk_data> x = _in_cube->read_chunk(i);
    for (uint16_t ib=0; ib < x->size()[0]; ++ib) {
      for (uint32_t it = 0; it < x->size()[1]; ++it) {
        for (uint32_t ixy = 0; ixy < x->size()[2] * x->size()[3]; ++ixy) {
          ((double*)inbuf->buf())[ib * _in_cube->size_t() *  x->size()[2] *  x->size()[3] + (it +  ichunk*_in_cube->chunk_size()[0]) * x->size()[2] *  x->size()[3] + ixy] = 
            ((double*)x->buf())[ib *  x->size()[1] *  x->size()[2] *  x->size()[3] + it * x->size()[2] *  x->size()[3] + ixy];
        }
      }
    }
    ++ichunk;
  }
  
  
  
  // generate in and out filename
  std::string f_in = filesystem::join(config::instance()->get_streaming_dir(), utils::generate_unique_filename(12, ".stream_", "_in"));
  std::string f_out = filesystem::join(config::instance()->get_streaming_dir(), utils::generate_unique_filename(12, ".stream_", "_out"));
  
  std::string errstr;  // capture error string
  
  // write input data
  std::ofstream f_in_stream(f_in, std::ios::out | std::ios::binary | std::ios::trunc);
  if (!f_in_stream.is_open()) {
    GCBS_ERROR("Cannot write streaming input data to file '" + f_in + "'");
    throw std::string("ERROR in stream_cube::stream_chunk_file(): cannot write streaming input data to file '" + f_in + "'");
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
  for (int i = 0; i < size[1]; ++i) {
    dims[i] = (_in_cube->st_reference()->t0() + _in_cube->st_reference()->dt() * i).to_double();
  }
  for (int i = size[1]; i < size[1] + size[2]; ++i) { // FIXME: coordinates are wrong, must start at chunk boundary
    dims[i] = _in_cube->st_reference()->win().bottom + i * _in_cube->st_reference()->dy();
  }
  for (int i = size[1] + size[2]; i < size[1] + size[2] + size[3]; ++i) { // FIXME: coordinates are wrong, must start at chunk boundary
    dims[i] = _in_cube->st_reference()->win().left + i * _in_cube->st_reference()->dx();
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
  mtx.lock();
#ifdef _WIN32
  _putenv("GDALCUBES_STREAMING=1");
  //_putenv((std::string("GDALCUBES_STREAMING_DIR") + "=" + config::instance()->get_streaming_dir().c_str()).c_str());
  _putenv((std::string("GDALCUBES_STREAMING_FILE_IN") + "=" + f_in.c_str()).c_str());
  _putenv((std::string("GDALCUBES_STREAMING_FILE_OUT") + "=" + f_out.c_str()).c_str());
#else
  setenv("GDALCUBES_STREAMING", "1", 1);
  // setenv("GDALCUBES_STREAMING_DIR", config::instance()->get_streaming_dir().c_str(), 1);
  setenv("GDALCUBES_STREAMING_FILE_IN", f_in.c_str(), 1);
  setenv("GDALCUBES_STREAMING_FILE_OUT", f_out.c_str(), 1);
#endif
  
  
  
  
  // start process
  TinyProcessLib::Process process(_cmd, "", [](const char *bytes, std::size_t n) {}, [&errstr](const char *bytes, std::size_t n) {
    errstr = std::string(bytes, n);
    GCBS_DEBUG(errstr); }, false);
  mtx.unlock();
  auto exit_status = process.get_exit_status();
  filesystem::remove(f_in);
  if (exit_status != 0) {
    GCBS_ERROR("Child process failed with exit code " + std::to_string(exit_status));
    GCBS_ERROR("Child process output: " + errstr);
    if (filesystem::exists(f_out)) {
      filesystem::remove(f_out);
    }
    throw std::string("ERROR in reduce_time_stream_cube::read_chunk(): external program returned exit code " + std::to_string(exit_status));
  }
  
  // read output data
  std::ifstream f_out_stream(f_out, std::ios::in | std::ios::binary);
  if (!f_out_stream.is_open()) {
    GCBS_ERROR("Cannot read streaming output data from file '" + f_out + "'");
    throw std::string("ERROR in reduce_time_stream_cube::read_chunk(): cannot read streaming output data from file '" + f_out + "'");
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