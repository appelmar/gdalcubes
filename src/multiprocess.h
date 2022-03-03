#ifndef GC_MULTIPROCESS_H
#define GC_MULTIPROCESS_H


#include "gdalcubes/src/cube.h"
#include <unordered_map>

namespace gdalcubes {
  class chunk_processor_multiprocess : public chunk_processor {
  public:
    
    chunk_processor_multiprocess() : _cmd(""), _interrupted(false), _nworker(1), _debug(false), _use_overviews(true), 
                                     _ncdf_compression_level(0), _gdal_options() {}
    
    uint32_t max_threads() override {
      return 1;
    }
    uint32_t nworker()  {
      return _nworker;
    }
    
    void set_cmd(std::string cmd) {
      _cmd = cmd;
    }
    
    void set_nworker(uint16_t n) {
      _nworker = n;
    }
    
    void set_debug(bool debug) {
      _debug = debug;
    }
    
    void set_ncdf_compression_level(int l) {
      _ncdf_compression_level = l;
    }
    
    void set_use_overviews(bool v) {
      _use_overviews = v;
    }
    
    void set_gdal_options(std::unordered_map<std::string, std::string> gdal_options) {
      _gdal_options = gdal_options;
    }
    
    void apply(std::shared_ptr<cube> c,
               std::function<void(chunkid_t, std::shared_ptr<chunk_data>, std::mutex &)> f) override;
    
    static void exec(std::string json_path, uint16_t pid, uint16_t nworker, std::string work_dir, int ncdf_compression_level = 0);
    void kill_all() {_interrupted = true;}
    
    
    
  private:
    std::string _cmd;
    bool _interrupted;
    uint16_t _nworker;
    bool _debug;
    bool _use_overviews;
    int _ncdf_compression_level;
    std::unordered_map<std::string, std::string> _gdal_options;
  };
}

#endif