#include "gdalcubes/src/cube.h"


namespace gdalcubes {
  class chunk_processor_multiprocess : public chunk_processor {
  public:
    
    chunk_processor_multiprocess() : _cmd(""), _interrupted(false), _nworker(1) {}
    
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
    
    void apply(std::shared_ptr<cube> c,
               std::function<void(chunkid_t, std::shared_ptr<chunk_data>, std::mutex &)> f) override;
    
    static void exec(std::string json_path, uint16_t pid, uint16_t nworker, std::string work_dir);
    void kill_all() {_interrupted = true;}
    
    
    
  private:
    std::string _cmd;
    bool _interrupted;
    uint16_t _nworker;
  
    
  };
}