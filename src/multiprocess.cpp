
#include "multiprocess.h"
#include "gdalcubes/src/cube_factory.h"
#include "gdalcubes/external/tiny-process-library/process.hpp"

#include <fstream>

namespace gdalcubes {

void chunk_processor_multiprocess::apply(std::shared_ptr<cube> c,
                                         std::function<void(chunkid_t, std::shared_ptr<chunk_data>, std::mutex &)> f) {
  
  GCBS_DEBUG("Using " + std::to_string(this->_nworker) + " worker");
  std::string cube_id = utils::generate_unique_filename();
  std::mutex mutex;
  std::string work_dir = filesystem::join(config::instance()->get_streaming_dir(), cube_id);
  if (filesystem::exists(work_dir)) {
    GCBS_ERROR("Directory '" + work_dir + "' for storing intermediate chunk data already exists");
    throw std::string("Directory '" + work_dir + "' for storing intermediate chunk data already exists");
  }
  filesystem::mkdir(work_dir);

  std::string json_path =filesystem::join(work_dir,"cube.json");
  std::ofstream jsonfile(json_path);
  jsonfile << c->make_constructible_json().dump();
  jsonfile.close();  
  
  uint16_t nworker = _nworker;
  std::vector<std::shared_ptr<TinyProcessLib::Process>> p;
  std::vector<std::string> err;
  err.resize(nworker,"");
  std::vector<bool> finished(nworker, false);
  bool all_finished = false;
  std::vector<int> exit_status(nworker, -1);
  for (uint16_t pid=0; pid < nworker; ++pid) {
    
    
#ifdef _WIN32
    _putenv("GDALCUBES_WORKER=1");
    _putenv((std::string("GDALCUBES_WORKER_JSON") + "=" + json_path.c_str()));
    _putenv((std::string("GDALCUBES_WORKER_ID") + "=" + std::to_string(pid).c_str()));
    _putenv((std::string("GDALCUBES_WORKER_N") + "=" + std::to_string(nworker).c_str()).c_str());
    _putenv((std::string("GDALCUBES_WORKER_DIR") + "=" + work_dir.c_str()).c_str());
#else
    setenv("GDALCUBES_WORKER", "1", 1);
    setenv("GDALCUBES_WORKER_JSON", json_path.c_str(), 1);
    setenv("GDALCUBES_WORKER_ID", std::to_string(pid).c_str(), 1);
    setenv("GDALCUBES_WORKER_N", std::to_string(nworker).c_str(), 1);
    setenv("GDALCUBES_WORKER_DIR", work_dir.c_str(), 1);
#endif
    
    // TODO: unset environment variables after process as been started
    auto pp = std::make_shared<TinyProcessLib::Process>(
      _cmd, "", [](const char *bytes, std::size_t n) {},
      [&err, &pid](const char *bytes, std::size_t n) {
        //err[pid] = std::string(bytes, n);
        //GCBS_DEBUG(err[pid]);
      },
      false);
    p.push_back(pp);
  }
  
  auto start = std::chrono::system_clock::now();
  while (!all_finished && !_interrupted) {
    
    
    // 1. check if processes are still running
    all_finished = true;
    for (uint16_t pid=0; pid < nworker; ++pid) {
      if(p[pid]->try_get_exit_status(exit_status[pid])) {
        finished[pid] = true;
      }
      else {
        all_finished = false;
      }
    }
    
    // 2. check whether there are finished chunks
    // Notice that this must be executed even if all_finished is true,
    // because there may be remaining chunks
    std::vector<std::pair<std::string, chunkid_t>> chunk_queue;
    filesystem::iterate_directory(work_dir, [&chunk_queue](const std::string& f) {
     
      // Consider files with name X.nc, where X is an integer number
      // Temporary files will start with a dot and are NOT considered here
      std::regex r("^([0-9]+)\\.nc$");
      std::smatch match;
      std::string basename  = filesystem::stem(f) + "." + filesystem::extension(f);
      if (std::regex_search(basename, match, r)) {
        if (match.size() == 2) {
          try {
            int chunkid = std::stoi(match[1].str());
            chunk_queue.push_back(std::make_pair<>(f, chunkid));
          }
          catch (...) {
            return;
          }
        }
      }
    });
    
    // add finished chunks to output
    for (auto it = chunk_queue.begin(); it != chunk_queue.end(); ++it) {
      try {
        GCBS_DEBUG("Merging chunk " + std::to_string(it->second) + " from " + it->first);
        std::shared_ptr<chunk_data> dat = std::make_shared<chunk_data>();
        dat->read_ncdf(it->first);
        f(it->second, dat, mutex);
        filesystem::remove(it->first);
        
        // for debugging only
        //filesystem::move(it->first,it->first + "DONE.nc");
        
        
      } catch (std::string s) {
        GCBS_ERROR(s);
        continue;
      } catch (...) {
        GCBS_ERROR("unexpected exception while processing chunk");
        continue;
      }
    }
    
    
    // 3. sleep and check user interrupt
    if (!all_finished) {
      if (_interrupted) {
        for (uint16_t pid=0; pid < nworker; ++pid) {
          if(p[pid]->try_get_exit_status(exit_status[pid])) {
            p[pid]->kill();
          }
        }
      }
      else {
        auto end = std::chrono::system_clock::now();
        auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
        if (elapsed > std::chrono::milliseconds(2000)) {
          // TODO: check interrupt
          start = std::chrono::system_clock::now();
        }
        std::this_thread::sleep_for(std::chrono::milliseconds(200));
      }
    }
  }
  if (_interrupted) {
    for (uint16_t pid=0; pid < nworker; ++pid) {
      if(p[pid]->try_get_exit_status(exit_status[pid])) {
        p[pid]->kill();
      }
    }
    GCBS_ERROR("computations have been interrupted by the user");
    
    // TODO: clean up potential intermediate files and similar?!
  }
  else {
    // check error status
    for (uint16_t pid=0; pid < nworker; ++pid) {
      if(exit_status[pid] != 0) {
        GCBS_ERROR("worker process " + std::to_string(pid) + " returned error code " + std::to_string(exit_status[pid]));
      }
    }
  }
  filesystem::remove(work_dir);
  _interrupted = false;
}

void chunk_processor_multiprocess::exec(std::string json_path, uint16_t pid, uint16_t nworker, std::string work_dir) {
  std::shared_ptr<cube> cube = cube_factory::instance()->create_from_json_file(json_path);
  
  // TODO: how to set important config options (e.g., netCDF compression, GDAL cache size, error handlers, ...)
  for (uint32_t i=pid; i<cube->count_chunks(); i+= nworker) {
    chunkid_t id = i;
    std::string outfile =  filesystem::join(work_dir, std::to_string(id) + ".nc");
    std::string outfile_temp =  filesystem::join(work_dir, "." + std::to_string(id) + ".nc");
    
    cube->read_chunk(id)->write_ncdf(outfile_temp); // TODO: add compression level and force?!
    filesystem::move(outfile_temp, outfile);
  }
  
}

}