
#include "multiprocess.h"
#include "gdalcubes/src/cube_factory.h"
#include "gdalcubes/src/external/tiny-process-library/process.hpp"
#include "error.h"

// [[Rcpp::plugins("cpp11")]]
#include <Rcpp.h>
#include <fstream>

namespace gdalcubes {

void chunk_processor_multiprocess::apply(std::shared_ptr<cube> c,
                                         std::function<void(chunkid_t, std::shared_ptr<chunk_data>, std::mutex &)> f) {
  
  GCBS_DEBUG("Using " + std::to_string(this->_nworker) + " worker");
  _interrupted = false;
  
  std::string job_id = utils::generate_unique_filename();
  std::mutex mutex;
  std::string work_dir = filesystem::join(config::instance()->get_streaming_dir(), job_id);
  if (filesystem::exists(work_dir)) {
    GCBS_ERROR("Directory '" + work_dir + "' for storing intermediate chunk data already exists");
    throw std::string("Directory '" + work_dir + "' for storing intermediate chunk data already exists");
  }
  filesystem::mkdir(work_dir);
  GCBS_DEBUG("Using '" + work_dir + "' as working directory for child processes");

  std::string json_path =filesystem::join(work_dir,"cube.json");
  std::ofstream jsonfile(json_path);
  jsonfile << c->make_constructible_json().dump();
  jsonfile.close();  
  
  uint16_t nworker = _nworker;
  
  std::vector<std::shared_ptr<TinyProcessLib::Process>> p;
  std::vector<bool> finished(nworker, false);
  bool all_finished = false;
  bool any_failed = false;
  std::vector<int> exit_status(nworker, -1);
  
  
  json11::Json::object j_gdal_options;
  for (auto it = _gdal_options.begin(); it != _gdal_options.end(); ++it) {
    j_gdal_options[it->first.c_str()] = it->second;
  }
  
  for (uint16_t pid=0; pid < nworker; ++pid) {
    json11::Json j = json11::Json::object{ 
      {"job_id", job_id},
      {"worker_id", pid},
      {"worker_count", nworker},
      {"job_start", ""}, // TODO
      {"workdir", work_dir},
      {"cube", filesystem::join(work_dir, "cube.json")},
      {"gdalcubes_options", json11::Json::object{
        {"debug", _debug}, 
        {"log_file", filesystem::join(work_dir, "worker_" + std::to_string(pid) + ".log")},
        {"ncdf_compression_level", _ncdf_compression_level}, 
        {"streaming_dir", work_dir},
        {"use_overview_images", _use_overviews}
      }},
      {"gdal_options",j_gdal_options}
    }; 
    
    // write json to file
    std::ofstream ojson;
    std::string worker_json_file = filesystem::join(work_dir, "worker_" + std::to_string(pid) + ".json");
    ojson.open(worker_json_file);
    ojson << j.dump();
    ojson.close();
    
    // start child process, with first argument being path to the json worker process description
    TinyProcessLib::Config pconf;
    pconf.show_window = TinyProcessLib::Config::ShowWindow::hide;
    auto pp = std::make_shared<TinyProcessLib::Process>(
      _cmd + " " + worker_json_file, "", nullptr,
      nullptr, false, pconf);
    p.push_back(pp);
    
  }
  
  auto start = std::chrono::system_clock::now();
  while (!all_finished) {
  
    // 1. check if processes are still running
    all_finished = true;
    any_failed = false;
    for (uint16_t pid=0; pid < nworker; ++pid) {
      if(p[pid]->try_get_exit_status(exit_status[pid])) {
        finished[pid] = true;
        if (exit_status[pid] > 0) {
          any_failed = true;
        }
      }
      else {
        all_finished = false;
      }
    }
    if (any_failed) {
      break;
    }
    
    
    // 2. check whether there are finished chunks
    // Notice that this must be executed even if all_finished is true,
    // because there may be remaining chunks
    std::vector<std::pair<std::string, chunkid_t>> chunk_queue;
    filesystem::iterate_directory(work_dir, [&chunk_queue](const std::string& f) {
      
      // Consider files with name X.nc, where X is an integer number
      // Temporary files will start with a dot and are NOT considered here
      std::string basename  = filesystem::stem(f) + "." + filesystem::extension(f);
      std::size_t pos = basename.find(".nc");
      if (pos > 0 && pos < std::string::npos) {
        try {
            int chunkid = std::stoi(basename.substr(0,pos));
            chunk_queue.push_back(std::make_pair<>(f, chunkid));
          }
          catch (...) {}
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
        Rcpp::warning("Chunk" + std::to_string(it->second) + " could not be added to output.");
        continue;
      } catch (...) {
        GCBS_ERROR("unexpected exception while processing chunk");
        Rcpp::warning("Chunk" + std::to_string(it->second) + " could not be added to output.");
        continue;
      }
    }
    
    
    // 3. sleep and check user interrupt
    if (!all_finished) {
      auto end = std::chrono::system_clock::now();
      auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
      if (elapsed > std::chrono::milliseconds(1500)) {
       try {
         Rcpp::checkUserInterrupt();
       }
       catch (...) {
         _interrupted = true;
         break;
       }
        start = std::chrono::system_clock::now();
      }
      std::this_thread::sleep_for(std::chrono::milliseconds(200));
    }
  }
  
  if (_interrupted) {
    for (uint16_t pid=0; pid < nworker; ++pid) {
      if(!p[pid]->try_get_exit_status(exit_status[pid])) {
        GCBS_DEBUG("killing worker process #" + std::to_string(pid) );
        p[pid]->kill(true);
      }
    }
    GCBS_ERROR("computations have been interrupted by the user");
    
    // TODO: clean up potential intermediate files and similar?!
  }
  else if (any_failed) {
    for (uint16_t pid=0; pid < nworker; ++pid) {
      if(!p[pid]->try_get_exit_status(exit_status[pid])) {
        GCBS_DEBUG("killing worker process #" + std::to_string(pid) + " because one or more other worker processes failed" );
        p[pid]->kill(true);
      }
      else {
        if (exit_status[pid] > 0) {
          GCBS_ERROR("worker process #" + std::to_string(pid) + " returned " + std::to_string(exit_status[pid]));
        }
      }
    }
    
  }
  else { // all_finished is true
    // check error status anyway
    for (uint16_t pid=0; pid < nworker; ++pid) {
      if(exit_status[pid] > 0) {
        GCBS_ERROR("worker process #" + std::to_string(pid) + " returned " + std::to_string(exit_status[pid]));
        any_failed = true;
      }
    }
  }
  
  // Print output from worker processes
  if (_debug || any_failed) {
    for (uint16_t pid=0; pid < nworker; ++pid) {
      std::ifstream wlogf(filesystem::join(work_dir, "worker_" + std::to_string(pid) + ".log")); // TODO: take path from JSON
      while(!wlogf.eof()) 
      {
        std::string msg;
        std::getline(wlogf, msg);
        if (!msg.empty()) {
          std::stringstream sss;
          sss << "[WORKER #" << std::to_string(pid) << "] " <<  msg << std::endl;
          r_stderr_buf::print(sss.str());
        }
      }
    }
  }
  
  filesystem::remove(work_dir);
  r_stderr_buf::print(); // make sure that deferred output is printed
  
  // Error message
  if(_interrupted) {
    _interrupted = false;
    Rcpp::stop("computations have been interrupted by the user");
  }
  else if (any_failed) {
    Rcpp::stop("one or more worker processes failed to compute data cube chunks");
  }
}

void chunk_processor_multiprocess::exec(std::string json_path, uint16_t pid, uint16_t nworker, std::string work_dir, int ncdf_compression_level) {
  std::shared_ptr<cube> cube = cube_factory::instance()->create_from_json_file(json_path);
  
  for (uint32_t i=pid; i<cube->count_chunks(); i+= nworker) {
    chunkid_t id = i;
    std::string outfile =  filesystem::join(work_dir, std::to_string(id) + ".nc");
    std::string outfile_temp =  filesystem::join(work_dir, "." + std::to_string(id) + ".nc");
    
    // TODO: exception handling?!
    cube->read_chunk(id)->write_ncdf(outfile_temp, ncdf_compression_level); 
    if (filesystem::exists(outfile_temp)) {
      filesystem::move(outfile_temp, outfile);
    }
    // TODO: error handling / exceptions
  }
  
}

}