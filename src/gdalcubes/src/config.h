/*
    MIT License

    Copyright (c) 2020 Marius Appel <marius.appel@hs-bochum.de>

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

#ifndef CONFIG_H
#define CONFIG_H

#include <memory>

#include "build_info.h"
#include "error.h"
#include "filesystem.h"
#include "progress.h"

namespace gdalcubes {

// forward declarations
class chunk_processor;
class chunk_processor_singlethread;

/**
 * A simple structure to store
 * version information of the library
 */
struct version_info {
    /**
     * @brief Major version number
     */
    uint16_t VERSION_MAJOR;

    /**
     * @brief Minor version number
     */
    uint16_t VERSION_MINOR;

    /**
     * @brief Patch version number
     */
    uint16_t VERSION_PATCH;

    /**
     * @brief Build data (from __DATE__ macro)
     */
    std::string BUILD_DATE;

    /**
    * @brief Build data (from __TIME__ macro)
    */
    std::string BUILD_TIME;

    /**
    * @brief git description (from git describe), i.e., tag and commit
    */
    std::string GIT_DESC;

    /**
    * @brief Last git commit (build might contain local changes anyway)
    */
    std::string GIT_COMMIT;
};

/**
 * @brief A singleton class to manage global configuration options
 */
class config {
   public:
    /**
     * Return the singleton instance
     */
    static config* instance() {
        static GC g;
        if (!_instance) {
            _instance = new config();
        }
        return _instance;
    }

    inline std::shared_ptr<chunk_processor> get_default_chunk_processor() {
        return _chunk_processor;
    }
    inline void set_default_chunk_processor(std::shared_ptr<chunk_processor> p) {
        _chunk_processor = p;
    }

    inline void set_default_progress_bar(std::shared_ptr<progress> p) {
        _progress_bar = p;
    }

    inline std::shared_ptr<progress> get_default_progress_bar() {
        return _progress_bar;
    }

    void set_gdal_cache_max(uint32_t size_bytes);

    inline void set_server_chunkcache_max(uint32_t size_bytes) {
        _server_chunkcache_max = size_bytes;
    }

    inline uint32_t get_server_chunkcache_max() {
        return _server_chunkcache_max;
    }

    inline void set_server_worker_threads_max(uint16_t max_threads) {
        _server_worker_threads_max = max_threads;
    }

    inline uint16_t get_server_worker_threads_max() {
        return _server_worker_threads_max;
    }

    inline bool get_gdal_use_overviews() { return _gdal_use_overviews; }
    inline void set_gdal_use_overviews(bool use_overviews) { _gdal_use_overviews = use_overviews; }

    inline bool get_swarm_curl_verbose() { return _swarm_curl_verbose; }
    inline void set_swarm_curl_verbose(bool verbose) { _swarm_curl_verbose = verbose; }

    // Get / set directory where to store chunk data for file-based streaming. This
    // should ideally fast storage like a ramdisk such as /dev/shm
    inline std::string get_streaming_dir() { return _streaming_dir; }
    inline void set_streaming_dir(std::string dir) { _streaming_dir = dir; }

    inline bool get_gdal_debug() { return _gdal_debug; }
    void set_gdal_debug(bool debug);

    void set_gdal_log(std::string logfile);
    void set_gdal_num_threads(uint16_t threads);

    inline uint16_t get_gdal_num_threads() { return _gdal_num_threads; }

    /**
     * @brief Global gdalcubes library initialization function
     */
    void gdalcubes_init();

    /**
   * @brief Global gdalcubes library cleanup function
   */
    void gdalcubes_cleanup();

    /**
   * @brief Global gdalcubes version information function
   */
    version_info get_version_info();

    inline void set_error_handler(error_action f) {
        _error_handler = f;
    }

    inline error_action get_error_handler() {
        return _error_handler;
    }

    inline std::vector<std::string> get_collection_format_preset_dirs() {
        return _collection_format_preset_dirs;
    }

    void add_collection_format_preset_dir(std::string dir);

    void set_gdal_option(std::string key, std::string value);

    std::string gdal_version_info();
    std::vector<std::string> gdal_formats();
    bool gdal_has_geos();

    static void gdal_err_handler_silent(CPLErr eErrClass, int err_no, const char *msg);
    static void gdal_err_handler_default(CPLErr eErrClass, int err_no, const char *msg);


   private:
    std::shared_ptr<chunk_processor> _chunk_processor;
    std::shared_ptr<progress> _progress_bar;
    error_action _error_handler;
    uint32_t _gdal_cache_max;
    uint32_t _server_chunkcache_max;
    uint16_t _server_worker_threads_max;  // number of threads for parallel chunk reads
    bool _swarm_curl_verbose;
    uint16_t _gdal_num_threads;
    bool _gdal_debug;
    bool _gdal_use_overviews;
    std::string _streaming_dir;
    std::vector<std::string> _collection_format_preset_dirs;

   private:
    config();
    ~config() {}
    config(const config&) = delete;
    static config* _instance;

    class GC {
       public:
        ~GC() {
            if (config::_instance) {
                delete config::_instance;
                config::_instance = nullptr;
            }
        }
    };
};

}  // namespace gdalcubes

#endif  //CONFIG_H
