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

#include "config.h"
#include "cube.h"

config* config::_instance = nullptr;

config::config() : _chunk_processor(std::make_shared<chunk_processor_singlethread>()),
                   _progress_bar(std::make_shared<progress_none>()),
                   _error_handler(error_handler::default_error_handler),
                   _gdal_cache_max(1024 * 1024 * 256),         // 256 MiB
                   _server_chunkcache_max(1024 * 1024 * 512),  // 512 MiB
                   _server_worker_threads_max(1),
                   _swarm_curl_verbose(false),
                   _gdal_num_threads(1),
                   _collection_format_preset_dirs() {
}

version_info config::get_version_info() {
    version_info v;
    v.VERSION_MAJOR = GDALCUBES_VERSION_MAJOR;
    v.VERSION_MINOR = GDALCUBES_VERSION_MINOR;
    v.VERSION_PATCH = GDALCUBES_VERSION_PATCH;
    v.BUILD_DATE = __DATE__;
    v.BUILD_TIME = __TIME__;
    v.GIT_COMMIT = GDALCUBES_GIT_COMMIT;
    v.GIT_DESC = GDALCUBES_GIT_DESC;
    return v;
}