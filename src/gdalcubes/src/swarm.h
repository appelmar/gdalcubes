/*
    MIT License

    Copyright (c) 2019 Marius Appel <marius.appel@uni-muenster.de>

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

#ifndef GDALCUBES_NO_SWARM

#ifndef SWARM_H
#define SWARM_H

#include <curl/curl.h>

#include "cube.h"

namespace gdalcubes {

/**
 * @brief Chunk processor implementation for distributed processing by connecting to gdalcubes_server instances
 *
 * This class connects the several gdalcubes_server instances in order to distribute read_chunk() operations
 *
 * @todo implement add / remove method for workers
 */
class gdalcubes_swarm : public chunk_processor {
   public:
    gdalcubes_swarm(std::vector<std::string> urls) : _cube(nullptr), _cube_ids(), _server_handles(), _server_uris(urls), _nthreads(1) {
        for (uint16_t i = 0; i < _server_uris.size(); ++i)
            _server_handles.push_back(curl_easy_init());
    }

    ~gdalcubes_swarm() {
        for (uint16_t i = 0; i < _server_uris.size(); ++i)
            curl_easy_cleanup(_server_handles[i]);
    }

    inline static std::shared_ptr<gdalcubes_swarm> from_urls(std::vector<std::string> urls) { return std::make_shared<gdalcubes_swarm>(urls); }

    // Create from txt file where each line is a uri
    static std::shared_ptr<gdalcubes_swarm> from_txtfile(std::string path);

    // upload execution context to all servers
    void push_execution_context(bool recursive = false);

    // create cube on all servers
    void push_cube(std::shared_ptr<cube> c);

    // Mimic cube::apply with distributed calls to cube::read_chunk()
    void apply(std::shared_ptr<cube> c, std::function<void(chunkid_t, std::shared_ptr<chunk_data>, std::mutex &)> f) override;

    /**
    * @copydoc chunk_processor::max_threads
    */
    uint32_t max_threads() override {
        return _nthreads;
    }
    inline void set_threads(uint16_t threads) { _nthreads = threads; }

   private:
    void post_file(std::string path, uint16_t server_index);

    void post_start(uint32_t chunk_id, uint16_t server_index);
    std::shared_ptr<chunk_data> get_download(uint32_t chunk_id, uint16_t server_index);

    uint32_t post_cube(std::string json, uint16_t server_index);

    std::shared_ptr<cube> _cube;
    std::vector<uint32_t> _cube_ids;  // IDs of the cube for each server

    std::vector<CURL *> _server_handles;
    std::vector<std::string> _server_uris;

    uint16_t _nthreads;
};

}  // namespace gdalcubes

#endif  //SWARM_H

#endif  // GDALCUBES_NO_SWARM