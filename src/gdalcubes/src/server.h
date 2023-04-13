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

#ifndef SERVER_H
#define SERVER_H

#include <cpprest/http_listener.h>
#include <cpprest/uri_builder.h>

#include <list>
#include <mutex>
#include <queue>
#include <thread>

#include "cube.h"

namespace gdalcubes {

/**
 * @brief An in-memory singleton cache for successfully read / computed chunks
 *
 * Chunks are identified by std::pair<cube_id, chunk_id>
 */
class server_chunk_cache {
   public:
    /**
     * @brief Get the singleton instance
     * @return pointer to the singleton instance
     */
    static server_chunk_cache* instance() {
        static GC g;
        _singleton_mutex.lock();
        if (!_instance) {
            _instance = new server_chunk_cache();
        }
        _singleton_mutex.unlock();
        return _instance;
    }

    /**
     * @brief Remove a chunk from the cache
     * @param key
     */
    void remove(std::pair<uint32_t, uint32_t> key) {
        _m.lock();
        auto it = _cache.find(key);
        if (it != _cache.end()) {
            _size_bytes -= it->second->total_size_bytes();
            _cache.erase(it);
        }

        auto it_2 = _prio_forward.find(key);
        uint32_t p = it_2->second;
        if (it_2 != _prio_forward.end()) {
            _prio_forward.erase(it_2);
        }

        auto it_3 = _prio_backward.find(p);
        if (it_3 != _prio_backward.end()) {
            _prio_backward.erase(it_3);
        }
        _m.unlock();
    }

    /**
     * Add chunk data to the cache
     * @param key chunk identifier (cube_id, chunk_id)
     * @param value chunk data to add
     */
    void add(std::pair<uint32_t, uint32_t> key, std::shared_ptr<chunk_data> value) {
        if (!has(key)) {
            _m.lock();
            while (_size_bytes + value->total_size_bytes() > config::instance()->get_server_chunkcache_max()) {
                auto it = _prio_backward.lower_bound(0);  // lowest value greater than or equal to 0

                // remove it from the data
                remove(it->second);
            }

            _cache[key] = value;
            _size_bytes += value->total_size_bytes();

            uint64_t p = inc_prio();  // synchronize?

            _prio_forward[key] = p;
            _prio_backward[p] = key;

            _size_bytes += value->total_size_bytes();
            _m.unlock();
        }
    }

    /**
     * Check whether a chunk is cached
     * @param key chunk key (cube_id, chunk_id)
     * @return true, if the chunk is cached
     */
    inline bool has(std::pair<uint32_t, uint32_t> key) {
        return _cache.find(key) != _cache.end();
    }

    /**
     * Get chunk data from the cache
     * @param key chunk key (cube_id, chunk_id)
     * @return chunk data as shared_ptr
     */
    std::shared_ptr<chunk_data> get(std::pair<uint32_t, uint32_t> key) {
        if (!has(key)) {
            throw std::string("ERROR: in server_chunk_cache::get(): requested chunk is not available");
        }
        _m.lock();
        uint64_t pnew = inc_prio();
        auto it_1 = _prio_forward.find(key);
        if (it_1 != _prio_forward.end()) {
            uint32_t pold = it_1->second;
            it_1->second = pnew;

            auto it_2 = _prio_backward.find(pold);
            if (it_2 != _prio_backward.end()) {
                _prio_backward.erase(it_2);  // remove old
                _prio_backward[pnew] = key;  // add new
            }
        }
        _m.unlock();
        return _cache[key];
    }

    /**
     * @brief Get the total amount of memory currently consumed by the cache
     * @return Size of the cache in bytes
     */
    inline uint32_t total_size_bytes() {
        return _size_bytes;
    }

   private:
    const uint64_t inc_prio() {
        if (_cur_priority == std::numeric_limits<uint64_t>::max()) {
            uint64_t newp = 0;
            for (auto it = _prio_forward.begin(); it != _prio_forward.end(); ++it) {
                it->second = newp++;
                _prio_backward.insert(std::make_pair(it->second, it->first));
            }
            _cur_priority = newp;
        } else {
            ++_cur_priority;
        }
        return _cur_priority;
    }

    std::map<std::pair<uint32_t, uint32_t>, std::shared_ptr<chunk_data>> _cache;
    std::map<uint64_t, std::pair<uint32_t, uint32_t>> _prio_backward;
    std::map<std::pair<uint32_t, uint32_t>, uint64_t> _prio_forward;

    std::mutex _m;
    uint64_t _size_bytes;

    static std::mutex _singleton_mutex;
    uint64_t _cur_priority;

   private:
    server_chunk_cache() : _cache(), _prio_backward(), _prio_forward(), _m(), _size_bytes(0), _cur_priority(0) {}
    ~server_chunk_cache() {}
    server_chunk_cache(const server_chunk_cache&) = delete;
    static server_chunk_cache* _instance;

    class GC {
       public:
        ~GC() {
            if (server_chunk_cache::_instance) {
                delete server_chunk_cache::_instance;
                server_chunk_cache::_instance = nullptr;
            }
        }
    };
};

/**
 * @brief Serve gdalcubes functionality as a REST-like API over HTTP on a provided host, port, and endpoint
 *
 * This class enables distributed processing by providing gdalcubes functionality of a simple HTTP REST-like API
 */
class gdalcubes_server {
   public:
    gdalcubes_server(std::string host, uint16_t port = 1111, std::string basepath = "gdalcubes/api/", bool ssl = false, std::string workdir = filesystem::join(filesystem::get_tempdir(), "gdalcubes_server"), std::set<std::string> whitelist = {}) : _listener(),
                                                                                                                                                                                                                                                       _port(port),
                                                                                                                                                                                                                                                       _host(host),
                                                                                                                                                                                                                                                       _basepath(basepath),
                                                                                                                                                                                                                                                       _ssl(ssl),
                                                                                                                                                                                                                                                       _workdir(workdir),
                                                                                                                                                                                                                                                       _cubestore(),
                                                                                                                                                                                                                                                       _cur_id(0),
                                                                                                                                                                                                                                                       _chunk_cond(),
                                                                                                                                                                                                                                                       _whitelist(whitelist) {
        if (filesystem::exists(_workdir) && filesystem::is_directory(_workdir)) {
            // boost::filesystem::remove_all(_workdir); // TODO: uncomment after testing
        } else if (filesystem::exists(_workdir) && !filesystem::is_directory(_workdir)) {
            throw std::string("ERROR in gdalcubes_server::gdalcubes_server(): working directory for gdalcubes_server is an existing file.");
        }

        filesystem::mkdir(_workdir);

        std::string url = ssl ? "https://" : "http://";
        url += host + ":" + std::to_string(port);
        web::uri_builder uri(url);
        uri.append_path(basepath);

        _listener = web::http::experimental::listener::http_listener(uri.to_string());

        _listener.support(web::http::methods::GET, std::bind(&gdalcubes_server::handle_get, this, std::placeholders::_1));
        _listener.support(web::http::methods::POST, std::bind(&gdalcubes_server::handle_post, this, std::placeholders::_1));
        _listener.support(web::http::methods::HEAD, std::bind(&gdalcubes_server::handle_head, this, std::placeholders::_1));
    }

   public:
    inline pplx::task<void> open() { return _listener.open(); }
    inline pplx::task<void> close() { return _listener.close(); }

    inline std::string get_service_url() { return _listener.uri().to_string(); }

   private:
    void handle_get(web::http::http_request req);
    void handle_post(web::http::http_request req);
    void handle_head(web::http::http_request req);

    web::http::experimental::listener::http_listener _listener;

    inline uint32_t get_unique_id() {
        _mutex_id.lock();
        uint32_t id = _cur_id++;
        _mutex_id.unlock();
        return id;
    }

    const uint16_t _port;
    const std::string _host;
    const std::string _basepath;
    const bool _ssl;
    const std::string _workdir;

    std::map<uint32_t, std::shared_ptr<cube>> _cubestore;

    uint16_t _cur_id;
    std::mutex _mutex_id;
    std::mutex _mutex_cubestore;

    uint16_t _worker_thread_count;

    std::mutex _mutex_chunk_read_requests;
    std::list<std::pair<uint32_t, uint32_t>> _chunk_read_requests;
    std::set<std::pair<uint32_t, uint32_t>> _chunk_read_requests_set;

    std::mutex _mutex_chunk_read_executing;
    std::set<std::pair<uint32_t, uint32_t>> _chunk_read_executing;

    std::vector<std::thread> _worker_threads;
    std::mutex _mutex_worker_threads;

    std::condition_variable _worker_cond;
    std::mutex _mutex_worker_cond;

    std::map<std::pair<uint32_t, uint32_t>, std::pair<std::condition_variable, std::mutex>> _chunk_cond;

    std::set<std::string> _whitelist;
};

}  // namespace gdalcubes

#endif  //SERVER_H
