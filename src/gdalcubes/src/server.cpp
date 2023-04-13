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

#include "server.h"

#include <cpprest/filestream.h>
#include <cpprest/rawptrstream.h>

#include <boost/program_options.hpp>
#include <condition_variable>

#include "build_info.h"
#include "cube_factory.h"
#include "image_collection.h"
#include "utils.h"
/**
GET  /version
POST /file (name query, body file)
POST /cube (json process descr), return cube_id
GET /cube/{cube_id}
POST /cube/{cube_id}/{chunk_id}/start
GET /cube/{cube_id}/{chunk_id}/status status= "notrequested", "submitted" "queued" "running" "canceled" "finished" "error"
GET /cube/{cube_id}/{chunk_id}/download


 TODO:
 - implement SSL
 - implement auth
 - implement whitelist to accept connections from selected clients only
**/

namespace gdalcubes {

server_chunk_cache* server_chunk_cache::_instance = nullptr;
std::mutex server_chunk_cache::_singleton_mutex;

void gdalcubes_server::handle_get(web::http::http_request req) {
    if (!_whitelist.empty()) {
        std::string remote = req.remote_address();
        if (_whitelist.find(remote) == _whitelist.end()) {
            GCBS_DEBUG("Incoming request from " + req.remote_address() + " has been blocked according to whitelist rule");
            req.reply(web::http::status_codes::NotFound);
        }
        GCBS_DEBUG("Incoming request from " + req.remote_address() + " has been accepted according to whitelist rule");
    } else {
        GCBS_DEBUG("Incoming request from " + req.remote_address());
    }

    std::vector<std::string> path = web::uri::split_path(web::uri::decode(req.relative_uri().path()));
    std::map<std::string, std::string> query_pars = web::uri::split_query(web::uri::decode(req.relative_uri().query()));
    //    std::for_each(path.begin(), path.end(), [](std::string s) { std::cout << s << std::endl; });
    //    std::for_each(query_pars.begin(), query_pars.end(), [](std::pair<std::string, std::string> s) { std::cout << s.first << ":" << s.second << std::endl; });
    if (!path.empty()) {
        if (path[0] == "version") {
            GCBS_DEBUG("GET /version");
            version_info v = config::instance()->get_version_info();
            std::stringstream ss;
            ss << "gdalcubes_server " << v.VERSION_MAJOR << "." << v.VERSION_MINOR << "." << v.VERSION_PATCH << " (" << v.GIT_COMMIT << ") built on " << v.BUILD_DATE << " " << v.BUILD_TIME;
            req.reply(web::http::status_codes::OK, ss.str().c_str(), "text/plain");
        } else if (path[0] == "cube") {
            if (path.size() == 2) {
                uint32_t cube_id = std::stoi(path[1]);
                GCBS_DEBUG("GET /cube" + std::to_string(cube_id));

                if (_cubestore.find(cube_id) == _cubestore.end()) {
                    req.reply(web::http::status_codes::NotFound, "ERROR in /GET /cube/{cube_id}: cube with given id is not available.", "text/plain");
                } else {
                    req.reply(web::http::status_codes::OK, _cubestore[cube_id]->make_constructible_json().dump().c_str(),
                              "application/json");
                }
            } else if (path.size() == 4) {
                uint32_t cube_id = std::stoi(path[1]);
                uint32_t chunk_id = std::stoi(path[2]);

                std::string cmd = path[3];
                if (cmd == "download") {
                    GCBS_DEBUG("GET /cube/" + std::to_string(cube_id) + "/" + std::to_string(chunk_id) + "/download");
                    if (_cubestore.find(cube_id) == _cubestore.end()) {
                        req.reply(web::http::status_codes::NotFound, "ERROR in /GET /cube/{cube_id}/{chunk_id}/download: cube is not available", "text/plain");
                    } else if (chunk_id >= _cubestore[cube_id]->count_chunks()) {
                        req.reply(web::http::status_codes::NotFound, "ERROR in /GET /cube/{cube_id}/{chunk_id}/download: invalid chunk_id given", "text/plain");
                    }
                    // if not in queue, executing, or finished, return 404
                    else if (!server_chunk_cache::instance()->has(std::make_pair(cube_id, chunk_id)) &&
                             _chunk_read_requests_set.find(std::make_pair(cube_id, chunk_id)) == _chunk_read_requests_set.end() &&
                             _chunk_read_executing.find(std::make_pair(cube_id, chunk_id)) == _chunk_read_executing.end()) {
                        req.reply(web::http::status_codes::NotFound, "ERROR in /GET /cube/{cube_id}/{chunk_id}/download: chunk read has not been requested yet", "text/plain");
                    } else {
                        while (!server_chunk_cache::instance()->has(std::make_pair(cube_id, chunk_id))) {
                            // req.reply(web::http::status_codes::NotFound, "ERROR in /GET /cube/{cube_id}/{chunk_id}/download: chunk is not available.", "text/plain");
                            std::unique_lock<std::mutex> lck(_chunk_cond[std::make_pair(cube_id, chunk_id)].second);
                            _chunk_cond[std::make_pair(cube_id, chunk_id)].first.wait(lck);
                        }
                        std::shared_ptr<chunk_data> dat = server_chunk_cache::instance()->get(std::make_pair(cube_id, chunk_id));

                        uint8_t* rawdata = (uint8_t*)std::malloc(4 * sizeof(uint32_t) + dat->total_size_bytes());
                        memcpy((void*)rawdata, (void*)(dat->size().data()), 4 * sizeof(uint32_t));
                        if (!dat->empty()) {
                            memcpy(rawdata + 4 * sizeof(uint32_t), dat->buf(), dat->total_size_bytes());
                        }

                        concurrency::streams::basic_istream<uint8_t> is = concurrency::streams::rawptr_stream<uint8_t>::open_istream(rawdata, 4 * sizeof(uint32_t) + dat->total_size_bytes());
                        req.reply(web::http::status_codes::OK, is, 4 * sizeof(uint32_t) + dat->total_size_bytes(), "application/octet-stream").then([&is, &rawdata]() {
                            is.close();
                            std::free(rawdata); });
                    }

                } else if (cmd == "status") {
                    GCBS_DEBUG("GET /cube/" + std::to_string(cube_id) + "/" + std::to_string(chunk_id) + "/status");
                    if (_cubestore.find(cube_id) == _cubestore.end()) {
                        req.reply(web::http::status_codes::NotFound, "ERROR in /GET /cube/{cube_id}/{chunk_id}/status: cube is not available", "text/plain");
                    } else if (chunk_id >= _cubestore[cube_id]->count_chunks()) {
                        req.reply(web::http::status_codes::NotFound, "ERROR in /GET /cube/{cube_id}/{chunk_id}/status: invalid chunk_id given", "text/plain");
                    }

                    if (server_chunk_cache::instance()->has(std::make_pair(cube_id, chunk_id))) {
                        req.reply(web::http::status_codes::OK, "finished", "text/plain");
                    } else if (_chunk_read_executing.find(std::make_pair(cube_id, chunk_id)) !=
                               _chunk_read_executing.end()) {
                        req.reply(web::http::status_codes::OK, "running", "text/plain");
                    } else if (_chunk_read_requests_set.find(std::make_pair(cube_id, chunk_id)) != _chunk_read_requests_set.end()) {
                        req.reply(web::http::status_codes::OK, "queued", "text/plain");
                    } else {
                        req.reply(web::http::status_codes::OK, "notrequested", "text/plain");
                    }

                } else {
                    req.reply(web::http::status_codes::NotFound);
                }
            } else {
                req.reply(web::http::status_codes::NotFound);
            }
        }

        else {
            req.reply(web::http::status_codes::NotFound);
        }
    } else {
        req.reply(web::http::status_codes::NotFound);
    }
}

void gdalcubes_server::handle_post(web::http::http_request req) {
    if (!_whitelist.empty()) {
        std::string remote = req.remote_address();
        if (_whitelist.find(remote) == _whitelist.end()) {
            GCBS_DEBUG("Incoming request from " + req.remote_address() + " has been blocked according to whitelist rule");
            req.reply(web::http::status_codes::NotFound);
        }
        GCBS_DEBUG("Incoming request from " + req.remote_address() + " has been accepted according to whitelist rule");
    } else {
        GCBS_DEBUG("Incoming request from " + req.remote_address());
    }

    std::vector<std::string> path = web::uri::split_path(web::uri::decode(req.relative_uri().path()));
    std::map<std::string, std::string> query_pars = web::uri::split_query(web::uri::decode(req.relative_uri().query()));
    //    std::for_each(path.begin(), path.end(), [](std::string s) { std::cout << s << std::endl; });
    //    std::for_each(query_pars.begin(), query_pars.end(), [](std::pair<std::string, std::string> s) { std::cout << s.first << ":" << s.second << std::endl; });

    if (!path.empty()) {
        if (path[0] == "file") {
            GCBS_DEBUG("POST /file " + req.remote_address());
            std::string fname;
            if (query_pars.find("name") != query_pars.end()) {
                fname = query_pars["name"];
            } else {
                fname = utils::generate_unique_filename();
            }
            fname = filesystem::join(_workdir, fname);
            if (filesystem::exists(fname)) {
                if (req.headers().has("Content-Length")) {
                    if (req.headers().content_length() == filesystem::file_size(fname)) {
                        req.reply(web::http::status_codes::OK);
                    } else {
                        req.reply(web::http::status_codes::Conflict);
                    }
                } else {
                    req.reply(web::http::status_codes::Conflict);
                }
            } else {
                filesystem::mkdir_recursive(filesystem::parent(fname));
                auto fstream = std::make_shared<concurrency::streams::ostream>();
                pplx::task<void> t = concurrency::streams::fstream::open_ostream(fname).then(
                    [fstream, &req](concurrency::streams::ostream outFile) {
                        *fstream = outFile;
                        concurrency::streams::istream indata = req.body();
                        uint32_t maxbytes = 1024 * 1024 * 8;  // Read up to 8MiB
                        while (!indata.is_eof()) {
                            indata.read(fstream->streambuf(), maxbytes).wait();
                        }
                        fstream->close().wait();
                    });
                t.wait();
                req.reply(web::http::status_codes::OK, fname.c_str(), "text/plain");
            }
        } else if (path[0] == "cube") {
            if (path.size() == 1) {
                GCBS_DEBUG("POST /cube");
                // we do not use cpprest JSON library here
                uint32_t id;
                std::string err;
                req.extract_string(true).then([&id, this, &err](std::string s) {
                                            std::shared_ptr<cube> c = cube_factory::instance()->create_from_json(json11::Json::parse(s, err));
                                            id = get_unique_id();
                                            _mutex_cubestore.lock();
                                            _cubestore.insert(std::make_pair(id, c));
                                            _mutex_cubestore.unlock();
                                        })
                    .wait();
                req.reply(web::http::status_codes::OK, std::to_string(id), "text/plain");
            } else if (path.size() == 4) {
                uint32_t cube_id = std::stoi(path[1]);
                uint32_t chunk_id = std::stoi(path[2]);

                std::string cmd = path[3];
                if (cmd == "start") {
                    GCBS_DEBUG("POST /cube/" + std::to_string(cube_id) + "/" + std::to_string(chunk_id) + "/start");

                    if (_cubestore.find(cube_id) == _cubestore.end()) {
                        req.reply(web::http::status_codes::NotFound, "ERROR in /POST /cube/{cube_id}/{chunk_id}/start: cube is not available", "text/plain");
                    } else if (chunk_id >= _cubestore[cube_id]->count_chunks()) {
                        req.reply(web::http::status_codes::NotFound, "ERROR in /POST /cube/{cube_id}/{chunk_id}/start: invalid chunk_id given", "text/plain");
                    }

                    // if already in queue, executing, or finished, do not compute again
                    else if (server_chunk_cache::instance()->has(std::make_pair(cube_id, chunk_id)) ||
                             _chunk_read_requests_set.find(std::make_pair(cube_id, chunk_id)) != _chunk_read_requests_set.end() ||
                             _chunk_read_executing.find(std::make_pair(cube_id, chunk_id)) != _chunk_read_executing.end()) {
                        req.reply(web::http::status_codes::OK);
                    } else {
                        _mutex_worker_threads.lock();
                        if (_worker_threads.size() < config::instance()->get_server_worker_threads_max()) {
                            // immediately start a thread

                            _worker_threads.push_back(std::thread([this, cube_id, chunk_id]() {
                                _mutex_chunk_read_requests.lock();
                                _chunk_read_requests.push_back(std::make_pair(cube_id, chunk_id));
                                _chunk_read_requests_set.insert(std::make_pair(cube_id, chunk_id));
                                _mutex_chunk_read_requests.unlock();
                                while (1) {
                                    _mutex_chunk_read_requests.lock();
                                    while (!_chunk_read_requests.empty()) {
                                        _mutex_chunk_read_executing.lock();
                                        uint32_t xcube_id = _chunk_read_requests.front().first;
                                        uint32_t xchunk_id = _chunk_read_requests.front().second;
                                        _chunk_read_executing.insert(_chunk_read_requests.front());
                                        _chunk_read_requests_set.erase(_chunk_read_requests.front());
                                        _chunk_read_requests.pop_front();
                                        _mutex_chunk_read_requests.unlock();
                                        _mutex_chunk_read_executing.unlock();

                                        _mutex_cubestore.lock();
                                        std::shared_ptr<cube> c = _cubestore[xcube_id];
                                        _mutex_cubestore.unlock();

                                        std::shared_ptr<chunk_data> dat = c->read_chunk(xchunk_id);

                                        server_chunk_cache::instance()->add(std::make_pair(xcube_id, xchunk_id), dat);

                                        _mutex_chunk_read_executing.lock();
                                        _chunk_read_executing.erase(std::make_pair(xcube_id, xchunk_id));
                                        _mutex_chunk_read_executing.unlock();

                                        // notify waiting download requests
                                        if (_chunk_cond.find(std::make_pair(xcube_id, xchunk_id)) != _chunk_cond.end()) {
                                            _chunk_cond[std::make_pair(xcube_id, xchunk_id)].first.notify_all();
                                        }

                                        _mutex_chunk_read_requests.lock();
                                    }
                                    _mutex_chunk_read_requests.unlock();

                                    std::unique_lock<std::mutex> lock(_mutex_worker_cond);
                                    _worker_cond.wait(lock);
                                }
                            }));

                            _mutex_worker_threads.unlock();
                            req.reply(web::http::status_codes::OK);
                        }

                        else {
                            _mutex_worker_threads.unlock();
                            std::lock_guard<std::mutex> lock(_mutex_worker_cond);
                            // add to queue
                            _mutex_chunk_read_requests.lock();
                            _chunk_read_requests.push_back(std::make_pair(cube_id, chunk_id));
                            _mutex_chunk_read_requests.unlock();
                            _worker_cond.notify_all();
                            req.reply(web::http::status_codes::OK);
                        }
                    }
                } else {
                    req.reply(web::http::status_codes::NotFound);
                }
            } else {
                req.reply(web::http::status_codes::NotFound);
            }

        } else {
            req.reply(web::http::status_codes::NotFound);
        }
    } else {
        req.reply(web::http::status_codes::NotFound);
    }
}

void gdalcubes_server::handle_head(web::http::http_request req) {
    if (!_whitelist.empty()) {
        std::string remote = req.remote_address();
        if (_whitelist.find(remote) == _whitelist.end()) {
            GCBS_DEBUG("Incoming request from " + req.remote_address() + " has been blocked according to whitelist rule");
            req.reply(web::http::status_codes::NotFound);
        }
        GCBS_DEBUG("Incoming request from " + req.remote_address() + " has been accepted according to whitelist rule");
    } else {
        GCBS_DEBUG("Incoming request from " + req.remote_address());
    }

    std::vector<std::string> path = web::uri::split_path(web::uri::decode(req.relative_uri().path()));
    std::map<std::string, std::string> query_pars = web::uri::split_query(web::uri::decode(req.relative_uri().query()));
    if (!path.empty() && path[0] == "file") {
        GCBS_DEBUG("HEAD /file");
        std::string fname;
        if (query_pars.find("name") != query_pars.end()) {
            fname = query_pars["name"];
            fname = filesystem::join(_workdir, fname);
            if (filesystem::exists(fname)) {
                if (query_pars.find("size") != query_pars.end()) {
                    if (std::stoi(query_pars["size"]) == (int)filesystem::file_size(fname)) {
                        req.reply(web::http::status_codes::OK);  // File exists and has the same size
                    } else {
                        req.reply(web::http::status_codes::Conflict);  // File exists but has different size
                    }
                } else {
                    req.reply(web::http::status_codes::BadRequest);  // File exists but Content Length missing for comparison
                }
            } else {
                req.reply(web::http::status_codes::NoContent);  // File does not exist yet
            }

        } else {
            req.reply(web::http::status_codes::BadRequest);  // no name given in request
        }
    } else {
        req.reply(web::http::status_codes::NotFound);  //
    }
}

}  // namespace gdalcubes

void print_usage() {
    std::cout << "Usage: gdalcubes_server [options]" << std::endl;
    std::cout << std::endl;
    std::cout << "Run gdalcubes_server and wait for incoming HTTP requests."
              << std::endl;
    std::cout << std::endl;
    std::cout << "Options:" << std::endl;
    std::cout << "  -b, --basepath              Base path for all API endpoints, defaults to /gdalcubes/api/" << std::endl;
    std::cout << "  -p, --port                  The port where gdalcubes_server is listening, defaults to 1111" << std::endl;
    std::cout << "  -t, --worker_threads        Maximum number of threads perfoming chunk reads, defaults to 1" << std::endl;
    std::cout << "  -D, --dir                   Working directory where files are stored, defaults to {TEMPDIR}/gdalcubes" << std::endl;
    std::cout << "      --ssl                   Use HTTPS (currently not implemented)" << std::endl;
    std::cout << "  -w, --whitelist             Optional path to a whitelist text file with a list of acceptable clients" << std::endl;
    std::cout << "  -d, --debug                 Print debug messages" << std::endl;
    std::cout << std::endl;
}

using namespace gdalcubes;

std::unique_ptr<gdalcubes_server> server;
int main(int argc, char* argv[]) {
    config::instance()->gdalcubes_init();

    namespace po = boost::program_options;
    // see https://stackoverflow.com/questions/15541498/how-to-implement-subcommands-using-boost-program-options

    po::options_description global_args("Options");
    global_args.add_options()("help,h", "")("version", "")("debug,d", "")("basepath,b", po::value<std::string>()->default_value("/gdalcubes/api"), "")("port,p", po::value<uint16_t>()->default_value(1111), "")("ssl", "")("worker_threads,t", po::value<uint16_t>()->default_value(1), "")("dir,D", po::value<std::string>()->default_value((filesystem::join(filesystem::get_tempdir(), "gdalcubes")), ""))("whitelist,w", po::value<std::string>(), "");

    po::variables_map vm;

    bool ssl;
    std::set<std::string> whitelist;
    try {
        po::parsed_options parsed = po::command_line_parser(argc, argv).options(global_args).allow_unregistered().run();
        po::store(parsed, vm);
        if (vm.count("version")) {
            version_info v = config::instance()->get_version_info();
            std::cout << "gdalcubes_server " << v.VERSION_MAJOR << "." << v.VERSION_MINOR << "." << v.VERSION_PATCH << " (" << v.GIT_COMMIT << ") built on " << v.BUILD_DATE << " " << v.BUILD_TIME << std::endl;
            return 0;
        }
        if (vm.count("help")) {
            print_usage();
            return 0;
        }
        if (vm.count("whitelist")) {
            std::string line;
            std::ifstream infile(vm["whitelist"].as<std::string>());
            while (std::getline(infile, line))
                whitelist.insert(line);
        }
        if (vm.count("debug")) {
            config::instance()->set_error_handler(error_handler::error_handler_debug_server);
        }

        ssl = vm.count("ssl") > 0;
    } catch (...) {
        std::cout << "ERROR in gdalcubes_server: cannot parse arguments." << std::endl;
        return 1;
    }

    std::unique_ptr<gdalcubes_server> srv = std::unique_ptr<gdalcubes_server>(
        new gdalcubes_server("0.0.0.0", vm["port"].as<uint16_t>(), vm["basepath"].as<std::string>(), ssl,
                             vm["dir"].as<std::string>(), whitelist));

    config::instance()->set_server_worker_threads_max(vm["worker_threads"].as<uint16_t>());

    srv->open().wait();
    std::cout << "gdalcubes_server waiting for incoming HTTP requests on " << srv->get_service_url() << "." << std::endl;
    std::cout << "Press ENTER to exit" << std::endl;
    std::string line;
    std::getline(std::cin, line);

    srv->close().wait();
    config::instance()->gdalcubes_cleanup();
    return 0;
}
