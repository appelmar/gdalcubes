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

#include "swarm.h"

#include <fstream>
#include <thread>

namespace gdalcubes {

size_t post_file_read_callback(char *buffer, size_t size, size_t nitems, void *userdata) {
    std::ifstream *is = ((std::ifstream *)userdata);
    is->read(buffer, size * nitems);  // or use readsome() ?
    return is->gcount();
}

std::shared_ptr<gdalcubes_swarm> gdalcubes_swarm::from_txtfile(std::string path) {
    std::ifstream file(path);
    std::vector<std::string> urllist;
    std::string line;
    while (std::getline(file, line)) {
        // TODO: check if valid URL
        if (!line.empty()) urllist.push_back(line);
    }

    file.close();
    return gdalcubes_swarm::from_urls(urllist);
}

void gdalcubes_swarm::post_file(std::string path, uint16_t server_index) {
    if (_server_handles[server_index]) {
        curl_easy_reset(_server_handles[server_index]);

        // Step 1: send a HEAD HTTP request to check whether the file already exists on the server
        curl_easy_setopt(_server_handles[server_index], CURLOPT_URL, (_server_uris[server_index] + "/file" + "?name=" + path + "&size=" + std::to_string(filesystem::file_size(path))).c_str());
        curl_easy_setopt(_server_handles[server_index], CURLOPT_CUSTOMREQUEST, "HEAD");
        curl_easy_setopt(_server_handles[server_index], CURLOPT_VERBOSE, config::instance()->get_swarm_curl_verbose() ? 1L : 0L);

        CURLcode res = curl_easy_perform(_server_handles[server_index]);
        long response_code;
        if (res != CURLE_OK) {
            throw std::string("ERROR in gdalcubes_swarm::post_file(): HEAD /file?name='" + path + "' to '" + _server_uris[server_index] + "' failed");
        } else {
            curl_easy_getinfo(_server_handles[server_index], CURLINFO_RESPONSE_CODE, &response_code);
        }

        if (response_code == 200) {
            // already exists on server
            return;
        } else if (response_code == 204 || response_code == 409) {
            // file does not exist or has different size

            // Step 2: if not upload

            std::ifstream is(path, std::ifstream::in | std::ifstream::binary);
            if (!is.is_open()) {
                throw std::string("ERROR in gdalcubes_swarm::push_execution_context(): cannot open file '" + path + "'");
            }

            curl_easy_setopt(_server_handles[server_index], CURLOPT_URL, (_server_uris[server_index] + "/file" + "?name=" + path).c_str());

            //struct curl_slist *header = NULL;
            //header = curl_slist_append(header, "Content-Type: application/octet-stream");
            // curl_easy_setopt(curl, CURLOPT_HTTPHEADER, header);

            curl_easy_setopt(_server_handles[server_index], CURLOPT_UPLOAD, 1L);
            curl_easy_setopt(_server_handles[server_index], CURLOPT_POST, 1L);
            curl_easy_setopt(_server_handles[server_index], CURLOPT_CUSTOMREQUEST, "POST");  // for whatever reason all requests are PUT if not set

            // This will most likely work only as long as curl_easy_perform is synchronous, otherwise &is might become invalid
            curl_easy_setopt(_server_handles[server_index], CURLOPT_READDATA, &is);
            curl_easy_setopt(_server_handles[server_index], CURLOPT_READFUNCTION, &post_file_read_callback);

            curl_easy_setopt(_server_handles[server_index], CURLOPT_INFILESIZE_LARGE, (curl_off_t)filesystem::file_size(path));

            curl_easy_setopt(_server_handles[server_index], CURLOPT_VERBOSE, config::instance()->get_swarm_curl_verbose() ? 1L : 0L);

            CURLcode res = curl_easy_perform(_server_handles[server_index]);
            if (res != CURLE_OK) {
                if (is.is_open()) is.close();
                throw std::string("ERROR in gdalcubes_swarm::post_file(): uploading '" + path + "' to '" + _server_uris[server_index] + "' failed");
            }
            if (is.is_open()) is.close();
        } else {
            throw std::string("ERROR in gdalcubes_swarm::post_file(): HEAD /file?name='" + path + "' to '" + _server_uris[server_index] + "' returned HTTP code " + std::to_string(response_code));
        }
    }
}

void gdalcubes_swarm::push_execution_context(bool recursive) {
    std::string p = filesystem::get_working_dir();

    std::vector<std::string> file_list;
    if (recursive) {
        filesystem::iterate_directory_recursive(p, [&file_list](const std::string &p) {
            if (filesystem::is_regular_file(p)) {
                file_list.push_back(p);
            }
        });
    } else {
        filesystem::iterate_directory(p, [&file_list](const std::string &p) {
            if (filesystem::is_regular_file(p)) {
                file_list.push_back(p);
            }
        });
    }

    // TODO: parallel upload
    for (auto it_f = file_list.begin(); it_f != file_list.end(); ++it_f) {
        for (uint16_t is = 0; is < _server_uris.size(); ++is) {
            post_file(*it_f, is);
        }
    }
}

size_t post_cube_read_callback(char *buffer, size_t size, size_t nitems, void *userdata) {
    std::istringstream *is = (std::istringstream *)userdata;
    is->readsome(buffer, size * nitems);  // or use readsome() ?
    return is->gcount();
}

size_t post_cube_write_callback(char *buffer, size_t size, size_t nitems, void *userdata) {
    std::vector<char> *x = (std::vector<char> *)userdata;
    x->insert(x->end(), buffer, buffer + (size * nitems));
    return size * nitems;
}

uint32_t gdalcubes_swarm::post_cube(std::string json, uint16_t server_index) {
    std::istringstream is(json);
    if (_server_handles[server_index]) {
        curl_easy_reset(_server_handles[server_index]);
        curl_easy_setopt(_server_handles[server_index], CURLOPT_URL, (_server_uris[server_index] + "/cube").c_str());

        // TODO: Disable Expect: 100-Continue header?

        curl_easy_setopt(_server_handles[server_index], CURLOPT_UPLOAD, 1L);
        curl_easy_setopt(_server_handles[server_index], CURLOPT_POST, 1L);
        curl_easy_setopt(_server_handles[server_index], CURLOPT_CUSTOMREQUEST, "POST");  // for whatever reason all requests are PUT if not set

        struct curl_slist *header = NULL;
        header = curl_slist_append(header, "Content-Type: application/json");
        header = curl_slist_append(header, "Expect:");
        curl_easy_setopt(_server_handles[server_index], CURLOPT_HTTPHEADER, header);

        curl_easy_setopt(_server_handles[server_index], CURLOPT_READFUNCTION, &post_cube_read_callback);

        // This will most likely work only as long as curl_easy_perform is synchronous, otherwise &is might become invalid
        curl_easy_setopt(_server_handles[server_index], CURLOPT_READDATA, &is);

        //curl_easy_setopt(_server_handles[server_index], CURLOPT_INFILESIZE, json.size() + 1);

        curl_easy_setopt(_server_handles[server_index], CURLOPT_WRITEFUNCTION, &post_cube_write_callback);

        // This will most likely work only as long as curl_easy_perform is synchronous, otherwise &response_body_bytes might become invalid
        std::vector<char> response_body_bytes;
        curl_easy_setopt(_server_handles[server_index], CURLOPT_WRITEDATA, &response_body_bytes);

        curl_easy_setopt(_server_handles[server_index], CURLOPT_VERBOSE, config::instance()->get_swarm_curl_verbose() ? 1L : 0L);

        CURLcode res = curl_easy_perform(_server_handles[server_index]);
        if (res != CURLE_OK) {
            throw std::string("ERROR in gdalcubes_swarm::post_cube(): POST /cube to '" + _server_uris[server_index] + "' failed");
        }
        std::string response_body(response_body_bytes.begin(), response_body_bytes.end());
        return std::stoi(response_body);
    } else {
        throw std::string("ERROR in gdalcubes_swarm::post_cube(): no connection with given server index available");
    }
}

void gdalcubes_swarm::push_cube(std::shared_ptr<cube> c) {
    _cube = c;  // Does this require a mutex?
    std::string json = c->make_constructible_json().dump();

    _cube_ids.clear();

    // TODO: parallelize requests
    for (uint16_t is = 0; is < _server_uris.size(); ++is) {
        _cube_ids.push_back(post_cube(json, is));
    }
}

void gdalcubes_swarm::post_start(uint32_t chunk_id, uint16_t server_index) {
    if (_server_handles[server_index]) {
        // TODO: URLencode?
        curl_easy_reset(_server_handles[server_index]);
        curl_easy_setopt(_server_handles[server_index], CURLOPT_URL, (_server_uris[server_index] + "/cube/" + std::to_string(_cube_ids[server_index]) + "/" + std::to_string(chunk_id) + "/start").c_str());
        curl_easy_setopt(_server_handles[server_index], CURLOPT_POST, 1L);
        curl_easy_setopt(_server_handles[server_index], CURLOPT_CUSTOMREQUEST, "POST");
        curl_easy_setopt(_server_handles[server_index], CURLOPT_VERBOSE, config::instance()->get_swarm_curl_verbose() ? 1L : 0L);

        CURLcode res = curl_easy_perform(_server_handles[server_index]);
        if (res != CURLE_OK) {
            throw std::string("ERROR in gdalcubes_swarm::post_start(): POST /cube/{cube_id}/{chunk_id}/start to '" + _server_uris[server_index] + "' failed");
        }
        return;
    }
}

size_t get_download_callback(char *buffer, size_t size, size_t nitems, void *userdata) {
    std::vector<char> *x = (std::vector<char> *)userdata;
    x->insert(x->end(), buffer, buffer + (size * nitems));
    return size * nitems;
}

std::shared_ptr<chunk_data> gdalcubes_swarm::get_download(uint32_t chunk_id, uint16_t server_index) {
    if (_server_handles[server_index]) {
        // TODO: URLencode?
        curl_easy_reset(_server_handles[server_index]);
        curl_easy_setopt(_server_handles[server_index], CURLOPT_URL, (_server_uris[server_index] + "/cube/" + std::to_string(_cube_ids[server_index]) + "/" + std::to_string(chunk_id) + "/download").c_str());
        curl_easy_setopt(_server_handles[server_index], CURLOPT_HTTPGET, 1L);
        curl_easy_setopt(_server_handles[server_index], CURLOPT_CUSTOMREQUEST, "GET");
        curl_easy_setopt(_server_handles[server_index], CURLOPT_WRITEFUNCTION, &get_download_callback);

        // This will most likely work only as long as curl_easy_perform is synchronous, otherwise &response_body_bytes might become invalid
        std::vector<char> response_body_bytes;
        curl_easy_setopt(_server_handles[server_index], CURLOPT_WRITEDATA, &response_body_bytes);

        curl_easy_setopt(_server_handles[server_index], CURLOPT_VERBOSE, config::instance()->get_swarm_curl_verbose() ? 1L : 0L);

        CURLcode res = curl_easy_perform(_server_handles[server_index]);
        if (res != CURLE_OK) {
            throw std::string("ERROR in gdalcubes_swarm::get_download(): GET /cube/{cube_id}/{chunk_id}/download to '" + _server_uris[server_index] + "' failed");
        }

        // construct chunk_data from byte vector
        std::shared_ptr<chunk_data> out = std::make_shared<chunk_data>();
        std::array<uint32_t, 4> size = {((uint32_t *)response_body_bytes.data())[0], ((uint32_t *)response_body_bytes.data())[1], ((uint32_t *)response_body_bytes.data())[2], ((uint32_t *)response_body_bytes.data())[3]};
        out->size(size);
        if (size[0] * size[1] * size[2] * size[3] > 0) {
            out->buf(std::malloc(out->size()[0] * out->size()[1] * out->size()[2] * out->size()[3] * sizeof(double)));
            std::copy(response_body_bytes.begin() + sizeof(std::array<uint32_t, 4>), response_body_bytes.end(),
                      (char *)out->buf());
        }

        return out;
    } else {
        throw std::string("ERROR in gdalcubes_swarm::get_download(): no connection with given server index available");
    }
}

void gdalcubes_swarm::apply(std::shared_ptr<cube> c, std::function<void(chunkid_t, std::shared_ptr<chunk_data>, std::mutex &)> f) {
    uint32_t nthreads = 1;
    // try whether default chunk processor is multithread and read number of threads if successful
    if (std::dynamic_pointer_cast<chunk_processor_multithread>(config::instance()->get_default_chunk_processor())) {
        nthreads = std::dynamic_pointer_cast<chunk_processor_multithread>(config::instance()->get_default_chunk_processor())->get_threads();
    }

    push_execution_context(false);
    push_cube(c);

    // TODO: what if servers fail?
    std::map<uint16_t, std::vector<uint32_t>> chunk_distr;
    for (uint32_t i = 0; i < c->count_chunks(); ++i) {
        post_start(i, i % _server_uris.size());  // start computations on remote servers
        chunk_distr[i % _server_uris.size()].push_back(i);
    }
    std::mutex mutex;
    std::vector<std::thread> workers;
    for (uint16_t it = 0; it < nthreads; ++it) {
        workers.push_back(std::thread([this, c, &chunk_distr, &f, it, &mutex](void) {
            // One remote server is processed by a single thread for now, changing this will require to using the multi interface of curl
            for (uint32_t iserver = it; iserver < _server_uris.size(); iserver += _server_uris.size()) {
                for (uint32_t ichunk = 0; ichunk < chunk_distr[iserver].size(); ++ichunk) {
                    std::shared_ptr<chunk_data> dat = get_download(chunk_distr[iserver][ichunk], iserver);
                    f(chunk_distr[iserver][ichunk], dat, mutex);
                }
            }
        }));
    }
    for (uint16_t it = 0; it < nthreads; ++it) {
        workers[it].join();
    }
}

}  // namespace gdalcubes

#endif  // GDALCUBES_NO_SWARM
