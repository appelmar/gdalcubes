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

#include "cube_factory.h"

#include "apply_pixel.h"
#include "external/json.hpp"
#include "filesystem.h"
#include "image_collection_cube.h"
#include "join_bands.h"
#include "reduce.h"
#include "select_bands.h"
#include "stream.h"

std::shared_ptr<cube> cube_factory::create_from_json(nlohmann::json j) {
    // TODO: move the map to somwhere else but not in this function, alternatives could be a singleton implementation of this class
    std::map<std::string, std::function<std::shared_ptr<cube>(nlohmann::json&)>> cube_generators;
    cube_generators.insert(std::make_pair<std::string, std::function<std::shared_ptr<cube>(nlohmann::json&)>>(
        "reduce", [](nlohmann::json& j) {
            auto x = reduce_cube::create(create_from_json(j["in_cube"]), j["reducer"].get<std::string>());
            return x;
        }));
    cube_generators.insert(std::make_pair<std::string, std::function<std::shared_ptr<cube>(nlohmann::json&)>>(
        "select_bands", [](nlohmann::json& j) {
            auto x = select_bands_cube::create(create_from_json(j["in_cube"]), j["bands"].get<std::vector<std::string>>());
            return x;
        }));

    cube_generators.insert(std::make_pair<std::string, std::function<std::shared_ptr<cube>(nlohmann::json&)>>(
        "apply_pixel", [](nlohmann::json& j) {
            if (j.count("band_names") > 0) {
                auto x = apply_pixel_cube::create(create_from_json(j["in_cube"]), j["expr"].get<std::vector<std::string>>(), j["band_names"].get<std::vector<std::string>>());
                return x;
            } else {
                auto x = apply_pixel_cube::create(create_from_json(j["in_cube"]), j["expr"].get<std::vector<std::string>>());
                return x;
            }
        }));
    cube_generators.insert(std::make_pair<std::string, std::function<std::shared_ptr<cube>(nlohmann::json&)>>(
        "join_bands", [](nlohmann::json& j) {
            auto x = join_bands_cube::create(create_from_json(j["A"]), create_from_json(j["B"]));
            return x;
        }));

    cube_generators.insert(std::make_pair<std::string, std::function<std::shared_ptr<cube>(nlohmann::json&)>>(
        "stream", [](nlohmann::json& j) {
            auto x = stream_cube::create(create_from_json(j["in_cube"]), j["command"].get<std::string>());
            return x;
        }));

    cube_generators.insert(std::make_pair<std::string, std::function<std::shared_ptr<cube>(nlohmann::json&)>>(
        "image_collection", [](nlohmann::json& j) {
            if (!filesystem::exists(j["file"].get<std::string>())) {
                throw std::string("ERROR in cube_generators[\"image_collection\"](): image collection file does not exist.");
            }
            cube_view v = cube_view::read_json_string(j["view"].dump());
            auto x = image_collection_cube::create(j["file"].get<std::string>(), v);
            x->set_chunk_size(j["chunk_size"][0].get<uint32_t>(), j["chunk_size"][1].get<uint32_t>(), j["chunk_size"][2].get<uint32_t>());
            return x;
        }));

    if (!j.count("cube_type")) {
        throw std::string("ERROR in cube_factory::create_from_json(): invalid object, missing cube_type key.");
    }

    std::string cube_type = j["cube_type"];

    return (cube_generators[cube_type](j));  //recursive creation
}
