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

#include "view.h"
#include "filesystem.h"

cube_view cube_view::read(nlohmann::json j) {
    cube_view v;

    v._dt = duration::from_string(j.at("time").at("dt").get<std::string>());

    std::string st0 = j.at("time").at("t0").get<std::string>();
    std::string st1 = j.at("time").at("t1").get<std::string>();

    v._t0 = datetime::from_string(st0);
    v._t1 = datetime::from_string(st1);

    // No matter how the start and end datetime are given, use the unit of the datetime interval!
    v._t0.unit() = v._dt.dt_unit;
    v._t1.unit() = v._dt.dt_unit;

    if (j.count("space")) {
        v._win.left = j.at("space").at("left").get<double>();
        v._win.right = j.at("space").at("right").get<double>();
        v._win.top = j.at("space").at("top").get<double>();
        v._win.bottom = j.at("space").at("bottom").get<double>();

        if (j.at("space").count("nx") && (j.at("space").count("ny"))) {
            v._nx = j.at("space").at("nx").get<uint32_t>();
            v._ny = j.at("space").at("ny").get<uint32_t>();
        }
        // if only one of nx and ny is given, use the aspect ratio of the spatial extent to derive the other
        else if (j.at("space").count("nx") && !(j.at("space").count("ny"))) {
            v._nx = j.at("space").at("nx").get<uint32_t>();
            v._ny = (uint32_t)((double)v._nx * (v._win.top - v._win.bottom) / (v._win.right - v._win.left));
        } else if (!j.at("space").count("nx") && (j.at("space").count("ny"))) {
            v._ny = j.at("space").at("ny").get<uint32_t>();
            v._nx = (uint32_t)((double)v._ny * (v._win.right - v._win.left) / (v._win.top - v._win.bottom));
        } else {
            throw std::string("ERROR in cube_view::read_json_string(): at least one of nx or ny must be given.");
        }

        v._proj = j.at("space").at("proj").get<std::string>();

    } else if (j.count("tile")) {
        const double EARTH_RADIUS_METERS = 6378137;
        const uint16_t TILE_SIZE_PX = 256;
        const double EARTH_CIRCUMFERENCE_METERS = 2 * M_PI * EARTH_RADIUS_METERS;
        const double WEBMERCATOR_BOUNDS_LEFT = -EARTH_CIRCUMFERENCE_METERS / 2.0;  // -20037508.342789244
        //const double WEBMERCATOR_BOUNDS_LOWER = -EARTH_CIRCUMFERENCE_METERS / 2.0;  // -20037508.342789244
        //const double WEBMERCATOR_BOUNDS_RIGHT = EARTH_CIRCUMFERENCE_METERS / 2.0;   // 20037508.342789244
        const double WEBMERCATOR_BOUNDS_UPPER = EARTH_CIRCUMFERENCE_METERS / 2.0;  // 20037508.342789244

        uint32_t x = j.at("tile").at("x").get<uint32_t>();
        uint32_t y = j.at("tile").at("y").get<uint32_t>();
        uint32_t z = j.at("tile").at("z").get<uint32_t>();

        bounds_2d<double> win;
        win.left = WEBMERCATOR_BOUNDS_LEFT + x * EARTH_CIRCUMFERENCE_METERS / std::pow(2, z);
        win.right = WEBMERCATOR_BOUNDS_LEFT + (x + 1) * EARTH_CIRCUMFERENCE_METERS / std::pow(2, z);
        win.top = WEBMERCATOR_BOUNDS_UPPER - y * EARTH_CIRCUMFERENCE_METERS / std::pow(2, z);
        win.bottom = WEBMERCATOR_BOUNDS_UPPER - (y + 1) * EARTH_CIRCUMFERENCE_METERS / std::pow(2, z);

        v._win = win;
        v._nx = TILE_SIZE_PX;
        v._ny = TILE_SIZE_PX;

        v._proj = "EPSG:3857";
    } else {
        throw std::string("ERROR in cube_view::read(): expected either 'space' or 'tile' in JSON cube view");
    }

    if (j.count("resampling") == 0) {
        v._resampling = resampling::resampling_type::RSMPL_NEAR;
    } else {
        v._resampling = resampling::from_string(j.at("resampling").get<std::string>());
    }

    if (j.count("aggregation") == 0) {
        v._aggregation = aggregation::aggregation_type::AGG_NONE;
    } else {
        v._aggregation = aggregation::from_string(j.at("aggregation").get<std::string>());
    }

    return v;
}

cube_view cube_view::read_json(std::string filename) {
    if (!filesystem::exists(filename))
        throw std::string("ERROR in cube_view::read_json(): image_collection_cube view file does not exist.");
    std::ifstream i(filename);
    nlohmann::json j;
    i >> j;
    return read(j);
}

cube_view cube_view::read_json_string(std::string str) {
    std::istringstream i(str);
    nlohmann::json j;
    i >> j;
    cube_view v;

    return read(j);
}

void cube_view::write_json(std::string filename) {
    nlohmann::json j = nlohmann::json{
        {"space", {{"nx", _nx}, {"ny", _ny}, {"left", _win.left}, {"right", _win.right}, {"top", _win.top}, {"bottom", _win.bottom}, {"proj", _proj}}},
        {"time", {{"dt", dt().to_string()}, {"t0", _t0.to_string()}, {"t1", _t1.to_string()}}},
        {"aggregation", aggregation::to_string(_aggregation)},
        {"resampling", resampling::to_string(_resampling)}};

    std::ofstream o(filename, std::ofstream::out);
    if (!o.good()) {
        throw std::string("ERROR in cube_view::write_json(): cannot write to file.");
    }

    o << std::setw(4) << j << std::endl;
    o.close();
}

std::string cube_view::write_json_string() {
    nlohmann::json j = nlohmann::json{
        {"space", {{"nx", _nx}, {"ny", _ny}, {"left", _win.left}, {"right", _win.right}, {"top", _win.top}, {"bottom", _win.bottom}, {"proj", _proj}}},
        {"time", {{"dt", dt().to_string()}, {"t0", _t0.to_string()}, {"t1", _t1.to_string()}}},
        {"aggregation", aggregation::to_string(_aggregation)},
        {"resampling", resampling::to_string(_resampling)}};
    std::ostringstream o;
    o << j;
    return o.str();
}