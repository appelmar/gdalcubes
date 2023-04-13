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

#include "view.h"

#include <fstream>

#include "filesystem.h"

namespace gdalcubes {

cube_view cube_view::read(json11::Json j) {
    cube_view v;

    v._dt = duration::from_string(j["time"]["dt"].string_value());

    std::string st0 = j["time"]["t0"].string_value();
    std::string st1 = j["time"]["t1"].string_value();

    v._t0 = datetime::from_string(st0);
    v._t1 = datetime::from_string(st1);

    // No matter how the start and end datetime are given, use the unit of the datetime interval!
    v._t0.unit(v._dt.dt_unit);
    v._t1.unit(v._dt.dt_unit);

    if (!j["space"].is_null()) {
        v._win.left = j["space"]["left"].number_value();
        v._win.right = j["space"]["right"].number_value();
        v._win.top = j["space"]["top"].number_value();
        v._win.bottom = j["space"]["bottom"].number_value();

        if (!j["space"]["nx"].is_null() && !j["space"]["ny"].is_null()) {
            v._nx = j["space"]["nx"].int_value();
            v._ny = j["space"]["ny"].int_value();
        }
        // if only one of nx and ny is given, use the aspect ratio of the spatial extent to derive the other
        else if (!j["space"]["nx"].is_null() && j["space"]["ny"].is_null()) {
            v._nx = j["space"]["nx"].int_value();
            v._ny = (uint32_t)((double)v._nx * (v._win.top - v._win.bottom) / (v._win.right - v._win.left));
        } else if (j["space"]["nx"].is_null() && !j["space"]["ny"].is_null()) {
            v._ny = j["space"]["ny"].int_value();
            v._nx = (uint32_t)((double)v._ny * (v._win.right - v._win.left) / (v._win.top - v._win.bottom));
        } else {
            throw std::string("ERROR in cube_view::read_json_string(): at least one of nx or ny must be given.");
        }

        v._srs = j["space"]["srs"].string_value();

    } else if (!j["tile"].is_null()) {
        const double EARTH_RADIUS_METERS = 6378137;
        const uint16_t TILE_SIZE_PX = 256;
        const double EARTH_CIRCUMFERENCE_METERS = 2 * M_PI * EARTH_RADIUS_METERS;
        const double WEBMERCATOR_BOUNDS_LEFT = -EARTH_CIRCUMFERENCE_METERS / 2.0;  // -20037508.342789244
        //const double WEBMERCATOR_BOUNDS_LOWER = -EARTH_CIRCUMFERENCE_METERS / 2.0;  // -20037508.342789244
        //const double WEBMERCATOR_BOUNDS_RIGHT = EARTH_CIRCUMFERENCE_METERS / 2.0;   // 20037508.342789244
        const double WEBMERCATOR_BOUNDS_UPPER = EARTH_CIRCUMFERENCE_METERS / 2.0;  // 20037508.342789244

        uint32_t x = j["tile"]["x"].int_value();
        uint32_t y = j["tile"]["y"].int_value();
        uint32_t z = j["tile"]["z"].int_value();

        bounds_2d<double> win;
        win.left = WEBMERCATOR_BOUNDS_LEFT + x * EARTH_CIRCUMFERENCE_METERS / std::pow(2, z);
        win.right = WEBMERCATOR_BOUNDS_LEFT + (x + 1) * EARTH_CIRCUMFERENCE_METERS / std::pow(2, z);
        win.top = WEBMERCATOR_BOUNDS_UPPER - y * EARTH_CIRCUMFERENCE_METERS / std::pow(2, z);
        win.bottom = WEBMERCATOR_BOUNDS_UPPER - (y + 1) * EARTH_CIRCUMFERENCE_METERS / std::pow(2, z);

        v._win = win;
        v._nx = TILE_SIZE_PX;
        v._ny = TILE_SIZE_PX;

        v._srs = "EPSG:3857";
    } else {
        throw std::string("ERROR in cube_view::read(): expected either 'space' or 'tile' in JSON cube view");
    }

    if (j["resampling"].is_null()) {
        v._resampling = resampling::resampling_type::RSMPL_NEAR;
    } else {
        v._resampling = resampling::from_string(j["resampling"].string_value());
    }

    if (j["aggregation"].is_null()) {
        v._aggregation = aggregation::aggregation_type::AGG_NONE;
    } else {
        v._aggregation = aggregation::from_string(j["aggregation"].string_value());
    }

    return v;
}

cube_view cube_view::read_json(std::string filename) {
    if (!filesystem::exists(filename))
        throw std::string("ERROR in cube_view::read_json(): image_collection_cube view file does not exist.");
    std::ifstream i(filename);

    std::stringstream buf;
    buf << i.rdbuf();

    std::string err;  // TODO: do something with error
    json11::Json j = json11::Json::parse(buf.str(), err);
    return read(j);
}

cube_view cube_view::read_json_string(std::string str) {
    std::istringstream i(str);
    std::string err;  // TODO: do something with error
    json11::Json j = json11::Json::parse(str, err);
    return read(j);
}

void cube_view::write_json(std::string filename) {
    json11::Json j = json11::Json::object{
        {"space", json11::Json::object{{"nx", (int)_nx}, {"ny", (int)_ny}, {"left", _win.left}, {"right", _win.right}, {"top", _win.top}, {"bottom", _win.bottom}, {"srs", _srs}}},
        {"time", json11::Json::object{{"dt", dt().to_string()}, {"t0", _t0.to_string()}, {"t1", _t1.to_string()}}},
        {"aggregation", aggregation::to_string(_aggregation)},
        {"resampling", resampling::to_string(_resampling)}};

    std::ofstream o(filename, std::ofstream::out);
    if (!o.good()) {
        throw std::string("ERROR in cube_view::write_json(): cannot write to file.");
    }

    o << std::setw(4) << j.dump() << std::endl;
    o.close();
}

std::string cube_view::write_json_string() {
    json11::Json j = json11::Json::object{
        {"space", json11::Json::object{{"nx", (int)_nx}, {"ny", (int)_ny}, {"left", _win.left}, {"right", _win.right}, {"top", _win.top}, {"bottom", _win.bottom}, {"srs", _srs}}},
        {"time", json11::Json::object{{"dt", dt().to_string()}, {"t0", _t0.to_string()}, {"t1", _t1.to_string()}}},
        {"aggregation", aggregation::to_string(_aggregation)},
        {"resampling", resampling::to_string(_resampling)}};
    std::ostringstream o;
    return j.dump();
}

}  // namespace gdalcubes