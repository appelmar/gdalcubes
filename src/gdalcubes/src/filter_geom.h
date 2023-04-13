/*
    MIT License

    Copyright (c) 2020 Marius Appel <marius.appel@uni-muenster.de>

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

#ifndef GDALCUBES_FILTER_GEOM_H
#define GDALCUBES_FILTER_GEOM_H

#include "cube.h"

namespace gdalcubes {

/**
     * @brief A data cube where space is cropped by a polygon
     */
class filter_geom_cube : public cube {
   public:
    /**
         * @brief Create a data cube from another cube where space is cropped by a polygon
         * @note This static creation method should preferably be used instead of the constructors as
         * the constructors will not set connections between cubes properly.
         * @param in source data cube
         * @param wkt WKT of the crop polygon
         * @param srs spatial reference system in any form that is readable by GDAL
         * @return a shared pointer to the created data cube instance
         */
    static std::shared_ptr<filter_geom_cube> create(std::shared_ptr<cube> in, std::string wkt, std::string srs) {
        std::shared_ptr<filter_geom_cube> out = std::make_shared<filter_geom_cube>(in, wkt, srs);
        return out;
    }

   public:
    filter_geom_cube(std::shared_ptr<cube> in, std::string wkt, std::string srs);
    ~filter_geom_cube() {
        if (!_ogr_dataset.empty()) {
            if (filesystem::exists(_ogr_dataset)) {
                filesystem::remove(_ogr_dataset);
            }
        }
    }

    std::shared_ptr<chunk_data> read_chunk(chunkid_t id) override;

    json11::Json make_constructible_json() override {
        json11::Json::object out;
        out["cube_type"] = "filter_geom";
        out["in_cube"] = _in_cube->make_constructible_json();
        out["wkt"] = _wkt;
        out["srs"] = _srs;
        return out;
    }

   private:
    std::shared_ptr<cube> _in_cube;
    std::string _wkt;
    std::string _srs;

    std::string _ogr_dataset;  // filename of transformed input polygon
    uint32_t _min_chunk_x;
    uint32_t _max_chunk_x;
    uint32_t _min_chunk_y;
    uint32_t _max_chunk_y;
};
}  // namespace gdalcubes

#endif  //GDALCUBES_FILTER_GEOM_H
