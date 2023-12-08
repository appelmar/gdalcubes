/*
    MIT License

    Copyright (c) 2022 Marius Appel <marius.appel@hs-bochum.de>

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

#ifndef GDALCUBES_EXTRACT_GEOM_H
#define GDALCUBES_EXTRACT_GEOM_H

#include "cube.h"

namespace gdalcubes {

/**
 * @brief A data cube where that extracts pixels from vector geometries
 * resulting in a sparse cube
 */
class extract_geom : public cube {
   public:
    /**
     * @brief Create a data cube from another cube where pixels are selected by vector geometries
     * @note This static creation method should preferably be used instead of the constructors as
     * the constructors will not set connections between cubes properly.
     * @param in source data cube
     * @param ogr_dataset Path to an existing OGR dataset
     * @param time_column optional name of the column in ogr_dataset containing time information
     * @param, datetime optional datetime strings of features, length must be identical to the number of features in ogr_dataset; only one of time_column and datetime should be provided
     * @return a shared pointer to the created data cube instance
     */
    static std::shared_ptr<extract_geom> create(std::shared_ptr<cube> in, std::string ogr_dataset, std::string time_column = "", std::string ogr_layer = "") {
        std::shared_ptr<extract_geom> out = std::make_shared<extract_geom>(in, ogr_dataset, time_column,  ogr_layer);
        return out;
    }

   public:
    extract_geom(std::shared_ptr<cube> in, std::string ogr_dataset, std::string time_column = "", std::string ogr_layer = "");
    ~extract_geom() {
        if (!_ogr_dataset.empty()) {
            if (_in_ogr_was_transformed) {
                filesystem::remove(_ogr_dataset);
            }
        }
        if (!_chunkmask_dataset.empty()) {
            if (filesystem::exists(_chunkmask_dataset)) {
                filesystem::remove(_chunkmask_dataset);
            }
        }
    }

    std::shared_ptr<chunk_data> read_chunk(chunkid_t id) override;

    json11::Json make_constructible_json() override {
        json11::Json::object out;
        out["cube_type"] = "extract";
        out["in_cube"] = _in_cube->make_constructible_json();
        out["ogr_dataset"] = _in_ogr_dataset;
        out["time_column"] = _in_time_column;
        out["ogr_layer"] = _in_ogr_layer;
        return out;
    }

   private:

    // Input variables used to create the object
    std::shared_ptr<cube> _in_cube;
    std::string _in_ogr_dataset;
    std::string _in_time_column;
    std::string _in_ogr_layer;

    // State variables
    std::string _ogr_dataset; // can be different from _in_ogr_dataset if reprojection / transformation was needed
    std::string _ogr_layer;
    std::string _fid_column;
    bool _in_ogr_was_transformed;
    bool _is_point;
    std::string _chunkmask_dataset;
};

}  // namespace gdalcubes

#endif  // GDALCUBES_EXTRACT_GEOM_H
