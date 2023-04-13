/*
   Copyright 2019 Marius Appel <marius.appel@uni-muenster.de>

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

#ifndef IMAGE_COLLECTION_OPS_H
#define IMAGE_COLLECTION_OPS_H

#include <cstdint> // 2023-01-12: GCC 13 compatibility
#include "image_collection.h"

namespace gdalcubes {

/**
     * Batch processing operations over all GDAL datasets of a collection
     */
class image_collection_ops {
   public:
    static void translate_gtiff(std::shared_ptr<gdalcubes::image_collection> in, std::string out_dir, uint16_t nthreads = 1, bool overwrite = true, std::vector<std::string> creation_options = {});

    static void translate_cog(std::shared_ptr<gdalcubes::image_collection> in, std::string out_dir, uint16_t nthreads = 1, bool overwrite = true, std::vector<std::string> creation_options = {});

    static void create_overviews(std::shared_ptr<image_collection> in, std::vector<int> levels = std::vector<int>{2, 4, 8, 16, 32}, std::string resampling = "NEAREST", uint16_t nthreads = 1);
};

}  // namespace gdalcubes

#endif
