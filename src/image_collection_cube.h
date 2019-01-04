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
#ifndef IMAGE_COLLECTION_CUBE_H
#define IMAGE_COLLECTION_CUBE_H

#include "cube.h"

/**
 * @brief A data cube that reads data from an image collection
 *
 * An image collection cube is created from a cube view and reads data from an image collection.
 * The cube view defines the shape of the cube (size, extent, etc.) and automatically derives
 * which images of the collection are relevant for which chunks. To transform images to the shape of the cube,
 * [gdalwarp](https://www.gdal.org/gdalwarp.html) is applied on each image and images that fall into the same temporal slice
 * are aggregated with an aggregation function.
 *
 * @see image_collection
 * @see cube_view
 */
class image_collection_cube : public cube {
   public:
    /**
     * @brief Create a data cube from an image collection
     * @note This static creation method should preferably be used instead of the constructors as
     * the constructors will not set connections between cubes properly.
     * @param ic input image collection
     * @param v data cube view
     * @return a shared pointer to the created data cube instance
     */
    static std::shared_ptr<image_collection_cube> create(std::shared_ptr<image_collection> ic, cube_view v) {
        return std::make_shared<image_collection_cube>(ic, v);
    }

    /**
      * @brief Create a data cube from an image collection
      * @note This static creation method should preferably be used instead of the constructors as
      * the constructors will not set connections between cubes properly.
      * @param icfile filename of the input image collection
      * @param v data cube view
      * @return a shared pointer to the created data cube instance
      */
    static std::shared_ptr<image_collection_cube> create(std::string icfile, cube_view v) {
        return std::make_shared<image_collection_cube>(icfile, v);
    }

    /**
     * @brief Create a data cube from an image collection
     * @note This static creation method should preferably be used instead of the constructors as
     * the constructors will not set connections between cubes properly.
     * @param ic input image collection
     * @param vfile filename of the data cube view json description
     * @return a shared pointer to the created data cube instance
     */
    static std::shared_ptr<image_collection_cube> create(std::shared_ptr<image_collection> ic, std::string vfile) {
        return std::make_shared<image_collection_cube>(ic, vfile);
    }

    /**
     * @brief Create a data cube from an image collection
     * @note This static creation method should preferably be used instead of the constructors as
     * the constructors will not set connections between cubes properly.
     * @param icfile filename of the input image collection
     * @param vfile filename of the data cube view json description
     * @return a shared pointer to the created data cube instance
     */
    static std::shared_ptr<image_collection_cube> create(std::string icfile, std::string vfile) {
        return std::make_shared<image_collection_cube>(icfile, vfile);
    }

    /**
      * @brief Create a data cube from an image collection
      * This function will create a data cube from a given image collection and automatically derive a default data cube view.
      * @note This static creation method should preferably be used instead of the constructors as
      * the constructors will not set connections between cubes properly.
      * @param ic input image collection
      * @return a shared pointer to the created data cube instance
      */
    static std::shared_ptr<image_collection_cube> create(std::shared_ptr<image_collection> ic) {
        return std::make_shared<image_collection_cube>(ic);
    }

    /**
     * @brief Create a data cube from an image collection
     * This function will create a data cube from a given image collection and automatically derive a default data cube view.
     * @note This static creation method should preferably be used instead of the constructors as
     * the constructors will not set connections between cubes properly.
     * @param icfile filename of the input image collection
     * @return a shared pointer to the created data cube instance
     */
    static std::shared_ptr<image_collection_cube> create(std::string icfile) {
        return std::make_shared<image_collection_cube>(icfile);
    }

   public:
    image_collection_cube(std::shared_ptr<image_collection> ic, cube_view v);
    image_collection_cube(std::string icfile, cube_view v);
    image_collection_cube(std::shared_ptr<image_collection> ic, std::string vfile);
    image_collection_cube(std::string icfile, std::string vfile);
    image_collection_cube(std::shared_ptr<image_collection> ic);
    image_collection_cube(std::string icfile);

   public:
    ~image_collection_cube() {}

    inline const std::shared_ptr<image_collection> collection() { return _collection; }
    inline std::shared_ptr<cube_view> view() { return std::dynamic_pointer_cast<cube_view>(_st_ref); }

    std::string to_string() override;

    /**
     * @brief Select bands by names
     * @param bands vector of bands to be considered in the cube, if empty, all bands will be selected
     */
    void select_bands(std::vector<std::string> bands);

    /**
     * @brief Select bands by indexes
     * @param bands vector of bands to be considered in the cube, if empty, all bands will be selected
     */
    void select_bands(std::vector<uint16_t> bands);

    std::shared_ptr<chunk_data> read_chunk(chunkid_t id) override;

    // image_collection_cube is the only class that supports changing chunk sizes from outside!
    // This is important for e.g. streaming.
    void set_chunk_size(uint32_t t, uint32_t y, uint32_t x) {
        _chunk_size = {t, y, x};
    }

    nlohmann::json make_constructible_json() override {
        if (_collection->is_temporary()) {
            throw std::string("ERROR in image_collection_cube::make_constructible_json(): image collection is temporary, please export as file using write() first.");
        }
        nlohmann::json out;
        out["cube_type"] = "image_collection";
        out["chunk_size"] = {_chunk_size[0], _chunk_size[1], _chunk_size[2]};
        out["view"] = nlohmann::json::parse(std::dynamic_pointer_cast<cube_view>(_st_ref)->write_json_string());
        out["file"] = _collection->get_filename();
        return out;
    }

    static cube_view default_view(std::shared_ptr<image_collection> ic);

   protected:
    virtual void set_st_reference(std::shared_ptr<cube_st_reference> stref) override {
        // _st_ref has always type cube_view, whereas stref can have types st_reference or cube_view

        // copy fields from st_reference type
        _st_ref->win() = stref->win();
        _st_ref->proj() = stref->proj();
        _st_ref->ny() = stref->ny();
        _st_ref->nx() = stref->nx();
        _st_ref->t0() = stref->t0();
        _st_ref->t1() = stref->t1();
        _st_ref->dt() = stref->dt();

        // if view: copy aggregation and resampling
        std::shared_ptr<cube_view> v = std::dynamic_pointer_cast<cube_view>(stref);
        if (v) {
            std::dynamic_pointer_cast<cube_view>(_st_ref)->aggregation_method() = v->aggregation_method();
            std::dynamic_pointer_cast<cube_view>(_st_ref)->resampling_method() = v->resampling_method();
        }
    }

   private:
    const std::shared_ptr<image_collection> _collection;

    void load_bands();

    band_collection _input_bands;
};

#endif  //IMAGE_COLLECTION_CUBE_H
