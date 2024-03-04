/*
    MIT License

    Copyright (c) 2020 Marius Appel <marius.appel@hs-bochum.de>

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
#ifndef IMAGE_COLLECTION_CUBE_H
#define IMAGE_COLLECTION_CUBE_H

#include <unordered_set>

#include "cube.h"
#include "image_collection.h"

namespace gdalcubes {

struct image_mask {
    virtual ~image_mask() {}
    virtual void apply(double *mask_buf, double *pixel_buf, uint32_t nb, uint32_t ny, uint32_t nx) = 0;
    virtual json11::Json as_json() = 0;
};

struct value_mask : public image_mask {
   public:
    value_mask(std::unordered_set<double> mask_values, bool invert = false, std::vector<uint8_t> bits = std::vector<uint8_t>()) : _mask_values(mask_values), _invert(invert), _bits(bits) {}
    void apply(double *mask_buf, double *pixel_buf, uint32_t nb, uint32_t ny, uint32_t nx) override {
        uint32_t bitmask = 0;
        if (!_bits.empty()) {
            for (uint8_t ib = 0; ib < _bits.size(); ++ib) {
                bitmask += (uint32_t)std::pow(2.0, (double)_bits[ib]);
            }
        }
        if (!_invert) {
            for (uint32_t ixy = 0; ixy < ny * nx; ++ixy) {
                if (!_bits.empty()) {
                    mask_buf[ixy] = (uint32_t)(mask_buf[ixy]) & bitmask;
                }
                if (_mask_values.count(mask_buf[ixy]) == 1) {  // if mask has value
                    // set all bands to NAN
                    for (uint32_t ib = 0; ib < nb; ++ib) {
                        pixel_buf[ib * nx * ny + ixy] = NAN;
                    }
                }
            }
        } else {
            for (uint32_t ixy = 0; ixy < ny * nx; ++ixy) {
                if (!_bits.empty()) {
                    mask_buf[ixy] = (uint32_t)(mask_buf[ixy]) & bitmask;
                }
                if (_mask_values.count(mask_buf[ixy]) == 0) {
                    // set all bands to NAN
                    for (uint32_t ib = 0; ib < nb; ++ib) {
                        pixel_buf[ib * nx * ny + ixy] = NAN;
                    }
                }
            }
        }
    }

    json11::Json as_json() override {
        json11::Json::object out;
        out["mask_type"] = "value_mask";
        out["values"] = _mask_values;
        out["invert"] = _invert;
        out["bits"] = _bits;
        return out;
    }

   private:
    std::unordered_set<double> _mask_values;
    bool _invert;
    std::vector<uint8_t> _bits;
};

struct range_mask : public image_mask {
   public:
    range_mask(double min, double max, bool invert = false, std::vector<uint8_t> bits = std::vector<uint8_t>()) : _min(min), _max(max), _invert(invert), _bits(bits) {}

    void apply(double *mask_buf, double *pixel_buf, uint32_t nb, uint32_t ny, uint32_t nx) override {
        uint32_t bitmask = 0;
        if (!_bits.empty()) {
            for (uint8_t ib = 0; ib < _bits.size(); ++ib) {
                bitmask += (uint32_t)std::pow(2.0, (double)_bits[ib]);
            }
        }
        if (!_invert) {
            for (uint32_t ixy = 0; ixy < ny * nx; ++ixy) {
                if (!_bits.empty()) {
                    mask_buf[ixy] = (uint32_t)(mask_buf[ixy]) & bitmask;
                }
                if (mask_buf[ixy] >= _min && mask_buf[ixy] <= _max) {
                    // set all bands to NAN
                    for (uint32_t ib = 0; ib < nb; ++ib) {
                        pixel_buf[ib * nx * ny + ixy] = NAN;
                    }
                }
            }
        } else {
            for (uint32_t ixy = 0; ixy < ny * nx; ++ixy) {
                if (!_bits.empty()) {
                    mask_buf[ixy] = (uint32_t)(mask_buf[ixy]) & bitmask;
                }
                if (mask_buf[ixy] < _min || mask_buf[ixy] > _max) {
                    // set all bands to NAN
                    for (uint32_t ib = 0; ib < nb; ++ib) {
                        pixel_buf[ib * nx * ny + ixy] = NAN;
                    }
                }
            }
        }
    }

    json11::Json as_json() override {
        json11::Json::object out;
        out["mask_type"] = "range_mask";
        out["min"] = _min;
        out["max"] = _max;
        out["invert"] = _invert;
        out["bits"] = _bits;
        return out;
    }

   private:
    double _min;
    double _max;
    bool _invert;
    std::vector<uint8_t> _bits;
};

// TODO: mask that applies a lambda expression / std::function on the mask band
//struct functor_mask : public image_mask {
//public:
//
//
//private:
//};

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

    void set_mask(std::string band, std::shared_ptr<image_mask> mask) {
        std::vector<image_collection::bands_row> bands = _collection->get_available_bands();
        for (uint16_t ib = 0; ib < bands.size(); ++ib) {
            if (bands[ib].name == band) {
                _mask = mask;
                _mask_band = band;
                return;
            }
        }
        GCBS_ERROR("Band '" + band + "' does not exist in image collection, image mask will not be modified.");
    }

    std::shared_ptr<chunk_data> read_chunk(chunkid_t id) override;

    // image_collection_cube allows changing chunk sizes from outside!
    // This is important for e.g. streaming.
    void set_chunk_size(uint32_t t, uint32_t y, uint32_t x) {
        _chunk_size = {t, y, x};
    }

    void set_strict(bool s) {
        _strict = s;
    }

    json11::Json make_constructible_json() override {
        if (_collection->is_temporary()) {
            throw std::string("ERROR in image_collection_cube::make_constructible_json(): image collection is temporary, please export as file using write() first.");
        }
        json11::Json::object out;
        out["cube_type"] = "image_collection";
        out["chunk_size"] = json11::Json::array({(int)_chunk_size[0], (int)_chunk_size[1], (int)_chunk_size[2]});
        std::string err;  // TODO: do something with err
        out["view"] = json11::Json::parse(std::dynamic_pointer_cast<cube_view>(_st_ref)->write_json_string(), err);
        out["file"] = _collection->get_filename();
        if (_mask) {
            out["mask"] = _mask->as_json();
            out["mask_band"] = _mask_band;
        }
        out["strict"] = _strict;
        return out;
    }

    static cube_view default_view(std::shared_ptr<image_collection> ic);

   private:
    const std::shared_ptr<image_collection> _collection;

    void load_bands();

    band_collection _input_bands;

    std::shared_ptr<image_mask> _mask;
    std::string _mask_band;

    bool _strict;
};

}  // namespace gdalcubes

#endif  //IMAGE_COLLECTION_CUBE_H
