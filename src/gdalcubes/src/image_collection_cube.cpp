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
#include "image_collection_cube.h"

#include <gdal_utils.h>

#include <map>
#include <unordered_map>

#include "error.h"
#include "utils.h"
#include "warp.h"

namespace gdalcubes {

image_collection_cube::image_collection_cube(std::shared_ptr<image_collection> ic, cube_view v) : cube(std::make_shared<cube_view>(v)), _collection(ic), _input_bands(), _mask(nullptr), _mask_band(""), _strict(true) { load_bands(); }
image_collection_cube::image_collection_cube(std::string icfile, cube_view v) : cube(std::make_shared<cube_view>(v)), _collection(std::make_shared<image_collection>(icfile)), _input_bands(), _mask(nullptr), _mask_band(""), _strict(true) { load_bands(); }
image_collection_cube::image_collection_cube(std::shared_ptr<image_collection> ic, std::string vfile) : cube(std::make_shared<cube_view>(cube_view::read_json(vfile))), _collection(ic), _input_bands(), _mask(nullptr), _mask_band(""), _strict(true) { load_bands(); }
image_collection_cube::image_collection_cube(std::string icfile, std::string vfile) : cube(std::make_shared<cube_view>(cube_view::read_json(vfile))), _collection(std::make_shared<image_collection>(icfile)), _input_bands(), _mask(nullptr), _mask_band(""), _strict(true) { load_bands(); }
image_collection_cube::image_collection_cube(std::shared_ptr<image_collection> ic) : cube(), _collection(ic), _input_bands(), _mask(nullptr), _mask_band(""), _strict(true) {
    st_reference(std::make_shared<cube_view>(image_collection_cube::default_view(_collection)));
    load_bands();
}

image_collection_cube::image_collection_cube(std::string icfile) : cube(), _collection(std::make_shared<image_collection>(icfile)), _input_bands(), _mask(nullptr), _mask_band(""), _strict(true) {
    st_reference(std::make_shared<cube_view>(image_collection_cube::default_view(_collection)));
    load_bands();
}

std::string image_collection_cube::to_string() {
    std::stringstream out;
    std::shared_ptr<cube_view> x = std::dynamic_pointer_cast<cube_view>(_st_ref);
    out << "GDAL IMAGE COLLECTION CUBE with (x,y,t)=(" << view()->nx() << "," << view()->ny() << "," << view()->nt() << ") cells in " << count_chunks() << " chunks." << std::endl;
    return out.str();
}

struct aggregation_state {
   public:
    aggregation_state(coords_nd<uint32_t, 4> size_btyx) : _size_btyx(size_btyx) {}
    virtual ~aggregation_state() {}

    virtual void init() = 0;
    virtual void update(void *chunk_buf, void *img_buf, uint32_t t) = 0;
    virtual void finalize(void *buf) = 0;

   protected:
    coords_nd<uint32_t, 4> _size_btyx;
};

struct aggregation_state_mean : public aggregation_state {
    aggregation_state_mean(coords_nd<uint32_t, 4> size_btyx) : aggregation_state(size_btyx), _m_count() {}

    ~aggregation_state_mean() {}

    void init() override {
        _m_count.resize(_size_btyx[0] * _size_btyx[1] * _size_btyx[2] * _size_btyx[3]);
    }

    void update(void *chunk_buf, void *img_buf, uint32_t t) override {
        for (uint32_t ib = 0; ib < _size_btyx[0]; ++ib) {
            uint32_t chunk_buf_offset = ib * _size_btyx[1] * _size_btyx[2] * _size_btyx[3] + t * _size_btyx[2] * _size_btyx[3];
            uint32_t img_buf_offset = ib * _size_btyx[2] * _size_btyx[3];

            // iterate over all pixels
            for (uint32_t i = 0; i < _size_btyx[2] * _size_btyx[3]; ++i) {
                if (std::isnan(((double *)img_buf)[img_buf_offset + i])) continue;
                if (std::isnan(((double *)chunk_buf)[chunk_buf_offset + i])) {
                    ((double *)chunk_buf)[chunk_buf_offset + i] = ((double *)img_buf)[img_buf_offset + i];
                    _m_count[chunk_buf_offset + i] = 1;
                } else {
                    ((double *)chunk_buf)[chunk_buf_offset + i] += ((double *)img_buf)[img_buf_offset + i];
                    _m_count[chunk_buf_offset + i] += 1;
                }
            }
        }
    }

    void finalize(void *buf) override {
        for (uint32_t i = 0; i < _size_btyx[0] * _size_btyx[1] * _size_btyx[2] * _size_btyx[3]; ++i) {
            if (!std::isnan(((double *)buf)[i])) {
                ((double *)buf)[i] /= (double)(_m_count[i]);
            }
        }
        _m_count.clear();
    }

   private:
    std::vector<uint32_t> _m_count;
};

struct aggregation_state_median : public aggregation_state {
    aggregation_state_median(coords_nd<uint32_t, 4> size_btyx) : aggregation_state(size_btyx) {}

    void init() override {
        _m_buckets.resize(_size_btyx[0] * _size_btyx[1] * _size_btyx[2] * _size_btyx[3]);
    }

    void update(void *chunk_buf, void *img_buf, uint32_t t) override {
        // iterate over all pixels
        for (uint32_t ib = 0; ib < _size_btyx[0]; ++ib) {
            // uint32_t chunk_buf_offset = ib * _size_btyx[1] * _size_btyx[2] * _size_btyx[3] + t * _size_btyx[2] * _size_btyx[3];
            uint32_t img_buf_offset = ib * _size_btyx[2] * _size_btyx[3];
            for (uint32_t i = 0; i < _size_btyx[2] * _size_btyx[3]; ++i) {
                if (std::isnan(((double *)img_buf)[img_buf_offset + i]))
                    continue;
                else {
                    _m_buckets[ib * _size_btyx[1] * _size_btyx[2] * _size_btyx[3] + t * _size_btyx[2] * _size_btyx[3] +
                               i]
                        .push_back(((double *)img_buf)[img_buf_offset + i]);
                }
            }
        }
    }

    void finalize(void *buf) override {
        for (uint32_t i = 0; i < _size_btyx[0] * _size_btyx[1] * _size_btyx[2] * _size_btyx[3]; ++i) {
            std::vector<double> &list = _m_buckets[i];
            std::sort(list.begin(), list.end());
            if (list.size() == 0) {
                ((double *)buf)[i] = NAN;
            } else if (list.size() % 2 == 1) {
                ((double *)buf)[i] = list[list.size() / 2];
            } else {
                ((double *)buf)[i] = (list[list.size() / 2] + list[list.size() / 2 - 1]) / ((double)2);
            }
        }
    }

   private:
    std::vector<std::vector<double>> _m_buckets;
};

struct aggregation_state_first : public aggregation_state {
    aggregation_state_first(coords_nd<uint32_t, 4> size_btyx) : aggregation_state(size_btyx) {}

    void init() override {}

    void update(void *chunk_buf, void *img_buf, uint32_t t) override {
        for (uint32_t ib = 0; ib < _size_btyx[0]; ++ib) {
            uint32_t chunk_buf_offset =
                ib * _size_btyx[1] * _size_btyx[2] * _size_btyx[3] + t * _size_btyx[2] * _size_btyx[3];
            uint32_t img_buf_offset = ib * _size_btyx[2] * _size_btyx[3];
            // iterate over all pixels
            for (uint32_t i = 0; i < _size_btyx[2] * _size_btyx[3]; ++i) {
                if (std::isnan(((double *)img_buf)[img_buf_offset + i])) continue;
                if (!std::isnan(((double *)chunk_buf)[chunk_buf_offset + i]))
                    continue;
                else {
                    ((double *)chunk_buf)[chunk_buf_offset + i] = ((double *)img_buf)[img_buf_offset + i];
                }
            }
        }
    }

    void finalize(void *buf) override {}
};

struct aggregation_state_count_values : public aggregation_state {
    aggregation_state_count_values(coords_nd<uint32_t, 4> size_btyx) : aggregation_state(size_btyx) {}

    void init() override {}

    void update(void *chunk_buf, void *img_buf, uint32_t t) override {
        for (uint32_t ib = 0; ib < _size_btyx[0]; ++ib) {
            uint32_t chunk_buf_offset =
                ib * _size_btyx[1] * _size_btyx[2] * _size_btyx[3] + t * _size_btyx[2] * _size_btyx[3];
            uint32_t img_buf_offset = ib * _size_btyx[2] * _size_btyx[3];
            // iterate over all pixels
            for (uint32_t i = 0; i < _size_btyx[2] * _size_btyx[3]; ++i) {
                if (std::isnan(((double *)chunk_buf)[chunk_buf_offset + i])) {
                    ((double *)chunk_buf)[chunk_buf_offset + i] = 0;
                }
                if (std::isnan(((double *)img_buf)[img_buf_offset + i])) continue;
                ((double *)chunk_buf)[chunk_buf_offset + i] += 1;
            }
        }
    }

    void finalize(void *buf) override {}
};

struct aggregation_state_count_images : public aggregation_state {
    aggregation_state_count_images(coords_nd<uint32_t, 4> size_btyx) : aggregation_state(size_btyx) {}

    void init() override {}

    void update(void *chunk_buf, void *img_buf, uint32_t t) override {
        for (uint32_t ib = 0; ib < _size_btyx[0]; ++ib) {
            uint32_t chunk_buf_offset =
                ib * _size_btyx[1] * _size_btyx[2] * _size_btyx[3] + t * _size_btyx[2] * _size_btyx[3];
            //uint32_t img_buf_offset = ib * _size_btyx[2] * _size_btyx[3];
            // iterate over all pixels
            for (uint32_t i = 0; i < _size_btyx[2] * _size_btyx[3]; ++i) {
                if (std::isnan(((double *)chunk_buf)[chunk_buf_offset + i])) {
                    ((double *)chunk_buf)[chunk_buf_offset + i] = 0;
                }
                ((double *)chunk_buf)[chunk_buf_offset + i] += 1;
            }
        }
    }

    void finalize(void *buf) override {}
};

struct aggregation_state_last : public aggregation_state {
    aggregation_state_last(coords_nd<uint32_t, 4> size_btyx) : aggregation_state(size_btyx) {}

    void init() override {}

    void update(void *chunk_buf, void *img_buf, uint32_t t) override {
        for (uint32_t ib = 0; ib < _size_btyx[0]; ++ib) {
            uint32_t chunk_buf_offset = ib * _size_btyx[1] * _size_btyx[2] * _size_btyx[3] + t * _size_btyx[2] * _size_btyx[3];
            uint32_t img_buf_offset = ib * _size_btyx[2] * _size_btyx[3];
            // iterate over all pixels
            for (uint32_t i = 0; i < _size_btyx[2] * _size_btyx[3]; ++i) {
                if (std::isnan(((double *)img_buf)[img_buf_offset + i])) continue;
                ((double *)chunk_buf)[chunk_buf_offset + i] = ((double *)img_buf)[img_buf_offset + i];
            }
        }
    }

    void finalize(void *buf) override {}
};

struct aggregation_state_min : public aggregation_state {
    aggregation_state_min(coords_nd<uint32_t, 4> size_btyx) : aggregation_state(size_btyx) {}

    void init() override {}

    void update(void *chunk_buf, void *img_buf, uint32_t t) override {
        for (uint32_t ib = 0; ib < _size_btyx[0]; ++ib) {
            uint32_t chunk_buf_offset = ib * _size_btyx[1] * _size_btyx[2] * _size_btyx[3] + t * _size_btyx[2] * _size_btyx[3];
            uint32_t img_buf_offset = ib * _size_btyx[2] * _size_btyx[3];
            // iterate over all pixels
            for (uint32_t i = 0; i < _size_btyx[2] * _size_btyx[3]; ++i) {
                if (std::isnan(((double *)img_buf)[img_buf_offset + i])) continue;
                if (std::isnan(((double *)chunk_buf)[chunk_buf_offset + i])) {
                    ((double *)chunk_buf)[chunk_buf_offset + i] = ((double *)img_buf)[img_buf_offset + i];
                } else {
                    ((double *)chunk_buf)[chunk_buf_offset + i] = std::min(((double *)chunk_buf)[chunk_buf_offset + i], ((double *)img_buf)[img_buf_offset + i]);
                }
            }
        }
    }

    void finalize(void *buf) override {}
};

struct aggregation_state_max : public aggregation_state {
    aggregation_state_max(coords_nd<uint32_t, 4> size_btyx) : aggregation_state(size_btyx) {}

    void init() override {}

    void update(void *chunk_buf, void *img_buf, uint32_t t) override {
        for (uint32_t ib = 0; ib < _size_btyx[0]; ++ib) {
            uint32_t chunk_buf_offset = ib * _size_btyx[1] * _size_btyx[2] * _size_btyx[3] + t * _size_btyx[2] * _size_btyx[3];
            uint32_t img_buf_offset = ib * _size_btyx[2] * _size_btyx[3];
            // iterate over all pixels
            for (uint32_t i = 0; i < _size_btyx[2] * _size_btyx[3]; ++i) {
                if (std::isnan(((double *)img_buf)[img_buf_offset + i])) continue;
                if (std::isnan(((double *)chunk_buf)[chunk_buf_offset + i])) {
                    ((double *)chunk_buf)[chunk_buf_offset + i] = ((double *)img_buf)[img_buf_offset + i];
                } else {
                    ((double *)chunk_buf)[chunk_buf_offset + i] = std::max(((double *)chunk_buf)[chunk_buf_offset + i], ((double *)img_buf)[img_buf_offset + i]);
                }
            }
        }
    }

    void finalize(void *buf) override {}
};

struct aggregation_state_none : public aggregation_state {
    aggregation_state_none(coords_nd<uint32_t, 4> size_btyx) : aggregation_state(size_btyx) {}

    void init() override {}
    void update(void *chunk_buf, void *img_buf, uint32_t t) override {
        for (uint32_t ib = 0; ib < _size_btyx[0]; ++ib) {
            uint32_t chunk_buf_offset = ib * _size_btyx[1] * _size_btyx[2] * _size_btyx[3] + t * _size_btyx[2] * _size_btyx[3];
            uint32_t img_buf_offset = ib * _size_btyx[2] * _size_btyx[3];
            memcpy(((double *)chunk_buf) + chunk_buf_offset, ((double *)img_buf) + img_buf_offset, sizeof(double) * _size_btyx[2] * _size_btyx[3]);
        }
    }
    void finalize(void *buf) override {}
};

/*
 * The procedure to read data for a chunk is the following:
 * 1. Exclude images that are completely ouside the spatiotemporal chunk boundaries
 * 2. create a temporary in-memory VRT dataset which crops images at the boundary of the corresponding chunks and selects its bands
 * 3. use gdal warp to reproject the VRT dataset to an in-memory GDAL dataset (this will take most of the time)
 * 4. use RasterIO to read from the dataset
 */
std::shared_ptr<chunk_data> image_collection_cube::read_chunk(chunkid_t id) {
    GCBS_TRACE("image_collection_cube::read_chunk(" + std::to_string(id) + ")");
    std::shared_ptr<chunk_data> out = std::make_shared<chunk_data>();
    if (id >= count_chunks()) {
        // chunk is outside of the cube, we don't need to read anything.
        GCBS_DEBUG("Chunk id " + std::to_string(id) + " is out of range");
        return out;
    }

    // Find intersecting images from collection and iterate over these
    // Note that these are ordered by image id and descriptor
    bounds_st cextent = bounds_from_chunk(id);
    std::vector<image_collection::find_range_st_row> datasets = _collection->find_range_st(cextent, _st_ref->srs(), std::vector<std::string>(), std::vector<std::string>{"gdalrefs.image_id", "gdalrefs.descriptor"});

    if (datasets.empty()) {
        //GCBS_DEBUG("Chunk " + std::to_string(id) + " does not intersect with any image from the image_collection_cube");
        return out;  // empty chunk data
    }

    // Derive how many pixels the chunk has (this varies for chunks at the boundary of the view)
    coords_nd<uint32_t, 3> size_tyx = chunk_size(id);
    coords_nd<uint32_t, 4> size_btyx = {_bands.count(), size_tyx[0], size_tyx[1], size_tyx[2]};
    out->size(size_btyx);

    if (size_btyx[0] * size_btyx[1] * size_btyx[2] * size_btyx[3] == 0)
        return out;

    // Fill buffers accordingly
    out->buf(std::calloc(size_btyx[0] * size_btyx[1] * size_btyx[2] * size_btyx[3], sizeof(double)));
    double *begin = (double *)out->buf();
    double *end = ((double *)out->buf()) + size_btyx[0] * size_btyx[1] * size_btyx[2] * size_btyx[3];
    std::fill(begin, end, NAN);

    OGRSpatialReference proj_out;
    proj_out.SetFromUserInput(_st_ref->srs().c_str());

    aggregation_state *agg = nullptr;
    if (view()->aggregation_method() == aggregation::aggregation_type::AGG_MEAN) {
        agg = new aggregation_state_mean(size_btyx);
    } else if (view()->aggregation_method() == aggregation::aggregation_type::AGG_MIN) {
        agg = new aggregation_state_min(size_btyx);
    } else if (view()->aggregation_method() == aggregation::aggregation_type::AGG_MAX) {
        agg = new aggregation_state_max(size_btyx);
    } else if (view()->aggregation_method() == aggregation::aggregation_type::AGG_FIRST) {
        agg = new aggregation_state_first(size_btyx);
    } else if (view()->aggregation_method() == aggregation::aggregation_type::AGG_LAST) {
        agg = new aggregation_state_last(size_btyx);
    } else if (view()->aggregation_method() == aggregation::aggregation_type::AGG_MEDIAN) {
        agg = new aggregation_state_median(size_btyx);
    } else if (view()->aggregation_method() == aggregation::aggregation_type::AGG_IMAGE_COUNT) {
        agg = new aggregation_state_count_images(size_btyx);
    } else if (view()->aggregation_method() == aggregation::aggregation_type::AGG_VALUE_COUNT) {
        agg = new aggregation_state_count_values(size_btyx);
    } else
        agg = new aggregation_state_none(size_btyx);

    agg->init();

    void *img_buf = std::calloc(size_btyx[0] * size_btyx[3] * size_btyx[2], sizeof(double));
    void *mask_buf = nullptr;
    if (_mask) {
        mask_buf = std::calloc(size_btyx[3] * size_btyx[2], sizeof(double));
    }

    uint32_t i = 0;
    uint32_t count_success = 0; // count successful image reads
    while (i < datasets.size()) {
        std::pair<std::string, uint16_t> mask_dataset_band;
        mask_dataset_band.first = "";
        mask_dataset_band.second = 0;

        uint32_t image_id = datasets[i].image_id;
        std::string image_name = datasets[i].image_name;
        std::string src_srs = datasets[i].srs;
        datetime dt = datetime::from_string(datasets[i].datetime);
        dt.unit(_st_ref->dt_unit());  // explicit datetime unit cast
        duration temp_dt = _st_ref->dt();
        int itime = (dt - cextent.t0) / temp_dt;  // time index, at which time slice of the chunk buffer will this image be written?

        // map: gdal dataset descriptor -> list of contained bands (name and number)
        std::unordered_map<std::string, std::vector<std::tuple<std::string, uint16_t>>> image_datasets;
        while (i < datasets.size() && datasets[i].image_id == image_id) {
            std::string descriptor_name = datasets[i].descriptor;
            while (i < datasets.size() && datasets[i].image_id == image_id && datasets[i].descriptor == descriptor_name) {
                if (_mask) {
                    if (datasets[i].band_name == _mask_band) {
                        mask_dataset_band.first = descriptor_name;
                        mask_dataset_band.second = datasets[i].band_num;
                    }
                }
                if (_bands.has(datasets[i].band_name)) {
                    image_datasets[descriptor_name].push_back(std::tuple<std::string, uint16_t>(datasets[i].band_name, datasets[i].band_num));
                }
                ++i;
            }
        }
        if (image_datasets.empty()) {
            continue;
        }
        if (itime < 0 || itime >= (int)(out->size()[1])) {
            continue;  // image would be written outside of the chunk buffer
        }

        // refill for all images
        std::fill((double *)img_buf, ((double *)img_buf) + size_btyx[0] * size_btyx[3] * size_btyx[2], NAN);

        for (auto it = image_datasets.begin(); it != image_datasets.end(); ++it) {
            GDALDataset *bandsel_vrt = nullptr;
            std::string bandsel_vrt_name = "";
            GDALDataset *g = (GDALDataset *)GDALOpen(it->first.c_str(), GA_ReadOnly);
            if (!g) {
                GCBS_WARN("GDAL could not open '" + it->first + "':  ERROR (" + std::to_string(CPLGetLastErrorNo()) + "): " + CPLGetLastErrorMsg());
                if (_strict) {
                    if (img_buf) std::free(img_buf);
                    if (mask_buf) std::free(mask_buf);
                    delete agg;
                    out = std::make_shared<chunk_data>();
                    out->set_status(chunk_data::chunk_status::ERROR);
                    return out;
                }
                GCBS_WARN("Dataset '" + it->first + "' will be ignored.");
                out->set_status(chunk_data::chunk_status::INCOMPLETE);
                continue;
            }

            if (g->GetRasterCount() == 0) {
                GCBS_DEBUG("GDAL dataset'" + it->first + "' does not contain any raster bands and will be ignored.");
                continue;
            }

            // If input dataset has more bands than requested
            bool create_band_subset_vrt = false;
            if (g->GetRasterCount() > int(it->second.size())) {
                create_band_subset_vrt = true;
                // create temporary VRT dataset

                CPLStringList translate_args;
                translate_args.AddString("-of");
                translate_args.AddString("VRT");

                for (uint16_t b = 0; b < it->second.size(); ++b) {
                    translate_args.AddString("-b");
                    translate_args.AddString(std::to_string(std::get<1>(it->second[b])).c_str());
                }

                GDALTranslateOptions *trans_options = GDALTranslateOptionsNew(translate_args.List(), NULL);
                if (trans_options == NULL) {
                    GCBS_ERROR("Cannot create gdal_translate options");
                    throw std::string("Cannot create gdal_translate options");
                }
                bandsel_vrt_name = "/vsimem/" + utils::generate_unique_filename() + ".vrt";
                bandsel_vrt = (GDALDataset *)GDALTranslate(bandsel_vrt_name.c_str(), (GDALDatasetH)g, trans_options, NULL);
                if (bandsel_vrt == NULL) {
                    create_band_subset_vrt = false;
                }
                GDALTranslateOptionsFree(trans_options);
            }

            std::vector<double> nodata_value_list;
            //std::string nodata_value_list = "";
            for (uint16_t b = 0; b < it->second.size(); ++b) {
                if (!_input_bands.get(std::get<0>(it->second[b])).no_data_value.empty()) {
                    //nodata_value_list += _input_bands.get(std::get<0>(it->second[b])).no_data_value;
                    nodata_value_list.push_back(std::stod(_input_bands.get(std::get<0>(it->second[b])).no_data_value));
                    //if (b < it->second.size() - 1) nodata_value_list += " ";
                }
            }
            if (nodata_value_list.empty()) {
                // try to derive nodata value from gdal dataset
                GDALDataset *d = (create_band_subset_vrt && bandsel_vrt != nullptr) ? bandsel_vrt : g;
                for (uint16_t b = 0; b < d->GetRasterCount(); ++b) {
                    int succ = 0;
                    double val = d->GetRasterBand(b + 1)->GetNoDataValue(&succ);
                    if (succ) {
                        nodata_value_list.push_back(val);
                    }
                }
                if ((int)nodata_value_list.size() != d->GetRasterCount()) {
                    nodata_value_list.clear();
                }
            }

            GDALDataset *gdal_out = nullptr;
            if (create_band_subset_vrt && bandsel_vrt != nullptr) {
                //gdal_out = (GDALDataset *)GDALWarp("", NULL, 1, (GDALDatasetH *)(&bandsel_vrt), warp_opts, NULL);
                gdal_out = gdalwarp_client::warp(bandsel_vrt, src_srs.c_str(), _st_ref->srs().c_str(), cextent.s.left, cextent.s.right,
                                                 cextent.s.top, cextent.s.bottom, size_btyx[3], size_btyx[2],
                                                 resampling::to_string(view()->resampling_method()), nodata_value_list);
            } else {
                //gdal_out = (GDALDataset *)GDALWarp("", NULL, 1, (GDALDatasetH *)(&g), warp_opts, NULL);
                gdal_out = gdalwarp_client::warp(g, src_srs.c_str(), _st_ref->srs().c_str(), cextent.s.left, cextent.s.right,
                                                 cextent.s.top, cextent.s.bottom, size_btyx[3], size_btyx[2],
                                                 resampling::to_string(view()->resampling_method()), nodata_value_list);
            }
            if (!gdal_out) {
                GCBS_WARN("GDAL could not warp '" + it->first + "':  ERROR (" + std::to_string(CPLGetLastErrorNo()) + "): " + CPLGetLastErrorMsg());
                if (_strict) {
                    if (img_buf) std::free(img_buf);
                    if (mask_buf) std::free(mask_buf);
                    delete agg;
                    out = std::make_shared<chunk_data>();
                    out->set_status(chunk_data::chunk_status::ERROR);
                    return out;
                }
                GCBS_WARN("Dataset '" + it->first + "' will be ignored.");
                out->set_status(chunk_data::chunk_status::INCOMPLETE);
                continue;
            }

            // For each band, call RasterIO to read and copy data to the right position in the buffers
            for (uint16_t b = 0; b < it->second.size(); ++b) {
                uint16_t b_internal = _bands.get_index(std::get<0>(it->second[b]));

                // Make sure that b_internal is valid in order to prevent buffer overflows
                if (b_internal < 0 || b_internal >= out->size()[0])
                    continue;

                CPLErr res;
                if (create_band_subset_vrt) {  // bands have been renumbered / sorted according to order of it->second
                    res = gdal_out->GetRasterBand(b + 1)->RasterIO(GF_Read, 0, 0, size_btyx[3], size_btyx[2], ((double *)img_buf) + b_internal * size_btyx[2] * size_btyx[3], size_btyx[3], size_btyx[2], GDT_Float64, 0, 0, NULL);
                } else {
                    res = gdal_out->GetRasterBand(std::get<1>(it->second[b]))->RasterIO(GF_Read, 0, 0, size_btyx[3], size_btyx[2], ((double *)img_buf) + b_internal * size_btyx[2] * size_btyx[3], size_btyx[3], size_btyx[2], GDT_Float64, 0, 0, NULL);
                }
                if (res != CE_None) {
                    GCBS_WARN("RasterIO (read) failed for '" + std::string(gdal_out->GetDescription()) + "':  ERROR (" + std::to_string(CPLGetLastErrorNo()) + "): " + CPLGetLastErrorMsg());
                    if (_strict) {
                        if (img_buf) std::free(img_buf);
                        if (mask_buf) std::free(mask_buf);
                        delete agg;
                        out = std::make_shared<chunk_data>();
                        out->set_status(chunk_data::chunk_status::ERROR);
                        return out;
                    }
                    GCBS_WARN("Dataset '" + it->first + "' will be ignored.");
                    out->set_status(chunk_data::chunk_status::INCOMPLETE);
                    continue;
                }
            }
            if (!bandsel_vrt_name.empty()) {
                filesystem::remove(bandsel_vrt_name);
            }
            GDALClose(gdal_out);
        }

        // now, we have filled img_buf with data from all available bands
        if (_mask) {
            // if we apply a mask, we again read the mask band with NN / MODE resampling
            // read mask again (with NN

            // find out, which dataset has mask band
            if (mask_dataset_band.first.empty()) {
                GCBS_WARN("Missing mask band for image '" + image_name + "', mask will be ignored");
            } else {
                GDALDataset *bandsel_vrt = nullptr;
                GDALDataset *g = (GDALDataset *)GDALOpen(mask_dataset_band.first.c_str(), GA_ReadOnly);
                if (!g) {
                    GCBS_WARN("GDAL could not open '" + mask_dataset_band.first + "':  ERROR (" + std::to_string(CPLGetLastErrorNo()) + "): " + CPLGetLastErrorMsg());
                    if (_strict) {
                        if (img_buf) std::free(img_buf);
                        if (mask_buf) std::free(mask_buf);
                        delete agg;
                        return std::make_shared<chunk_data>();
                    }
                    GCBS_WARN("Mask dataset '" + mask_dataset_band.first + "' will be ignored.");
                    continue;
                }
                else {
                    // If input dataset has more bands than requested
                    bool create_band_subset_vrt = false;
                    if (g->GetRasterCount() > 1) {
                        create_band_subset_vrt = true;
                        // create temporary VRT dataset

                        CPLStringList translate_args;
                        translate_args.AddString("-of");
                        translate_args.AddString("VRT");

                        translate_args.AddString("-b");
                        translate_args.AddString(std::to_string(mask_dataset_band.second).c_str());

                        GDALTranslateOptions *trans_options = GDALTranslateOptionsNew(translate_args.List(), NULL);
                        if (trans_options == NULL) {
                            GCBS_ERROR("Cannot create gdal_translate options");
                            throw std::string("Cannot create gdal_translate options");
                        }

                        bandsel_vrt = (GDALDataset *)GDALTranslate("", (GDALDatasetH)g, trans_options, NULL);
                        if (bandsel_vrt == NULL) {
                            create_band_subset_vrt = false;
                        }
                        GDALTranslateOptionsFree(trans_options);
                    }

                    GDALDataset *gdal_out = nullptr;
                    if (create_band_subset_vrt && bandsel_vrt != nullptr) {
                        //gdal_out = (GDALDataset *)GDALWarp("", NULL, 1, (GDALDatasetH *)(&bandsel_vrt), warp_opts, NULL);
                        gdal_out = gdalwarp_client::warp(bandsel_vrt, src_srs.c_str(), _st_ref->srs().c_str(), cextent.s.left, cextent.s.right,
                                                         cextent.s.top, cextent.s.bottom, size_btyx[3], size_btyx[2],
                                                         "near", std::vector<double>());
                    } else {
                        //gdal_out = (GDALDataset *)GDALWarp("", NULL, 1, (GDALDatasetH *)(&g), warp_opts, NULL);
                        gdal_out = gdalwarp_client::warp(g, src_srs.c_str(), _st_ref->srs().c_str(), cextent.s.left, cextent.s.right,
                                                         cextent.s.top, cextent.s.bottom, size_btyx[3], size_btyx[2],
                                                         "near", std::vector<double>());
                    }
                    if (!gdal_out) {
                        GCBS_WARN("GDAL could not warp '" + mask_dataset_band.first + "':  ERROR (" + std::to_string(CPLGetLastErrorNo()) + "): " + CPLGetLastErrorMsg());
                        if (_strict) {
                            if (img_buf) std::free(img_buf);
                            if (mask_buf) std::free(mask_buf);
                            delete agg;
                            out = std::make_shared<chunk_data>();
                            out->set_status(chunk_data::chunk_status::ERROR);
                            return out;
                        }
                        GCBS_WARN("Mask dataset '" + mask_dataset_band.first + "' will be ignored.");
                        out->set_status(chunk_data::chunk_status::INCOMPLETE);
                        continue;
                    }
                    CPLErr res = gdal_out->GetRasterBand(mask_dataset_band.second)->RasterIO(GF_Read, 0, 0, size_btyx[3], size_btyx[2], mask_buf, size_btyx[3], size_btyx[2], GDT_Float64, 0, 0, NULL);
                    
                    if (res != CE_None) {
                        GCBS_WARN("RasterIO (read) failed for '" + std::string(gdal_out->GetDescription()) + "':  ERROR (" + std::to_string(CPLGetLastErrorNo()) + "): " + CPLGetLastErrorMsg());
                        if (_strict) {
                            if (img_buf) std::free(img_buf);
                            if (mask_buf) std::free(mask_buf);
                            delete agg;
                            out = std::make_shared<chunk_data>();
                            out->set_status(chunk_data::chunk_status::ERROR);
                            return out;
                        }
                        GCBS_WARN("Mask dataset '" + mask_dataset_band.first + "' will be ignored.");
                        out->set_status(chunk_data::chunk_status::INCOMPLETE);
                        continue;
                    }
                    GDALClose(gdal_out);
                    _mask->apply((double *)mask_buf, (double *)img_buf, size_btyx[0], size_btyx[2], size_btyx[3]);
                }
            }
        }

        // feed the aggregator
        agg->update(out->buf(), img_buf, itime);

        count_success++;
    }


    if (out->status() == chunk_data::chunk_status::INCOMPLETE && count_success == 0) {
        out->set_status(chunk_data::chunk_status::ERROR);
    }

    agg->finalize(out->buf());
    delete agg;

    std::free(img_buf);
    if (mask_buf) std::free(mask_buf);

    // check if chunk is completely NAN and if yes, return empty chunk
    if (out->all_nan()) {
        auto s = out->status();
        out = std::make_shared<chunk_data>();
        out->set_status(s);
    }

    //    CPLFree(srs_out_str);
    return out;
}

void image_collection_cube::load_bands() {
    // Access image collection and fetch band information
    std::vector<image_collection::bands_row> band_info = _collection->get_available_bands();

    // this is the band information of the cube, not of the original image bands
    for (uint16_t ib = 0; ib < band_info.size(); ++ib) {
        band bout(band_info[ib].name);
        band bin(band_info[ib].name);
        bout.unit = band_info[ib].unit;
        bin.unit = band_info[ib].unit;
        bout.type = "float64";
        bin.type = utils::string_from_gdal_type(band_info[ib].type);
        bout.scale = band_info[ib].scale;
        bin.scale = band_info[ib].scale;
        bout.offset = band_info[ib].offset;
        bin.offset = band_info[ib].offset;
        bout.no_data_value = std::to_string(NAN);
        bin.no_data_value = band_info[ib].nodata;
        _bands.add(bout);
        _input_bands.add(bin);
    }
}

cube_view image_collection_cube::default_view(std::shared_ptr<image_collection> ic) {
    bounds_st extent = ic->extent();

    cube_view out;

    std::string srs = ic->distinct_srs();
    if (srs.empty()) {
        out.srs("EPSG:3857");
    } else {
        out.srs(srs);
    }

    // Transform WGS84 boundaries to target srs
    bounds_2d<double> ext_transformed = extent.s.transform("EPSG:4326", out.srs().c_str());


    uint32_t ncells_space = 512 * 512;
    double asp_ratio = (ext_transformed.right - ext_transformed.left) / (ext_transformed.top - ext_transformed.bottom);

    out.set_x_axis(ext_transformed.left, ext_transformed.right, (uint32_t)std::fmax((uint32_t)sqrt(ncells_space * asp_ratio), 1.0));
    out.set_y_axis(ext_transformed.bottom, ext_transformed.top, (uint32_t)std::fmax((uint32_t)sqrt(ncells_space * 1 / asp_ratio), 1.0));


    duration d = extent.t1- extent.t0;

    if (extent.t0 == extent.t1) {
        out.set_t_axis(extent.t0, extent.t1, duration(1, datetime_unit::DAY));
    } else {
        datetime_unit u;
        if (d.convert(datetime_unit::YEAR).dt_interval > 4) {
           u = datetime_unit::YEAR;
        } else if (d.convert(datetime_unit::MONTH).dt_interval > 4) {
            u = datetime_unit::MONTH;
        } else if (d.convert(datetime_unit::DAY).dt_interval > 4) {
            u = datetime_unit::DAY;
        } else if (d.convert(datetime_unit::HOUR).dt_interval > 4) {
            u = datetime_unit::HOUR;
        } else if (d.convert(datetime_unit::MINUTE).dt_interval > 4) {
            u = datetime_unit::MINUTE;
        } else if (d.convert(datetime_unit::SECOND).dt_interval > 4) {
            u = datetime_unit::SECOND;
        } else {
            u = datetime_unit::DAY;
        }
        datetime t0 = extent.t0;
        t0.unit(u);
        datetime t1 = extent.t1;
        t1.unit(u);
        out.set_t_axis(t0, t1, 4);
    }

    out.aggregation_method() = aggregation::aggregation_type::AGG_FIRST;
    out.resampling_method() = resampling::resampling_type::RSMPL_NEAR;

    return out;
}

void image_collection_cube::select_bands(std::vector<std::string> bands) {
    if (bands.empty()) {
        load_bands();  // restore band selection from original image collection
        return;
    }
    // Check that all given bands exist
    for (uint16_t i = 0; i < bands.size(); ++i) {
        if (!_bands.has(bands[i])) {
            GCBS_ERROR("Band '" + bands[i] + "' does not exist in image collection");
            return;
        }
    }
    band_collection sel;
    for (uint16_t i = 0; i < bands.size(); ++i) {
        band b = _bands.get(bands[i]);
        sel.add(b);
    }
    _bands = sel;
}

void image_collection_cube::select_bands(std::vector<uint16_t> bands) {
    if (bands.empty()) {
        load_bands();  // restore band selection from original image collection
        return;
    }
    // Check that all given bands exist
    for (uint16_t i = 0; i < bands.size(); ++i) {
        if (!(bands[i] >= 0 && bands[i] < _bands.count())) {
            GCBS_ERROR("Band '" + std::to_string(bands[i]) + "' does not exist in image collection");
            return;
        }
    }
    band_collection sel;
    for (uint16_t i = 0; i < bands.size(); ++i) {
        band b = _bands.get(bands[i]);
        sel.add(b);
    }
    _bands = sel;
}

}  // namespace gdalcubes
