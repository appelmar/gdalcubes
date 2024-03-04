/*
    MIT License

    Copyright (c) 2021 Marius Appel <marius.appel@hs-bochum.de>

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
#include "simple_cube.h"

#include <set>

#include "datetime.h"

namespace gdalcubes {

simple_cube::simple_cube(std::vector<std::string> files, std::vector<std::string> datetime_values, std::vector<std::string> bands,
                         std::vector<std::string> band_names, double dx, double dy) : cube(), _in_files(files), _in_datetime(datetime_values), _in_bands(bands), _in_band_names(band_names), _in_dx(dx), _in_dy(dy), _strict(true), _orig_bands(), _band_selection() {
    if (files.size() != datetime_values.size()) {
        GCBS_ERROR("Number of files differs from number of provided datetime_values values");
        throw std::string("Number of files differs from number of provided datetime_values values");
    }

    std::shared_ptr<cube_stref_labeled_time> stref = std::make_shared<cube_stref_labeled_time>();

    if (bands.empty()) {
        // All files have the same band(s)

        // TODO: open first file and extract metadata (st_reference, bands)
        GDALDataset *dataset;

        bool success = false;
        for (uint16_t i=0; i<files.size() && !success; ++i) {
            dataset = (GDALDataset *)GDALOpen(files[i].c_str(), GA_ReadOnly);
            if (!dataset) {
                GCBS_DEBUG("GDAL failed to open '" + files[i] + "'");
                continue;
            }

            double affine_in[6] = {0, 0, 1, 0, 0, 1};
            if (dataset->GetGeoTransform(affine_in) != CE_None) {
                GDALClose(dataset);
                GCBS_DEBUG("GDAL failed to fetch geotransform parameters for '" + files[i] + "'");
                continue;
            }

            // Make sure that grid axes are aligned with south/ north and west/east direction of cooordinate reference system
            if (std::abs(affine_in[2]) > std::nextafter(0.0, 1.0) ||
                std::abs(affine_in[4]) > std::nextafter(0.0, 1.0)) {
                GCBS_ERROR("Simple cubes do not support rotated rasters, please use an image collection instead'" + files[0] + "'");
                throw std::string("Simple cubes do not support rotated rasters, please use an image collection instead");
            }

            bounds_2d<double> bbox;
            bbox.left = affine_in[0];
            bbox.right = affine_in[0] + affine_in[1] * dataset->GetRasterXSize() + affine_in[2] * dataset->GetRasterYSize();
            bbox.top = affine_in[3];
            bbox.bottom = affine_in[3] + affine_in[4] * dataset->GetRasterXSize() + affine_in[5] * dataset->GetRasterYSize();

            stref->srs(dataset->GetProjectionRef());
            if (dx <= 0) {
                stref->set_x_axis(bbox.left, bbox.right, (uint32_t)dataset->GetRasterXSize());
            } else {
                stref->set_x_axis(bbox.left, bbox.right, dx);
            }
            if (dy <= 0) {
                stref->set_y_axis(bbox.bottom, bbox.top, (uint32_t)dataset->GetRasterYSize());
            } else {
                stref->set_y_axis(bbox.bottom, bbox.top, dy);
            }
            if (!band_names.empty() && ((int32_t)band_names.size() != dataset->GetRasterCount())) {
                GDALClose(dataset);
                GCBS_ERROR("Number of provided band names does not match the number of bands in files");
                throw std::string("Number of provided band names does not match the number of bands in files");
            }
            for (uint16_t ib = 0; ib < dataset->GetRasterCount(); ++ib) {
                std::string bname = band_names.empty() ? "x" + std::to_string(ib + 1) : band_names[ib];
                band b(bname);
                int success = 0;
                b.no_data_value = std::to_string(dataset->GetRasterBand(ib + 1)->GetNoDataValue(&success));  // TODO: check for success
                b.offset = dataset->GetRasterBand(ib + 1)->GetOffset(&success);
                b.scale = dataset->GetRasterBand(ib + 1)->GetScale(&success);
                b.unit = dataset->GetRasterBand(ib + 1)->GetUnitType();
                b.type = utils::string_from_gdal_type(dataset->GetRasterBand(ib + 1)->GetRasterDataType());
                _bands.add(b);
                _orig_bands.add(b);
            } 
            success = true; // -> break
        }
        
        if (!success) {
            GCBS_ERROR("GDAL failed to open any dataset, all provided filenames / URLs seem to be invalid or inaccessible");
            throw std::string("GDAL failed to open any dataset, all provided filenames / URLs seem to be invalid or inaccessible");
        }
        GDALClose(dataset);
        
        
        

        // Create ordered time labels and index
        std::set<datetime> dtset;
        for (uint32_t i = 0; i < files.size(); ++i) {
            datetime d = datetime::from_string(datetime_values[i]);
            dtset.insert(d);
            for (uint16_t ib = 0; ib < _bands.count(); ++ib) {
                _index[d][_bands.get(ib).name] = std::pair<std::string, uint16_t>{files[i], ib + 1};
            }
        }
        std::vector<datetime> time_labels;
        std::copy(dtset.begin(), dtset.end(), std::back_inserter(time_labels));
        stref->set_time_labels(time_labels);
        stref->dt_interval(1);
        stref->dt_unit(time_labels[0].unit());

    } else {
        // Files contain different bands
        if (files.size() != bands.size()) {
            GCBS_ERROR("Number of files differs from number of provided datetime_values values");
            throw std::string("Number of files differs from number of provided datetime_values values");
        }

        // get unique bands
        std::set<std::string> band_set(bands.begin(), bands.end());

        // extract band metadata
        std::set<std::string> bands_processed;
        for (uint32_t i = 0; i < bands.size(); ++i) {
            if (bands_processed.count(bands[i]) == 0) {
                GDALDataset *dataset = (GDALDataset *)GDALOpen(files[i].c_str(), GA_ReadOnly);
                if (!dataset) {
                    GCBS_DEBUG("GDAL failed to open '" + files[i] + "'");
                    continue;
                }

                double affine_in[6] = {0, 0, 1, 0, 0, 1};
                if (dataset->GetGeoTransform(affine_in) != CE_None) {
                    GDALClose(dataset);
                    GCBS_DEBUG("GDAL failed to fetch geotransform parameters for '" + files[0] + "'");
                    continue;
                }

                // Make sure that grid axes are aligned with south/ north and west/east direction of cooordinate reference system
                if (std::abs(affine_in[2]) > std::nextafter(0.0, 1.0) ||
                    std::abs(affine_in[4]) > std::nextafter(0.0, 1.0)) {
                    GCBS_ERROR("Simple cubes do not support rotated rasters, please use an image collection instead'" + files[0] + "'");
                    throw std::string("Simple cubes do not support rotated rasters, please use an image collection instead");
                }

                bounds_2d<double> bbox;
                bbox.left = affine_in[0];
                bbox.right = affine_in[0] + affine_in[1] * dataset->GetRasterXSize() + affine_in[2] * dataset->GetRasterYSize();
                bbox.top = affine_in[3];
                bbox.bottom = affine_in[3] + affine_in[4] * dataset->GetRasterXSize() + affine_in[5] * dataset->GetRasterYSize();


                stref->srs(dataset->GetProjectionRef());
                if (dx <= 0) {
                    stref->set_x_axis(bbox.left, bbox.right, (uint32_t)dataset->GetRasterXSize());
                } else {
                    stref->set_x_axis(bbox.left, bbox.right, dx);
                }
                if (dy <= 0) {
                    stref->set_y_axis(bbox.bottom, bbox.top, (uint32_t)dataset->GetRasterYSize());
                } else {
                    stref->set_y_axis(bbox.bottom, bbox.top, dy);
                }
                if (dataset->GetRasterCount() > 1) {
                    GCBS_WARN("Assuming a 1:1 relationship between files and bands but at least one file contains > 1 band that will be ignored.'");
                }

                band b(bands[i]);
                int success = 0;
                b.no_data_value = std::to_string(dataset->GetRasterBand(1)->GetNoDataValue(&success));  // TODO: check for success
                b.offset = dataset->GetRasterBand(1)->GetOffset(&success);
                b.scale = dataset->GetRasterBand(1)->GetScale(&success);
                b.unit = dataset->GetRasterBand(1)->GetUnitType();
                b.type = utils::string_from_gdal_type(dataset->GetRasterBand(1)->GetRasterDataType());
                _bands.add(b);
                _orig_bands.add(b);

                GDALClose(dataset);
                bands_processed.insert(bands[i]);
            }

            // if metadata for all unique bands have been already read, skip the rest
            if (bands_processed.size() == band_set.size()) {
                break;
            }
        }

        // Create ordered time labels and index
        std::set<datetime> dtset;
        for (uint32_t i = 0; i < files.size(); ++i) {
            datetime d = datetime::from_string(datetime_values[i]);
            dtset.insert(d);
            _index[d][bands[i]] = std::pair<std::string, uint16_t>{files[i], 1};
        }
        std::vector<datetime> time_labels;
        std::copy(dtset.begin(), dtset.end(), std::back_inserter(time_labels));
        stref->set_time_labels(time_labels);
        stref->dt_interval(1);
        stref->dt_unit(time_labels[0].unit());
    }

    _st_ref = stref;
}

std::string simple_cube::to_string() {
    std::stringstream out;
    std::shared_ptr<cube_view> x = std::dynamic_pointer_cast<cube_view>(_st_ref);
    out << "CUBE with (x,y,t)=(" << st_reference()->nx() << "," << st_reference()->ny() << "," << st_reference()->nt() << ") cells in " << count_chunks() << " chunks." << std::endl;
    return out.str();
}

std::shared_ptr<chunk_data> simple_cube::read_chunk(chunkid_t id) {
    std::shared_ptr<chunk_data> out = std::make_shared<chunk_data>();
    if (id >= count_chunks()) {
        // chunk is outside of the cube, we don't need to read anything.
        GCBS_DEBUG("Chunk id " + std::to_string(id) + " is out of range");
        return out;
    }

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

    bounds_st cextent = bounds_from_chunk(id);

    auto climits = chunk_limits(id);
    //////
    /*
     * 1. if each image contains all bands: find out which images are between start and end date,
     * 2. read images with rasterIO and copy to correct slice in buffer, read band-wise
     *
     */
    //////

    //


    uint32_t count_success = 0; // count successful image reads
    for (uint32_t it = 0; it < size_btyx[1]; ++it) {
        datetime dt = _st_ref->datetime_at_index(climits.low[0] + it);
        auto iter = _index.find(dt);
        if (iter == _index.end()) {
            continue;
        }
        for (uint16_t ib = 0; ib < _bands.count(); ++ib) {
            if (iter->second.find(_bands.get(ib).name) == iter->second.end()) {
                continue;
            }
            auto file_band = iter->second[_bands.get(ib).name];

            std::string gdal_file = file_band.first;
            uint16_t gdal_band = file_band.second;

            GDALDataset *dataset = (GDALDataset *)GDALOpen(gdal_file.c_str(), GA_ReadOnly);
            if (!dataset) {
                GCBS_WARN("GDAL could not open '" + gdal_file + "':  ERROR no " + std::to_string(CPLGetLastErrorNo()) + ":" + CPLGetLastErrorMsg());
                if (_strict) {
                    GDALClose(dataset);
                    out = std::make_shared<chunk_data>();
                    out->set_status(chunk_data::chunk_status::ERROR);
                    return out;
                }
                GCBS_WARN("Dataset '" + gdal_file + "' will be ignored.");
                out->set_status(chunk_data::chunk_status::INCOMPLETE);
                continue;
            }

            double affine_in[6] = {0, 0, 1, 0, 0, 1};
            if (dataset->GetGeoTransform(affine_in) != CE_None) {
                GDALClose(dataset);
                GCBS_DEBUG("GDAL failed to fetch geotransform parameters for '" + gdal_file + "'");
                continue;
            }

            bounds_2d<double> bbox;
            bbox.left = affine_in[0];
            bbox.right = affine_in[0] + affine_in[1] * dataset->GetRasterXSize() + affine_in[2] * dataset->GetRasterYSize();
            bbox.top = affine_in[3];
            bbox.bottom = affine_in[3] + affine_in[4] * dataset->GetRasterXSize() + affine_in[5] * dataset->GetRasterYSize();

            int32_t left = (cextent.s.left - bbox.left) / affine_in[1];
            int32_t bottom = dataset->GetRasterYSize() - (cextent.s.bottom - bbox.bottom) / std::abs(affine_in[5]);
            int32_t top = dataset->GetRasterYSize() - (cextent.s.top - bbox.bottom) / std::abs(affine_in[5]);
            int32_t right = (cextent.s.right - bbox.left) / affine_in[1];

            // If dx and dy have been user-defined, we need to handle
            // cubes with extent larger than the images
            if (_in_dx > 0.0 ) {
                if (right > dataset->GetRasterXSize()) {
                    right = dataset->GetRasterXSize();
                }
                if (left < 0) {
                    left = 0;
                }

            }
            if (_in_dy > 0.0) {
                if (bottom > dataset->GetRasterYSize() ) {
                    bottom = dataset->GetRasterYSize();
                }
                if (top < 0) {
                   top = 0;
                }
            }



            CPLErr res;
            // TODO: add resampling parameter?!

            if (left < 0 || left > dataset->GetRasterXSize() ||
                right < 0 || right > dataset->GetRasterXSize() ||
                top < 0 || top > dataset->GetRasterYSize() ||
                bottom < 0 || bottom > dataset->GetRasterYSize()) {
                GCBS_WARN("RasterIO (read) request out of bounds, skipping " + gdal_file);
            } else {
                // TODO: time index?!
                res = dataset->GetRasterBand(gdal_band)->RasterIO(GF_Read, left, top, right - left, bottom - top,
                                                                  (double *)(out->buf()) + ib * size_btyx[1] * size_btyx[2] * size_btyx[3] +
                                                                      it * size_btyx[2] * size_btyx[3],
                                                                  size_btyx[3], size_btyx[2], GDT_Float64, 0, 0, NULL);
                if (res != CE_None) {
                    GCBS_WARN("RasterIO (read) failed for '" + gdal_file + "':  ERROR no " + std::to_string(CPLGetLastErrorNo()) + ":" + CPLGetLastErrorMsg());
                    if (_strict) {
                        GDALClose(dataset);
                        out = std::make_shared<chunk_data>();
                        out->set_status(chunk_data::chunk_status::ERROR);
                        return out;
                    }
                    GCBS_WARN("Dataset '" + gdal_file + "' will be ignored.");
                    out->set_status(chunk_data::chunk_status::INCOMPLETE);
                    continue;
                }
            }
            GDALClose(dataset);
        }
        count_success++;
    }
    if (out->status() == chunk_data::chunk_status::INCOMPLETE && count_success == 0) {
        out->set_status(chunk_data::chunk_status::ERROR);
    }
    return out;
}

}  // namespace gdalcubes
