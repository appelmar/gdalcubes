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

#include "extract_geom.h"

#include <gdal_utils.h>
#include <ogrsf_frmts.h>
#include <cstring>

    namespace gdalcubes {

extract_geom::extract_geom(std::shared_ptr<cube> in, std::string ogr_dataset,
                           std::string time_column, std::string ogr_layer) : cube(in->st_reference()->copy()),
                                                                             _in_cube(in), _in_ogr_dataset(ogr_dataset),
                                                                             _in_time_column(time_column),
                                                                             _in_ogr_layer(ogr_layer), _ogr_dataset(ogr_dataset),
                                                                             _ogr_layer(ogr_layer), _fid_column(""),
                                                                             _in_ogr_was_transformed(false), _is_point(false),
                                                                             _chunkmask_dataset() {
    _chunk_size[0] = _in_cube->chunk_size()[0];
    _chunk_size[1] = _in_cube->chunk_size()[1];
    _chunk_size[2] = _in_cube->chunk_size()[2];
    _bands.add(band("FID"));
    _bands.add(band("time"));
    for (uint16_t ib = 0; ib < in->size_bands(); ++ib) {
        band b = in->bands().get(ib);
        _bands.add(b);
    }



    if (!OGRGeometryFactory::haveGEOS()) {
        GCBS_ERROR("Missing GEOS support in GDAL installation");
        throw std::string("Missing GEOS support in GDAL installation");
    }


    // open input OGR dataset
    GDALDataset *in_ogr_dataset;
    in_ogr_dataset = (GDALDataset *)GDALOpenEx(ogr_dataset.c_str(), GDAL_OF_VECTOR | GDAL_OF_READONLY, NULL, NULL,
                                               NULL);
    if (in_ogr_dataset == NULL) {
        GCBS_ERROR("failed to open '" + ogr_dataset + "'");
        throw std::string("failed to open '" + ogr_dataset + "'");
    }

    OGRLayer *layer;
    if (in_ogr_dataset->GetLayerCount() > 1) {
        if (ogr_layer.empty()) {
            GCBS_WARN("input OGR dataset has multiple layers, using the first.");
            layer = in_ogr_dataset->GetLayer(0);
        } else {
            layer = in_ogr_dataset->GetLayerByName(ogr_layer.c_str());
        }
    } else {
        if (ogr_layer.empty()) {
            layer = in_ogr_dataset->GetLayer(0);
        } else {
            layer = in_ogr_dataset->GetLayerByName(ogr_layer.c_str());
        }
    }
    if (layer == NULL) {
        GDALClose(in_ogr_dataset);
        GCBS_ERROR("invalid OGR layer");
        throw std::string("invalid OGR layer");
    }
    _ogr_layer = layer->GetName();

    // If layer has more than one geometry column, only the first will be used.
    // Warn if there are more geometry columns
    if (layer->GetLayerDefn()->GetGeomFieldCount() > 1) {
        std::string geom_field_name = layer->GetLayerDefn()->GetGeomFieldDefn(0)->GetNameRef();
        GCBS_WARN("Found more than one geometry field for input features, using only the first ('" + geom_field_name + "'");
    }

    // check that cube and ogr dataset have same spatial reference system.
    OGRSpatialReference srs_cube = st_reference()->srs_ogr();
    srs_cube.AutoIdentifyEPSG();
    OGRSpatialReference srs_features = *(layer->GetSpatialRef());
    srs_features.AutoIdentifyEPSG();

    if (!srs_cube.IsSame(&srs_features)) {
        // Transform features to srs of cube

        CPLStringList translate_args;
        translate_args.AddString("-f");
        translate_args.AddString("GPKG");
        translate_args.AddString("-t_srs");
        translate_args.AddString(st_reference()->srs().c_str());
        GDALVectorTranslateOptions *opts = GDALVectorTranslateOptionsNew(translate_args.List(), NULL);
        if (opts == NULL) {
            GDALClose(in_ogr_dataset);
            GDALVectorTranslateOptionsFree(opts);
            throw std::string("ERROR in extract_geom::extract_geom(): cannot create ogr2ogr options.");
        }
        // remember filename of temporary copy with transformed input features
        _ogr_dataset = filesystem::join(filesystem::get_tempdir(), utils::generate_unique_filename() + ".gpkg");
        GDALDatasetH temp = GDALVectorTranslate(_ogr_dataset.c_str(), NULL, 1, (GDALDatasetH *)&in_ogr_dataset, opts, NULL);
        if (!temp) {
            GDALClose(in_ogr_dataset);
            GCBS_ERROR("Failed to transform input features to data cube SRS");
            throw std::string("Failed to transform input features to data cube SRS");
        }
        GDALClose(in_ogr_dataset);
        in_ogr_dataset = (GDALDataset *)temp;
        layer = in_ogr_dataset->GetLayerByName(_ogr_layer.c_str());
        GDALVectorTranslateOptionsFree(opts);
        _in_ogr_was_transformed = true;
    }

    /// Check assumption: Input dataset has FIDs
    _fid_column = layer->GetFIDColumn();
    if (_fid_column.empty()) {
        GDALClose(in_ogr_dataset);
        GCBS_ERROR("Input feature dataset must have FIDs");
        throw std::string("Input feature dataset must have FIDs");
    }

    if (!_in_time_column.empty()) {
        if (layer->GetLayerDefn()->GetFieldIndex(_in_time_column.c_str()) == -1){
            GDALClose(in_ogr_dataset);
            GCBS_ERROR("Provided time column does not exist in '" + _in_ogr_dataset + "'");
            throw std::string("Provided time column does not exist in '" + _in_ogr_dataset + "'");
        }
    }

    if (!layer->TestCapability(OLCRandomRead)) {
        GCBS_WARN("Input feature layer does not support efficient random reads; computations may take considerably longer.");
    }

    OGRwkbGeometryType geom_type = layer->GetGeomType();
    if (geom_type == wkbPoint || geom_type == wkbPoint25D || geom_type == wkbPointM ||
        geom_type == wkbPointZM) {
        _is_point = true;
    }
    else {
        _is_point = false;
    }




    CPLStringList rasterize_args;
    rasterize_args.AddString("-burn");
    rasterize_args.AddString("1");
    rasterize_args.AddString("-ot");
    rasterize_args.AddString("Byte");
    rasterize_args.AddString("-of");
    rasterize_args.AddString("GTiff");
    rasterize_args.AddString("-init");
    rasterize_args.AddString("0");
    rasterize_args.AddString("-at");
    rasterize_args.AddString("-tr");
    rasterize_args.AddString(utils::dbl_to_string(st_reference()->dx() * chunk_size()[2]).c_str());
    rasterize_args.AddString(utils::dbl_to_string(st_reference()->dy() * chunk_size()[1]).c_str());

    rasterize_args.AddString("-te");
    rasterize_args.AddString(utils::dbl_to_string(st_reference()->left()).c_str());  // xmin
    rasterize_args.AddString(utils::dbl_to_string(st_reference()->top() - chunk_size()[1] * st_reference()->dy() * count_chunks_y()).c_str());  // ymin
    rasterize_args.AddString(utils::dbl_to_string(st_reference()->left() + chunk_size()[2] * st_reference()->dx() * count_chunks_x()).c_str());  // xmax
    rasterize_args.AddString(utils::dbl_to_string(st_reference()->top()).c_str()); // ymax
    rasterize_args.AddString("-l");
    rasterize_args.AddString(_ogr_layer.c_str());

    // TODO: avoid recreating this for every worker process
    if (count_chunks_x() * count_chunks_y() < 1e6) {
        _chunkmask_dataset = "/vsimem/" + utils::generate_unique_filename(8, "chunkmask_", ".tif");
    }
    else {
        _chunkmask_dataset = filesystem::join(filesystem::get_tempdir(),utils::generate_unique_filename(8, "chunkmask_", ".tif"));
    }
    // log gdal_rasterize call
    //         std::stringstream ss;
    //         ss << "Running gdal_rasterize ";
    //         for (uint16_t iws = 0; iws < rasterize_args.size(); ++iws) {
    //             ss << rasterize_args[iws] << " ";
    //         }
    //         ss << _ogr_dataset;
    //         GCBS_DEBUG(ss.str());
    GDALRasterizeOptions *rasterize_opts = GDALRasterizeOptionsNew(rasterize_args.List(), NULL);
    if (rasterize_opts == NULL) {
        GDALRasterizeOptionsFree(rasterize_opts);
        GDALClose(in_ogr_dataset);
        throw std::string("ERROR in extract_geom::extract_geom(): cannot create gdal_rasterize options.");
    }
    int err = 0;
    GDALDataset *gdal_rasterized = (GDALDataset *)GDALRasterize(_chunkmask_dataset.c_str(), NULL, (GDALDatasetH)in_ogr_dataset, rasterize_opts, &err);
    if (gdal_rasterized == NULL) {
        GDALRasterizeOptionsFree(rasterize_opts);
        GDALClose(in_ogr_dataset);
        GCBS_ERROR("gdal_rasterize failed (error code " + std::to_string(err) + ")");
        throw std::string("gdal_rasterize failed");
    }
    GDALRasterizeOptionsFree(rasterize_opts);
    GDALClose(gdal_rasterized);
    GDALClose(in_ogr_dataset);
}

std::shared_ptr<chunk_data> extract_geom::read_chunk(chunkid_t id) {
    GCBS_TRACE("extract_geom::read_chunk(" + std::to_string(id) + ")");
    std::shared_ptr<chunk_data> out = std::make_shared<chunk_data>();

    if (id >= count_chunks()) {
        return std::make_shared<chunk_data>();  // chunk is outside of the view, we don't need to read anything.
    }


    GDALDataset *gdal_rasterized_chunkmask = (GDALDataset *)GDALOpen(_chunkmask_dataset.c_str(), GA_ReadOnly);
    if (!gdal_rasterized_chunkmask) {
        GCBS_DEBUG("GDAL cannot open '" + _chunkmask_dataset + "', chunk filtering will be ignored, which can lead to significantly longer computation times in some cases");
    }
    else {
        uint8_t chunk_has_data = 1;
        auto ccoords = chunk_coords_from_id(id);

        if (gdal_rasterized_chunkmask->GetRasterBand(1)->RasterIO(GF_Read, ccoords[2], ccoords[1], 1, 1, &chunk_has_data, 1, 1, GDT_Byte, 0, 0, NULL) == CE_None) {
            if (chunk_has_data == 0) {
                GDALClose(gdal_rasterized_chunkmask);
                return std::make_shared<chunk_data>();
            }
        }
        else {
            GCBS_DEBUG("GDAL failed to read from '" + _chunkmask_dataset + "', chunk filtering will be ignored, which can lead to significantly longer computation times in some cases");
        }
        GDALClose(gdal_rasterized_chunkmask);
    }

    bounds_st cbounds = bounds_from_chunk(id);
    auto csize = chunk_size(id);
    auto ccoords = chunk_coords_from_id(id);
    // gdal_rasterize -a "feat_id" -sql "SELECT FID as feat_id,* FROM nrw" -tr 500 500 nrw.gpkg test.tif

    // 1. Rasterize feature dataset with regard to current chunk, using FID as burn value
    GDALDataset *in_ogr_dataset;
    in_ogr_dataset = (GDALDataset *)GDALOpenEx(_ogr_dataset.c_str(), GDAL_OF_VECTOR | GDAL_OF_READONLY,
                                               NULL, NULL,NULL);
    if (in_ogr_dataset == NULL) {
        GCBS_ERROR("failed to open '" + _ogr_dataset + "'");
        throw std::string("failed to open '" + _ogr_dataset + "'");
    }


    // compare extent to chunk extent
    OGRLayer *layer = in_ogr_dataset->GetLayerByName(_ogr_layer.c_str());
    OGREnvelope extent;
    if (layer->GetExtent(&extent, true) != OGRERR_NONE) {
        GCBS_WARN("Could not derive spatial extent of input feature layer; computations may take considerably longer.");
    }
    else {
        bool outside = extent.MaxX <= cbounds.s.left ||
                       extent.MinX >= cbounds.s.right ||
                       extent.MinY >= cbounds.s.top ||
                       extent.MaxY <= cbounds.s.bottom;
        if (outside) {
            GDALClose(in_ogr_dataset);
            return out;
        }
    }


    layer->SetSpatialFilterRect(cbounds.s.left,cbounds.s.bottom, cbounds.s.right, cbounds.s.top);
    layer->ResetReading();

    std::vector<uint32_t> fids;
    std::vector<datetime> t;
    std::vector<OGREnvelope> fbbox;
    OGRFeature *cur_feature = layer->GetNextFeature();

    while (cur_feature != NULL) {
        int32_t fid = cur_feature->GetFID();
        if (fid != OGRNullFID) {

            if (!_in_time_column.empty()) {
                // Notice that we read time from all(!) features that spatially intersect with the chunk,
                // even if they do not intersect in the time dimension.
                // A more efficient way would be to apply an additional attribute filter before, but
                // since the time unit of the cube might be different, this might be less reliable.
                datetime tt = datetime::from_YmdHMS_digits(cur_feature->GetFieldAsString(_in_time_column.c_str()));
                tt.unit(st_reference()->dt_unit());
                // if outside chunk
                if (tt < cbounds.t0 || tt > cbounds.t1) {
                    OGRFeature::DestroyFeature(cur_feature);
                    cur_feature = layer->GetNextFeature();
                    continue;
                }
                t.push_back(tt);
            }
            fids.push_back(fid);
            OGREnvelope feature_bbox;
            cur_feature->GetGeometryRef()->getEnvelope(&feature_bbox);
            fbbox.push_back(feature_bbox);
        }
        OGRFeature::DestroyFeature(cur_feature);
        cur_feature = layer->GetNextFeature();
    }

    // Initialize variable-length output
    std::vector<std::vector<double>> data_frame_out;
    // colums: FID, time, cube bands
    data_frame_out.resize(_in_cube->size_bands() + 2);

    std::shared_ptr<chunk_data> dat; // Input chunk

    // Iterate over all features
    bool initialized = false;
    for (uint32_t ifeature = 0; ifeature < fids.size(); ++ifeature) {

        if (!initialized) {
            // Read input chunk and return if empty
            dat = _in_cube->read_chunk(id);
            if (dat->empty()) {
                OGRFeature::DestroyFeature(cur_feature);
                GDALClose(in_ogr_dataset);
                return out;
            }
            initialized = true;
        }


        int32_t x_start = (int32_t)std::floor((fbbox[ifeature].MinX - cbounds.s.left) /st_reference()->dx());
        int32_t x_end   = (int32_t)std::ceil((fbbox[ifeature].MaxX - cbounds.s.left) / st_reference()->dx());
        int32_t y_start = (int32_t)std::floor((cbounds.s.top - fbbox[ifeature].MaxY) / st_reference()->dy());
        int32_t y_end   = (int32_t)std::ceil((cbounds.s.top - fbbox[ifeature].MinY) / st_reference()->dy());
        if (x_end == x_start) { // point exactly on cell boundary (unlikely!)
            ++x_end;
        }
        if (y_end == y_start) { // point exactly on cell boundary (unlikely!)
            ++y_end;
        }


        x_start = std::min(std::max(x_start, 0),(int32_t)csize[2] - 1);
        x_end = std::min(std::max(x_end, 0),(int32_t)csize[2]);
        y_start = std::min(std::max(y_start, 0),(int32_t)csize[1] - 1);
        y_end = std::min(std::max(y_end, 0),(int32_t)csize[1]);

        assert(y_end > y_start);
        assert(x_end > x_start);

        // rasterize
        CPLStringList rasterize_args;
        rasterize_args.AddString("-burn");
        rasterize_args.AddString("1");
        rasterize_args.AddString("-ot");
        rasterize_args.AddString("Byte");
        rasterize_args.AddString("-of");
        rasterize_args.AddString("MEM");
        rasterize_args.AddString("-init");
        rasterize_args.AddString("0");
        rasterize_args.AddString("-tr");
        rasterize_args.AddString(utils::dbl_to_string(st_reference()->dx()).c_str());
        rasterize_args.AddString(utils::dbl_to_string(st_reference()->dy()).c_str());
        //            rasterize_args.AddString("-te");
        //            rasterize_args.AddString(std::to_string(cbounds.s.left).c_str());  // xmin
        //            rasterize_args.AddString(std::to_string(cbounds.s.bottom).c_str());  // ymin
        //            rasterize_args.AddString(std::to_string(cbounds.s.right).c_str()); // xmax
        //            rasterize_args.AddString(std::to_string(cbounds.s.top).c_str());// ymax
        rasterize_args.AddString("-te");
        rasterize_args.AddString(utils::dbl_to_string((cbounds.s.left + x_start * st_reference()->dx())).c_str());  // xmin
        rasterize_args.AddString(utils::dbl_to_string((cbounds.s.top - (y_end) * st_reference()->dy())).c_str());  // ymin
        rasterize_args.AddString(utils::dbl_to_string((cbounds.s.left + (x_end) * st_reference()->dx())).c_str());  // xmax
        rasterize_args.AddString(utils::dbl_to_string((cbounds.s.top - y_start * st_reference()->dy())).c_str());  // ymax
        rasterize_args.AddString("-where");
        std::string where = _fid_column + "=" + std::to_string(fids[ifeature]);
        rasterize_args.AddString(where.c_str());
        rasterize_args.AddString("-l");
        rasterize_args.AddString(layer->GetName());

          // log gdal_rasterize call
//         std::stringstream ss;
//         ss << "Running gdal_rasterize ";
//         for (uint16_t iws = 0; iws < rasterize_args.size(); ++iws) {
//             ss << rasterize_args[iws] << " ";
//         }
//         ss << _ogr_dataset;
//         GCBS_DEBUG(ss.str());

        GDALRasterizeOptions *rasterize_opts = GDALRasterizeOptionsNew(rasterize_args.List(), NULL);
        if (rasterize_opts == NULL) {
            GDALRasterizeOptionsFree(rasterize_opts);
            GDALClose(in_ogr_dataset);
            throw std::string("ERROR in extract_geom::read_chunk(): cannot create gdal_rasterize options.");
        }

        int err = 0;
        GDALDataset *gdal_rasterized = (GDALDataset *)GDALRasterize("", NULL, (GDALDatasetH)in_ogr_dataset, rasterize_opts, &err);
        if (gdal_rasterized == NULL) {
            GDALRasterizeOptionsFree(rasterize_opts);
            GDALClose(in_ogr_dataset);
            GCBS_ERROR("gdal_rasterize failed for feature with FID " + std::to_string(fids[ifeature]) + "(error code " + std::to_string(err) + ")");
            throw std::string("gdal_rasterize failed for feature with FID " + std::to_string(fids[ifeature]));
        }
        GDALRasterizeOptionsFree(rasterize_opts);
        uint8_t *geom_mask = (uint8_t *)std::malloc(sizeof(uint8_t) * (x_end - x_start) * (y_end - y_start));
        if (gdal_rasterized->GetRasterBand(1)->RasterIO(GF_Read, 0, 0, x_end - x_start, y_end - y_start, geom_mask, x_end - x_start, y_end - y_start, GDT_Byte, 0, 0, NULL) != CE_None) {
            GDALClose(gdal_rasterized);
            GDALClose(in_ogr_dataset);
            GCBS_ERROR("RasterIO failed" ); // TODO improve error message
            throw std::string("RasterIO failed" ); // TODO improve error message
        }
        GDALClose(gdal_rasterized);

        if (_in_time_column.empty()) { // full time series
            for (int32_t iy = y_start; iy < y_end; ++iy) {
                for (int32_t ix = x_start; ix < x_end; ++ix) {
                    // if mask is 1
                    if (geom_mask[(iy - y_start) * (x_end - x_start) + ix - x_start] == 1) {
                        for (uint32_t it = 0; it < csize[0]; ++it) {
                            bool allna = true;
                            std::vector<double> vv(dat->count_bands());
                            for (uint16_t ib = 0; ib < dat->count_bands(); ++ib) {
                                double v = ((double *)dat->buf())[ib * csize[0] * csize[1] * csize[2] +
                                                                  it * csize[1] * csize[2] +
                                                                  iy * csize[2] +
                                                                  ix];
                                vv[ib] = v;
                                if (!std::isnan(v)) {
                                    allna = false;
                                }

                            }
                            if (!allna) {
                                data_frame_out[0].push_back(fids[ifeature]); // fid
                                data_frame_out[1].push_back(ccoords[0] * _chunk_size[0] + it); // t
                                for (uint16_t ib = 0; ib < dat->count_bands(); ++ib) {
                                    data_frame_out[2+ib].push_back(vv[ib]);
                                }
                            }
                        }
                    }
                }
            }
        }
        else {

            // Derive it
            uint32_t it = st_reference()->index_at_datetime(t[ifeature]);
            // Make sure, t is inside chunk
            if (it >= ccoords[0] * _chunk_size[0] && it < ccoords[0] * _chunk_size[0] + csize[0]) {
                it = it % _chunk_size[0];
                for (int32_t iy = y_start; iy < y_end; ++iy) {
                    for (int32_t ix = x_start; ix < x_end; ++ix) {
                        // if mask is 1
                        if (geom_mask[(iy - y_start) * (x_end - x_start) + ix - x_start] == 1) {
                            bool allna = true;
                            std::vector<double> vv(dat->count_bands());
                            for (uint16_t ib = 0; ib < dat->count_bands(); ++ib) {
                                double v = ((double *)dat->buf())[ib * csize[0] * csize[1] * csize[2] +
                                                                  it * csize[1] * csize[2] +
                                                                  iy * csize[2] +
                                                                  ix];
                                vv[ib] = v;
                                if (!std::isnan(v)) {
                                    allna = false;
                                }
                            }
                            if (!allna) {
                                data_frame_out[0].push_back(fids[ifeature]); // fid
                                data_frame_out[1].push_back(ccoords[0] * _chunk_size[0] + it); // t
                                for (uint16_t ib = 0; ib < dat->count_bands(); ++ib) {
                                    data_frame_out[2+ib].push_back(vv[ib]);
                                }
                            }

                        }
                    }
                }
            }
        }
        std::free(geom_mask);
        //}
    }
    GDALClose(in_ogr_dataset);

    uint32_t nrow = 0;
    nrow = data_frame_out[0].size();
    for (uint32_t i=1; i< data_frame_out.size(); ++i) {
        if (data_frame_out[i].size() != nrow) {
            GCBS_ERROR("Invalid output in extract_geom()");
            throw std::string("Invalid output in extract_geom()");
        }
    }

    // WARNING: The following approach misuses the chunk idea, because every chunk may produce
    // a buffer of different size...
    // Make sure to NEVER export result as a cube (e.g. using write_tif_collection or write_netcdf_file)
    if (nrow > 0) {
        out->size({(uint32_t)data_frame_out.size(), nrow, 1, 1});
        out->buf(std::calloc(data_frame_out.size() * nrow, sizeof(double)));

        for (uint32_t i=0; i < data_frame_out.size(); ++i) {
            std::memcpy(&((double*)out->buf())[i * nrow], data_frame_out[i].data(), sizeof(double) * nrow);
        }
    }

    // TODO: can we always output an empty chunk and write results somehwere else?!
    // This would make sure that the result is not "mis-used"
    return out;
}

}  // namespace gdalcubes
