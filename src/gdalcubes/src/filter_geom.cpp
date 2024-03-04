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

#include "filter_geom.h"

#include <gdal_utils.h>
#include <ogrsf_frmts.h>

namespace gdalcubes {

filter_geom_cube::filter_geom_cube(std::shared_ptr<cube> in, std::string wkt, std::string srs) : cube(in->st_reference()->copy()),
                                                                                                 _in_cube(in),
                                                                                                 _wkt(wkt),
                                                                                                 _srs(srs),
                                                                                                 _ogr_dataset(""),
                                                                                                 _min_chunk_x(0),
                                                                                                 _max_chunk_x(0),
                                                                                                 _min_chunk_y(0),
                                                                                                 _max_chunk_y(0) {
    _chunk_size[0] = _in_cube->chunk_size()[0];
    _chunk_size[1] = _in_cube->chunk_size()[1];
    _chunk_size[2] = _in_cube->chunk_size()[2];
    for (uint16_t ib = 0; ib < in->size_bands(); ++ib) {
        band b = in->bands().get(ib);
        _bands.add(b);
    }

    // Derive new spatial reference (extent) based on polygon

    if (srs.empty()) {
        GCBS_ERROR("Missing spatial reference for polygon");
        throw std::string("Missing spatial reference for polygon");
    }
    OGRSpatialReference srsogr;
    srsogr.SetFromUserInput(srs.c_str());

    // 1. get bounding box of polygon
    OGRGeometry *p = nullptr;
#if GDAL_VERSION_MAJOR < 2 || (GDAL_VERSION_MAJOR == 2 && GDAL_VERSION_MINOR < 3)
    char *wktstr = (char *)std::malloc(sizeof(char) * (_wkt.length() + 1));
    _wkt.copy(wktstr, _wkt.length());
    _wkt[_wkt.length()] = '\0';
    OGRErr res = OGRGeometryFactory::createFromWkt(&wktstr, &srsogr, &p);
    //std::free(wktstr); // It seems that GDAL 2.2.3 moves the wktstr pointer, which causes segmentation faults
#else
    OGRErr res = OGRGeometryFactory::createFromWkt(_wkt.c_str(), &srsogr, &p);
#endif
    if (res != OGRERR_NONE) {
        //            if (res == OGRERR_NOT_ENOUGH_DATA) {}
        //            if (res == OGRERR_UNSUPPORTED_GEOMETRY_TYPE) {}
        //            if (res == OGRERR_CORRUPT_DATA) {}
        GCBS_ERROR("Geometry creation from WKT failed");
        throw std::string("Geometry creation from WKT failed");
    }

    // This seems to cause trouble with old GDAL versions on Travis
    //    if (p->getSpatialReference()->IsEmpty()) {
    //        GCBS_ERROR("Missing spatial reference for polygon");
    //        throw std::string("Missing spatial reference for polygon");
    //    }
    OGRwkbGeometryType geom_type = p->getGeometryType();
    if (geom_type != wkbPolygon && geom_type != wkbPolygonM && geom_type != wkbPolygon25D &&
        geom_type != wkbPolygonZM && geom_type != wkbMultiPolygon) {
        GCBS_ERROR("Input geometry is not a (multi)polygon");
        throw std::string("Input geometry is not a (multi)polygon");
    }

    // 2. If needed transform to data cube reference system
    OGRSpatialReference srs_in = _in_cube->st_reference()->srs_ogr();
    if (!srs_in.IsSame(p->getSpatialReference())) {
        p->transformTo(&srs_in);  // TODO: Error checking if transformation fails
    }

    // 3. Derive new extent, and store x and y offsets
    OGREnvelope poly_extent;
    p->getEnvelope(&poly_extent);

    int32_t iminx = std::floor(poly_extent.MinX - _in_cube->st_reference()->left()) / _in_cube->st_reference()->dx();
    int32_t imaxx = std::floor(poly_extent.MaxX - _in_cube->st_reference()->left()) / _in_cube->st_reference()->dx();

    int32_t iminy = std::floor(_in_cube->st_reference()->top() - poly_extent.MaxY) / _in_cube->st_reference()->dy();
    int32_t imaxy = std::floor(_in_cube->st_reference()->top() - poly_extent.MinY) / _in_cube->st_reference()->dy();

    bool within =
        iminx >= 0 && uint32_t(iminx) < _in_cube->st_reference()->nx() &&
        imaxx >= 0 && uint32_t(imaxx) < _in_cube->st_reference()->nx() &&
        iminy >= 0 && uint32_t(iminy) < _in_cube->st_reference()->ny() &&
        imaxy >= 0 && uint32_t(imaxy) < _in_cube->st_reference()->ny();

    if (!within) {
        GCBS_ERROR("Polygon must be located completely within the data cube");
        throw std::string("Polygon must be located completely within the data cube");
    }

    _min_chunk_x = iminx / _in_cube->chunk_size()[2];
    _max_chunk_x = imaxx / _in_cube->chunk_size()[2];
    _min_chunk_y = iminy / _in_cube->chunk_size()[1];
    _max_chunk_y = imaxy / _in_cube->chunk_size()[1];

    // Create new OGR dataset with single feature...
    //std::string output_file = filesystem::join(filesystem::get_tempdir(), utils::generate_unique_filename(8, "crop_", ".gpkg"));
    // TODO: avoid recreating this for every worker process
    std::string output_file = "/vsimem/" + utils::generate_unique_filename(8, "crop_", ".gpkg");
    _ogr_dataset = output_file;
    GDALDriver *gpkg_driver = GetGDALDriverManager()->GetDriverByName("GPKG");
    GDALDataset *gpkg_out = gpkg_driver->Create(output_file.c_str(), 0, 0, 0, GDT_Unknown, NULL);
    if (gpkg_out == NULL) {
        GCBS_ERROR("Creation of GPKG file '" + output_file + "' failed");
        throw std::string("Creation of GPKG file '" + output_file + "' failed");
    }
    OGRLayer *layer_out = gpkg_out->CreateLayer("temp", &srs_in, geom_type, NULL);
    if (layer_out == NULL) {
        GCBS_ERROR("Failed to create output layer in '" + output_file + "'");
        throw std::string("Failed to create output layer in  '" + output_file + "'");
    }

    OGRFeature *geom_feature_out = OGRFeature::CreateFeature(layer_out->GetLayerDefn());
    geom_feature_out->SetGeometry(p);
    geom_feature_out->SetFID(1);
    if (layer_out->CreateFeature(geom_feature_out) != OGRERR_NONE) {
        GCBS_ERROR("Failed to create output feature with FID '" + std::to_string(1) + "' in  '" + output_file + "'");
        throw std::string("Failed to create output feature with FID '" + std::to_string(1) + "' in  '" + output_file + "'");
    }
    OGRFeature::DestroyFeature(geom_feature_out);
    GDALClose(gpkg_out);
    OGRGeometryFactory::destroyGeometry(p);
}

std::shared_ptr<chunk_data> filter_geom_cube::read_chunk(chunkid_t id) {
    GCBS_TRACE("filter_geom_cube::read_chunk(" + std::to_string(id) + ")");

    std::shared_ptr<chunk_data> out = std::make_shared<chunk_data>();

    if (id >= count_chunks())
        return out;  // chunk is outside of the view, we don't need to read anything.

    auto chunkcoords = chunk_coords_from_id(id);
    if (chunkcoords[2] < _min_chunk_x || chunkcoords[2] > _max_chunk_x || chunkcoords[1] < _min_chunk_y || chunkcoords[1] > _max_chunk_y) {
        return out;  // chunk does not intersect with polygon
    }



    bounds_st chunk_bounds = bounds_from_chunk(id);

    GDALDataset *in_ogr_dataset;
    in_ogr_dataset = (GDALDataset *)GDALOpenEx(_ogr_dataset.c_str(), GDAL_OF_VECTOR | GDAL_OF_READONLY, NULL, NULL,
                                               NULL);
    if (in_ogr_dataset == NULL) {
        GCBS_ERROR("failed to open '" + _ogr_dataset + "'");
        throw std::string("failed to open '" + _ogr_dataset + "'");
    }

    OGRLayer *layer;
    layer = in_ogr_dataset->GetLayer(0);
    if (layer == NULL) {
        GDALClose(in_ogr_dataset);
        GCBS_ERROR("invalid OGR layer");
        throw std::string("invalid OGR layer");
    }
    OGRPolygon pp;
    OGRLinearRing a;
    bounds_2d<double> sextent = bounds_from_chunk(id).s;
    OGRSpatialReference srs_cube = st_reference()->srs_ogr();

    a.addPoint(sextent.left, sextent.bottom);
    a.addPoint(sextent.right, sextent.bottom);
    a.addPoint(sextent.right, sextent.top);
    a.addPoint(sextent.left, sextent.top);
    a.addPoint(sextent.left, sextent.bottom);
    pp.addRing(&a);
    pp.assignSpatialReference(&srs_cube);

    // iterate over all features
    bool chunk_within_polygon = false;
    bool outside = false;
    layer->ResetReading();
    OGRFeature *cur_feature = layer->GetNextFeature();  // assumption, there is only one feature
    OGRGeometry *geom = cur_feature->GetGeometryRef();
    if (geom != NULL) {
        if (geom->Contains(&pp)) {
            chunk_within_polygon = true;
        }
        else if (!geom->Intersects(&pp)) {
            outside = true;
        }
    }
    OGRFeature::DestroyFeature(cur_feature);
    layer->ResetReading();

    if (outside) {
        GDALClose(in_ogr_dataset);
        return out;
    }

    std::shared_ptr<chunk_data> in = _in_cube->read_chunk(id);
    // propagate chunk status
    if (in->status() == chunk_data::chunk_status::ERROR) {
        out->set_status(chunk_data::chunk_status::ERROR);
    }
    else if (in->status() == chunk_data::chunk_status::INCOMPLETE && out->status() != chunk_data::chunk_status::ERROR) {
        out->set_status(chunk_data::chunk_status::INCOMPLETE);
    }

    if (in->empty()) {
        GDALClose(in_ogr_dataset);
        return out;
    }

    if (chunk_within_polygon) {
        //  special case: chunk is completely within polygon -> do not rasterize but simply copy buffers
        //std::memcpy(out->buf(), in->buf(), _bands.count() * in->size()[1] * in->size()[2] * in->size()[3] * sizeof(double));
        out = in;
    } else {
        out->size({_bands.count(), in->size()[1], in->size()[2], in->size()[3]});
        out->buf(std::calloc(_bands.count() * in->size()[1] * in->size()[2] * in->size()[3], sizeof(double)));
        double *begin = (double *)out->buf();
        double *end = ((double *)out->buf()) + _bands.count() * in->size()[1] * in->size()[2] * in->size()[3];
        std::fill(begin, end, NAN);

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
        rasterize_args.AddString(std::to_string(_in_cube->st_reference()->dx()).c_str());
        rasterize_args.AddString(std::to_string(_in_cube->st_reference()->dy()).c_str());
        rasterize_args.AddString("-te");
        rasterize_args.AddString(std::to_string(chunk_bounds.s.left).c_str());    // xmin
        rasterize_args.AddString(std::to_string(chunk_bounds.s.bottom).c_str());  // ymin
        rasterize_args.AddString(std::to_string(chunk_bounds.s.right).c_str());   // xmax
        rasterize_args.AddString(std::to_string(chunk_bounds.s.top).c_str());     // ymax
        rasterize_args.AddString("-l");
        rasterize_args.AddString("temp");

        //        //log gdal_rasterize call
        //        std::stringstream ss;
        //        ss << "Running gdal_rasterize ";
        //        for (uint16_t iws = 0; iws < rasterize_args.size(); ++iws) {
        //            ss << rasterize_args[iws] << " ";
        //        }
        //        ss << _ogr_dataset;
        //        GCBS_DEBUG(ss.str());

        GDALRasterizeOptions *rasterize_opts = GDALRasterizeOptionsNew(rasterize_args.List(), NULL);
        if (rasterize_opts == NULL) {
            GDALClose(in_ogr_dataset);
            GDALRasterizeOptionsFree(rasterize_opts);
            throw std::string("ERROR in filter_geom_cub::read_chunk(): cannot create gdal_rasterize options.");
        }
        int err = 0;
        GDALDataset *gdal_rasterized = (GDALDataset *)GDALRasterize("", NULL, (GDALDatasetH)in_ogr_dataset, rasterize_opts, &err);
        if (gdal_rasterized == NULL) {
            GCBS_ERROR("gdal_rasterize failed ");
        }

        GDALRasterizeOptionsFree(rasterize_opts);

        uint8_t *geom_mask = (uint8_t *)std::malloc(sizeof(uint8_t) * in->size()[3] * in->size()[2]);
        if (gdal_rasterized->GetRasterBand(1)->RasterIO(GF_Read, 0, 0, in->size()[3], in->size()[2], geom_mask, in->size()[3], in->size()[2], GDT_Byte, 0, 0, NULL) != CE_None) {
            GCBS_ERROR("RasterIO failed");
        }
        for (uint32_t iy = 0; iy < in->size()[2]; ++iy) {
            for (uint32_t ix = 0; ix < in->size()[3]; ++ix) {
                if (geom_mask[iy * in->size()[3] + ix] != 0) {
                    for (uint32_t ib = 0; ib < in->size()[0]; ++ib) {
                        for (uint32_t it = 0; it < in->size()[1]; ++it) {
                            uint32_t idx = ib * in->size()[1] * in->size()[2] * in->size()[3] + it * in->size()[2] * in->size()[3] + iy * in->size()[3] + ix;
                            ((double *)out->buf())[idx] = ((double *)in->buf())[idx];
                        }
                    }
                }
            }
        }
        std::free(geom_mask);
        GDALClose(gdal_rasterized);
    }
    GDALClose(in_ogr_dataset);
    return out;
}

}  // namespace gdalcubes
