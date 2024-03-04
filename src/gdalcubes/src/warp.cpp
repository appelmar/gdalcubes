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

#include "warp.h"

#include <gdalwarper.h>

#include "config.h"

namespace gdalcubes {

gdalwarp_client::gdalcubes_reprojection_info *gdalwarp_client::gdal_transformation_cache::get(std::string srs_in_str, std::string srs_out_str) {
    auto q = std::pair<std::string, std::string>(srs_in_str, srs_out_str);

    _mutex.lock();
    auto x = _cache.find(q);
    if (x != _cache.end()) {
        _mutex.unlock();
        return x->second;
    }

    gdalwarp_client::gdalcubes_reprojection_info *pt = create_reprojection(srs_in_str, srs_out_str);
    _cache.insert(std::make_pair(q, pt));

    _mutex.unlock();
    return pt;
}

gdalwarp_client::gdal_transformation_cache::~gdal_transformation_cache() {
    //GCBS_INFO("CACHE HAS " + std::to_string(_cache.size()) + " reprojections");
    for (auto it = _cache.begin(); it != _cache.end(); ++it) {
        destroy_reprojection(it->second);
    }
}


GDALDataset *gdalwarp_client::warp(GDALDataset *in, std::string s_srs, std::string t_srs, double te_left,
                                   double te_right, double te_top, double te_bottom, uint32_t ts_x, uint32_t ts_y,
                                   std::string resampling, std::vector<double> srcnodata) {
    double temp[6];                                
    if (in->GetGeoTransform(temp) == CE_None) {
        return warp_simple(in, s_srs, t_srs, te_left, te_right, te_top, te_bottom, ts_x, ts_y, resampling, srcnodata);
    }
    else  { // GCPs, RCPs, or GeolocationArrays
        return warp_complex(in, s_srs, t_srs, te_left, te_right, te_top, te_bottom, ts_x, ts_y, resampling, srcnodata);
    }
}




GDALDataset *gdalwarp_client::warp_simple(GDALDataset *in, std::string s_srs, std::string t_srs, double te_left,
                                   double te_right, double te_top, double te_bottom, uint32_t ts_x, uint32_t ts_y,
                                   std::string resampling, std::vector<double> srcnodata) {
    char *wkt_out = NULL;

    OGRSpatialReference srs_out;
    srs_out.SetFromUserInput(t_srs.c_str());
    srs_out.exportToWkt(&wkt_out);

    double dst_geotransform[6];
    dst_geotransform[0] = te_left;
    dst_geotransform[1] = (te_right - te_left) / double(ts_x);
    dst_geotransform[2] = 0.0;
    dst_geotransform[3] = te_top;
    dst_geotransform[4] = 0.0;
    dst_geotransform[5] = (te_bottom - te_top) / double(ts_y);

    GDALDriver *mem_driver = (GDALDriver *)GDALGetDriverByName("MEM");
    if (mem_driver == NULL) {
        GCBS_ERROR("Cannot find GDAL MEM driver");
        throw std::string("Cannot find GDAL MEM driver");
    }

    GDALDataset *out = mem_driver->Create("", ts_x, ts_y, in->GetRasterCount(), GDT_Float64, NULL);

    out->SetProjection(wkt_out);
    out->SetGeoTransform(dst_geotransform);

    // Setup warp options.
    GDALWarpOptions *psWarpOptions = GDALCreateWarpOptions();
    if (s_srs.empty()) {
        // try reading from dataset
        s_srs = std::string(in->GetProjectionRef());
        if (s_srs.empty()) {
            GCBS_DEBUG("Failed to read input crs.");
        }
    }
    //    if (std::string(in->GetProjectionRef()).empty()) {
    //        if (in->SetProjection(s_srs.c_str()) != CE_None) {
    //            GCBS_DEBUG("Failed to overwrite projection for dataset.");
    //        }
    //    }
    //    std::string s = (in->GetProjectionRef());
    //    GCBS_DEBUG(s);
    psWarpOptions->hSrcDS = in;
    psWarpOptions->hDstDS = out;
    psWarpOptions->pfnProgress = GDALDummyProgress;
    psWarpOptions->pTransformerArg = create_transform(in, out, s_srs, t_srs);
    psWarpOptions->pfnTransformer = transform;

    // Derive best overview level to use
    int n_ov = in->GetRasterBand(1)->GetOverviewCount();
    if (config::instance()->get_gdal_use_overviews() && n_ov > 0) {
        double *x = (double *)std::malloc(sizeof(double) * 4);
        double *y = (double *)std::malloc(sizeof(double) * 4);
        int *succ = (int *)std::malloc(sizeof(int) * 4);
        x[0] = 0;
        x[1] = 0;
        x[2] = ts_x;
        x[3] = ts_x;
        y[0] = ts_y;
        y[1] = 0;
        y[2] = ts_y;
        y[3] = 0;
        transform(psWarpOptions->pTransformerArg, 1, 4, x, y, NULL, succ);

        double minx = std::min(std::min(x[0], x[1]), std::min(x[2], x[3]));
        double maxx = std::max(std::max(x[0], x[1]), std::max(x[2], x[3]));
        double target_ratio = (maxx - minx) / double(ts_x);

        int16_t ilevel = 0;
        while (ilevel < n_ov) {
            double ov_ratio = double(in->GetRasterBand(1)->GetXSize()) / double(in->GetRasterBand(1)->GetOverview(ilevel)->GetXSize());
            if (ov_ratio > target_ratio) {
                --ilevel;
                break;
            }
            ++ilevel;
        }
        if (ilevel >= n_ov) {
            ilevel = n_ov - 1;
        }
        if (ilevel >= 0) {
            //GCBS_TRACE("Using overview level" + std::to_string(ilevel));
            char **oo = nullptr;
            oo = CSLAddString(oo, ("OVERVIEW_LEVEL=" + std::to_string(ilevel)).c_str());

            std::string descr = in->GetDescription();
            GDALClose(in);
            in = (GDALDataset *)GDALOpenEx(descr.c_str(), GDAL_OF_RASTER | GDAL_OF_READONLY, NULL, oo, NULL);
            if (in != NULL) {
                destroy_transform((gdalwarp_client::gdalcubes_transform_info *)psWarpOptions->pTransformerArg);
                psWarpOptions->pTransformerArg = create_transform(in, out, s_srs, t_srs);
                psWarpOptions->hSrcDS = in;  // TODO: close in_ov
            } else {
                GCBS_WARN("Failed to open GDAL overview dataset for '" + descr + "', using original full resolution image.");
            }
            CSLDestroy(oo);
        }

        std::free(x);
        std::free(y);
        std::free(succ);
    }

    psWarpOptions->eResampleAlg = GDALResampleAlg::GRA_NearestNeighbour;
    if (resampling == "bilinear") {
        psWarpOptions->eResampleAlg = GDALResampleAlg::GRA_Bilinear;
    } else if (resampling == "cubic") {
        psWarpOptions->eResampleAlg = GDALResampleAlg::GRA_Cubic;
    } else if (resampling == "cubicspline") {
        psWarpOptions->eResampleAlg = GDALResampleAlg::GRA_CubicSpline;
    } else if (resampling == "lanczos") {
        psWarpOptions->eResampleAlg = GDALResampleAlg::GRA_Lanczos;
    } else if (resampling == "average") {
        psWarpOptions->eResampleAlg = GDALResampleAlg::GRA_Average;
    } else if (resampling == "mode") {
        psWarpOptions->eResampleAlg = GDALResampleAlg::GRA_Mode;
    } else if (resampling == "max") {
        psWarpOptions->eResampleAlg = GDALResampleAlg::GRA_Max;
    } else if (resampling == "min") {
        psWarpOptions->eResampleAlg = GDALResampleAlg::GRA_Min;
    } else if (resampling == "med") {
        psWarpOptions->eResampleAlg = GDALResampleAlg::GRA_Med;
    } else if (resampling == "q1") {
        psWarpOptions->eResampleAlg = GDALResampleAlg::GRA_Q1;
    } else if (resampling == "q3") {
        psWarpOptions->eResampleAlg = GDALResampleAlg::GRA_Q3;
    }

    psWarpOptions->nBandCount = in->GetRasterCount();
    psWarpOptions->panSrcBands = (int *)CPLMalloc(sizeof(int) * psWarpOptions->nBandCount);
    psWarpOptions->panDstBands = (int *)CPLMalloc(sizeof(int) * psWarpOptions->nBandCount);
    double *dst_nodata = (double *)CPLMalloc(sizeof(double) * psWarpOptions->nBandCount);

    // Due to issues on Windows with GDAL 2.2.3 if imag no data values are not set,
    // we explicitly set then to 0 in the following. This might be removed in the future.
    double *dst_nodata_img = (double *)CPLMalloc(sizeof(double) * psWarpOptions->nBandCount);
    for (uint16_t i = 0; i < psWarpOptions->nBandCount; ++i) {
        dst_nodata[i] = NAN;
        dst_nodata_img[i] = 0.0;
        psWarpOptions->panSrcBands[i] = i + 1;
        psWarpOptions->panDstBands[i] = i + 1;
    }
    double *src_nodata = nullptr;
    double *src_nodata_img = nullptr;

    if (!srcnodata.empty()) {
        src_nodata = (double *)CPLMalloc(sizeof(double) * psWarpOptions->nBandCount);
        src_nodata_img = (double *)CPLMalloc(sizeof(double) * psWarpOptions->nBandCount);
        if (srcnodata.size() == 1) {
            for (uint16_t i = 0; i < psWarpOptions->nBandCount; ++i) {
                src_nodata[i] = srcnodata[0];
                src_nodata_img[i] = 0.0;
            }
            psWarpOptions->padfSrcNoDataReal = src_nodata;
            psWarpOptions->padfSrcNoDataImag = src_nodata_img;
        } else if (srcnodata.size() == (uint16_t)psWarpOptions->nBandCount) {
            for (uint16_t i = 0; i < psWarpOptions->nBandCount; ++i) {
                src_nodata[i] = srcnodata[i];
                src_nodata_img[i] = 0.0;
            }
            psWarpOptions->padfSrcNoDataReal = src_nodata;
            psWarpOptions->padfSrcNoDataImag = src_nodata_img;
        } else {
            GCBS_DEBUG("Number of no data values does not match number of bands of source dataset, no data values will be ignored");
            std::free(src_nodata);
            std::free(src_nodata_img);
        }
    }
    psWarpOptions->padfDstNoDataReal = dst_nodata;
    psWarpOptions->padfDstNoDataImag = dst_nodata_img;

    char **wo = nullptr;
    wo = CSLAddString(wo, "INIT_DEST=nan");
    wo = CSLAddString(wo, ("NUM_THREADS=" + std::to_string(config::instance()->get_gdal_num_threads())).c_str());
    psWarpOptions->papszWarpOptions = wo;

    // Initialize and execute the warp operation.
    GDALWarpOperation oOperation;
    if (oOperation.Initialize(psWarpOptions) != CE_None) {
        GCBS_ERROR("Initialization of gdalwarp failed");
    }
    oOperation.ChunkAndWarpImage(0, 0, GDALGetRasterXSize(out), GDALGetRasterYSize(out));

    destroy_transform((gdalwarp_client::gdalcubes_transform_info *)psWarpOptions->pTransformerArg);
    GDALDestroyWarpOptions(psWarpOptions);

    CPLFree(wkt_out);

    if (in) {
        GDALClose(in);
    }
    return out;
}

/*
     * Source code of this function has been adapted from original GDAL code starting at
     * https://github.com/OSGeo/gdal/blob/0bfd1bcb38b3fe321fd15f3c485cfb91537faf0e/gdal/alg/gdaltransformer.cpp#L1355
     */
gdalwarp_client::gdalcubes_transform_info *gdalwarp_client::create_transform(GDALDataset *in, GDALDataset *out, std::string srs_in_str, std::string srs_out_str) {
    gdalcubes_transform_info *res = new gdalcubes_transform_info();
    res->pReprojectArg = nullptr;
    in->GetGeoTransform(res->adfSrcGeoTransform);
    if (!GDALInvGeoTransform(res->adfSrcGeoTransform, res->adfSrcInvGeoTransform)) {
        GCBS_ERROR("Cannot invert affine transformation of source image");
        destroy_transform(res);
        throw std::string("Cannot invert affine transformation of source image");
    }

    out->GetGeoTransform(res->adfDstGeoTransform);
    if (!GDALInvGeoTransform(res->adfDstGeoTransform, res->adfDstInvGeoTransform)) {
        GCBS_ERROR("Cannot invert affine transformation of destination image");
        destroy_transform(res);
        throw std::string("Cannot invert affine transformation of destination image");
    }

    // Set reprojection transform if needed
    OGRSpatialReference srs_in;
    OGRSpatialReference srs_out;
    srs_in.SetFromUserInput(srs_in_str.c_str());
    srs_out.SetFromUserInput(srs_out_str.c_str());
    if (!srs_in.IsSame(&srs_out)) {
        res->pReprojectArg = gdal_transformation_cache::instance()->get(srs_in_str, srs_out_str);
        res->pReproject = reproject;
    }
    return res;
}

/*
     * Source code of this function has been adapted from original GDAL code starting at
     * https://github.com/OSGeo/gdal/blob/0bfd1bcb38b3fe321fd15f3c485cfb91537faf0e/gdal/alg/gdaltransformer.cpp#L2281
     */
int gdalwarp_client::transform(void *pTransformerArg, int bDstToSrc, int nPointCount, double *x, double *y, double *z, int *panSuccess) {
    gdalwarp_client::gdalcubes_transform_info *psInfo = static_cast<gdalwarp_client::gdalcubes_transform_info *>(pTransformerArg);

    if (panSuccess) {
        for (int i = 0; i < nPointCount; i++) {
            panSuccess[i] = (x[i] != HUGE_VAL && y[i] != HUGE_VAL);
        }
    }

    // 1. convert source / destination image coordinates to georeferenced coordinates
    double *padfGeoTransform = nullptr;
    if (bDstToSrc) {
        padfGeoTransform = psInfo->adfDstGeoTransform;
    } else {
        padfGeoTransform = psInfo->adfSrcGeoTransform;
    }

    for (int i = 0; i < nPointCount; i++) {
        if (panSuccess) {
            if (!panSuccess[i])
                continue;
        }

        const double dfNewX = padfGeoTransform[0] + x[i] * padfGeoTransform[1] + y[i] * padfGeoTransform[2];
        const double dfNewY = padfGeoTransform[3] + x[i] * padfGeoTransform[4] + y[i] * padfGeoTransform[5];

        x[i] = dfNewX;
        y[i] = dfNewY;
    }

    // 2. Reproject to destination / source CRS if needed
    if (psInfo->pReprojectArg) {
        if (!psInfo->pReproject(psInfo->pReprojectArg, bDstToSrc,
                                nPointCount, x, y, z,
                                panSuccess))
            return 0;
    }

    // 3. convert georeferenced coordinates to destination / source image coordinates
    if (bDstToSrc) {
        padfGeoTransform = psInfo->adfSrcInvGeoTransform;
    } else {
        padfGeoTransform = psInfo->adfDstInvGeoTransform;
    }

    for (int i = 0; i < nPointCount; i++) {
        if (panSuccess) {
            if (!panSuccess[i])
                continue;
        }

        const double dfNewX = padfGeoTransform[0] + x[i] * padfGeoTransform[1] + y[i] * padfGeoTransform[2];
        const double dfNewY = padfGeoTransform[3] + x[i] * padfGeoTransform[4] + y[i] * padfGeoTransform[5];
        x[i] = dfNewX;
        y[i] = dfNewY;
    }
    return 1;
}

void gdalwarp_client::destroy_transform(gdalcubes_transform_info *transform) {
    // IMPORTANT: do NOT destroy reproject transform and args
    if (transform) {
        delete transform;
    }
}

gdalwarp_client::gdalcubes_reprojection_info *gdalwarp_client::create_reprojection(std::string srs_in_str, std::string srs_out_str) {
    // TODO: add area of interest for GDAL >= 3.0
    // TODO: add further options, e.g. from global config options

    OGRSpatialReference srs_in;
    OGRSpatialReference srs_out;

    srs_in.SetFromUserInput(srs_in_str.c_str());
    srs_out.SetFromUserInput(srs_out_str.c_str());
    OGRCoordinateTransformation *poForwardTransform = OGRCreateCoordinateTransformation(&srs_in, &srs_out);
    if (poForwardTransform == nullptr)
        return nullptr;

    gdalcubes_reprojection_info *res = new gdalcubes_reprojection_info();
    res->poForwardTransform = poForwardTransform;
    res->poReverseTransform = nullptr;
    res->poReverseTransform = OGRCreateCoordinateTransformation(&srs_out, &srs_in);
    return (res);
}

/*
    * Source code of this function has been adapted from original GDAL code starting at
    * https://github.com/OSGeo/gdal/blob/0bfd1bcb38b3fe321fd15f3c485cfb91537faf0e/gdal/alg/gdaltransformer.cpp#L2912
    */
int gdalwarp_client::reproject(void *pTransformerArg, int bDstToSrc, int nPointCount, double *x, double *y, double *z, int *panSuccess) {
    gdalwarp_client::gdalcubes_reprojection_info *psInfo = static_cast<gdalwarp_client::gdalcubes_reprojection_info *>(pTransformerArg);
    int bSuccess;
    if (bDstToSrc) {
        if (psInfo->poReverseTransform == nullptr) {
            GCBS_ERROR("Inverse coordinate transformation cannot be instantiated");
            if (panSuccess) {
                for (int i = 0; i < nPointCount; i++)
                    panSuccess[i] = FALSE;
            }
            bSuccess = false;
        } else {
#if GDAL_VERSION_MAJOR >= 3
            bSuccess = psInfo->poReverseTransform->Transform(nPointCount, x, y, z, panSuccess);
#else
            bSuccess = psInfo->poReverseTransform->Transform(nPointCount, x, y, z);
#endif
        }
    } else {
#if GDAL_VERSION_MAJOR >= 3
        bSuccess = psInfo->poForwardTransform->Transform(nPointCount, x, y, z, panSuccess);
#else
        bSuccess = psInfo->poForwardTransform->Transform(nPointCount, x, y, z);
#endif
    }
    return bSuccess;
}

void gdalwarp_client::destroy_reprojection(gdalcubes_reprojection_info *reprojection) {
    if (reprojection) {
        if (reprojection->poForwardTransform) {
            delete reprojection->poForwardTransform;
        }
        if (reprojection->poReverseTransform) {
            delete reprojection->poReverseTransform;
        }
        delete reprojection;
    }
}

GDALDataset *gdalwarp_client::warp_complex(GDALDataset *in, std::string s_srs, std::string t_srs, double te_left,
                                   double te_right, double te_top, double te_bottom, uint32_t ts_x, uint32_t ts_y,
                                   std::string resampling, std::vector<double> srcnodata) {
    char *wkt_out = NULL;

    OGRSpatialReference srs_out;
    srs_out.SetFromUserInput(t_srs.c_str());
    srs_out.exportToWkt(&wkt_out);

    double dst_geotransform[6];
    dst_geotransform[0] = te_left;
    dst_geotransform[1] = (te_right - te_left) / double(ts_x);
    dst_geotransform[2] = 0.0;
    dst_geotransform[3] = te_top;
    dst_geotransform[4] = 0.0;
    dst_geotransform[5] = (te_bottom - te_top) / double(ts_y);

    GDALDriver *mem_driver = (GDALDriver *)GDALGetDriverByName("MEM");
    if (mem_driver == NULL) {
        GCBS_ERROR("Cannot find GDAL MEM driver");
        throw std::string("Cannot find GDAL MEM driver");
    }

    GDALDataset *out = mem_driver->Create("", ts_x, ts_y, in->GetRasterCount(), GDT_Float64, NULL);

    for (uint16_t i=0; i<in->GetRasterCount(); ++i) {
        out->GetRasterBand(i+1)->Fill(NAN); // Avoid spurious output when warping fails.
    }

    out->SetProjection(wkt_out);
    out->SetGeoTransform(dst_geotransform);


    // Setup warp options.
    GDALWarpOptions *psWarpOptions = GDALCreateWarpOptions();
    if (s_srs.empty()) {
        // try reading from dataset
        s_srs = std::string(in->GetProjectionRef());
        if (s_srs.empty()) {
            GCBS_DEBUG("Failed to read input crs; assuming EPSG:4326");
            s_srs = "EPSG:4326";
        }
    }
    //    if (std::string(in->GetProjectionRef()).empty()) {
    //        if (in->SetProjection(s_srs.c_str()) != CE_None) {
    //            GCBS_DEBUG("Failed to overwrite projection for dataset.");
    //        }
    //    }
    //    std::string s = (in->GetProjectionRef());
    //    GCBS_DEBUG(s);
    psWarpOptions->hSrcDS = in;
    psWarpOptions->hDstDS = out;
    psWarpOptions->pfnProgress = GDALDummyProgress;

    // see https://github.com/OSGeo/gdal/blob/64cf9b4e889c93e34177237665fe842186d1f581/alg/gdaltransformer.cpp#L1274C5-L1275C67

    CPLStringList trnsfrm_opts;
    trnsfrm_opts.AddNameValue("SRC_SRS", s_srs.c_str());
    trnsfrm_opts.AddNameValue("GEOLOC_USE_TEMP_DATASETS", "NO");
    GDALGenImgProjTransformInfo *psInfo = static_cast<GDALGenImgProjTransformInfo *>(GDALCreateGenImgProjTransformer2(in, out, trnsfrm_opts.List()));
    if (!psInfo) {
        GCBS_ERROR("Cannot find coordinate transformation from input image to target data cube");
        throw std::string("Cannot find coordinate transformation from input image to target data cube");
    }

    
    // Overwrite reprojection transform if needed, to benefit from cache
    OGRSpatialReference srs_in;
    srs_in.SetFromUserInput(s_srs.c_str());
    void *tmparg = psInfo->pReprojectArg; // must be freed later
    if (!srs_in.IsSame(&srs_out)) {
        psInfo->pReprojectArg = gdal_transformation_cache::instance()->get(s_srs, t_srs);
        psInfo->pReproject = reproject;
    }
    
    psWarpOptions->pTransformerArg =  (void*)psInfo;
    psWarpOptions->pfnTransformer = GDALGenImgProjTransform;

    // Derive best overview level to use
    int n_ov = in->GetRasterBand(1)->GetOverviewCount();
    if (config::instance()->get_gdal_use_overviews() && n_ov > 0) {
        double *x = (double *)std::malloc(sizeof(double) * 4);
        double *y = (double *)std::malloc(sizeof(double) * 4);
        int *succ = (int *)std::malloc(sizeof(int) * 4);
        x[0] = 0;
        x[1] = 0;
        x[2] = ts_x;
        x[3] = ts_x;
        y[0] = ts_y;
        y[1] = 0;
        y[2] = ts_y;
        y[3] = 0;

        // Transform corners
        GDALGenImgProjTransform(psWarpOptions->pTransformerArg, 1, 4, x, y, NULL, succ);

        double minx = std::min(std::min(x[0], x[1]), std::min(x[2], x[3]));
        double maxx = std::max(std::max(x[0], x[1]), std::max(x[2], x[3]));
        double target_ratio = (maxx - minx) / double(ts_x);

        int16_t ilevel = 0;
        while (ilevel < n_ov) {
            double ov_ratio = double(in->GetRasterBand(1)->GetXSize()) / double(in->GetRasterBand(1)->GetOverview(ilevel)->GetXSize());
            if (ov_ratio > target_ratio) {
                --ilevel;
                break;
            }
            ++ilevel;
        }
        if (ilevel >= n_ov) {
            ilevel = n_ov - 1;
        }
        if (ilevel >= 0) {
            //GCBS_TRACE("Using overview level" + std::to_string(ilevel));
            char **oo = nullptr;
            oo = CSLAddString(oo, ("OVERVIEW_LEVEL=" + std::to_string(ilevel)).c_str());

            std::string descr = in->GetDescription();
            GDALClose(in);
            in = (GDALDataset *)GDALOpenEx(descr.c_str(), GDAL_OF_RASTER | GDAL_OF_READONLY, NULL, oo, NULL);
            if (in != NULL) {
                destroy_transform((gdalwarp_client::gdalcubes_transform_info *)psWarpOptions->pTransformerArg);
                psWarpOptions->pTransformerArg = create_transform(in, out, s_srs, t_srs);
                psWarpOptions->hSrcDS = in;  // TODO: close in_ov
            } else {
                GCBS_WARN("Failed to open GDAL overview dataset for '" + descr + "', using original full resolution image.");
            }
            CSLDestroy(oo);
        }

        std::free(x);
        std::free(y);
        std::free(succ);
    }

    psWarpOptions->eResampleAlg = GDALResampleAlg::GRA_NearestNeighbour;
    if (resampling == "bilinear") {
        psWarpOptions->eResampleAlg = GDALResampleAlg::GRA_Bilinear;
    } else if (resampling == "cubic") {
        psWarpOptions->eResampleAlg = GDALResampleAlg::GRA_Cubic;
    } else if (resampling == "cubicspline") {
        psWarpOptions->eResampleAlg = GDALResampleAlg::GRA_CubicSpline;
    } else if (resampling == "lanczos") {
        psWarpOptions->eResampleAlg = GDALResampleAlg::GRA_Lanczos;
    } else if (resampling == "average") {
        psWarpOptions->eResampleAlg = GDALResampleAlg::GRA_Average;
    } else if (resampling == "mode") {
        psWarpOptions->eResampleAlg = GDALResampleAlg::GRA_Mode;
    } else if (resampling == "max") {
        psWarpOptions->eResampleAlg = GDALResampleAlg::GRA_Max;
    } else if (resampling == "min") {
        psWarpOptions->eResampleAlg = GDALResampleAlg::GRA_Min;
    } else if (resampling == "med") {
        psWarpOptions->eResampleAlg = GDALResampleAlg::GRA_Med;
    } else if (resampling == "q1") {
        psWarpOptions->eResampleAlg = GDALResampleAlg::GRA_Q1;
    } else if (resampling == "q3") {
        psWarpOptions->eResampleAlg = GDALResampleAlg::GRA_Q3;
    }

    psWarpOptions->nBandCount = in->GetRasterCount();
    psWarpOptions->panSrcBands = (int *)CPLMalloc(sizeof(int) * psWarpOptions->nBandCount);
    psWarpOptions->panDstBands = (int *)CPLMalloc(sizeof(int) * psWarpOptions->nBandCount);
    double *dst_nodata = (double *)CPLMalloc(sizeof(double) * psWarpOptions->nBandCount);

    // Due to issues on Windows with GDAL 2.2.3 if imag no data values are not set,
    // we explicitly set then to 0 in the following. This might be removed in the future.
    double *dst_nodata_img = (double *)CPLMalloc(sizeof(double) * psWarpOptions->nBandCount);
    for (uint16_t i = 0; i < psWarpOptions->nBandCount; ++i) {
        dst_nodata[i] = NAN;
        dst_nodata_img[i] = 0.0;
        psWarpOptions->panSrcBands[i] = i + 1;
        psWarpOptions->panDstBands[i] = i + 1;
    }
    double *src_nodata = nullptr;
    double *src_nodata_img = nullptr;

    if (!srcnodata.empty()) {
        src_nodata = (double *)CPLMalloc(sizeof(double) * psWarpOptions->nBandCount);
        src_nodata_img = (double *)CPLMalloc(sizeof(double) * psWarpOptions->nBandCount);
        if (srcnodata.size() == 1) {
            for (uint16_t i = 0; i < psWarpOptions->nBandCount; ++i) {
                src_nodata[i] = srcnodata[0];
                src_nodata_img[i] = 0.0;
            }
            psWarpOptions->padfSrcNoDataReal = src_nodata;
            psWarpOptions->padfSrcNoDataImag = src_nodata_img;
        } else if (srcnodata.size() == (uint16_t)psWarpOptions->nBandCount) {
            for (uint16_t i = 0; i < psWarpOptions->nBandCount; ++i) {
                src_nodata[i] = srcnodata[i];
                src_nodata_img[i] = 0.0;
            }
            psWarpOptions->padfSrcNoDataReal = src_nodata;
            psWarpOptions->padfSrcNoDataImag = src_nodata_img;
        } else {
            GCBS_DEBUG("Number of no data values does not match number of bands of source dataset, no data values will be ignored");
            std::free(src_nodata);
            std::free(src_nodata_img);
        }
    }
    psWarpOptions->padfDstNoDataReal = dst_nodata;
    psWarpOptions->padfDstNoDataImag = dst_nodata_img;

    char **wo = nullptr;
    wo = CSLAddString(wo, "INIT_DEST=nan");
    wo = CSLAddString(wo, ("NUM_THREADS=" + std::to_string(config::instance()->get_gdal_num_threads())).c_str());
    psWarpOptions->papszWarpOptions = wo;

    // Initialize and execute the warp operation.
    GDALWarpOperation oOperation;
    if (oOperation.Initialize(psWarpOptions) != CE_None) {
        GCBS_ERROR("Initialization of gdalwarp failed");
    }

    if(oOperation.ChunkAndWarpImage(0, 0,
                                 GDALGetRasterXSize(out),
                                 GDALGetRasterYSize(out)) != CE_None) {
        GCBS_DEBUG("gdalwarp returned error code #" + std::to_string(CPLGetLastErrorNo()) + ": " + std::string(CPLGetLastErrorMsg()) + 
                   "[" + std::string(in->GetDescription()) + "]");

    }

    static_cast<GDALGenImgProjTransformInfo *>(psWarpOptions->pTransformerArg)->pReprojectArg = tmparg; // safe release
    GDALDestroyGenImgProjTransformer(psWarpOptions->pTransformerArg); 
    GDALDestroyWarpOptions(psWarpOptions);

    CPLFree(wkt_out);

    if (in) {
        GDALClose(in);
    }
    return out;
}







}  // namespace gdalcubes






