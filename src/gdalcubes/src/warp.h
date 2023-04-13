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

#ifndef WARP_H
#define WARP_H

#include <gdal_alg.h>

#include <map>

#include "coord_types.h"

namespace gdalcubes {

/**
 * Minimal interface to the GDAL Warp API to reduce computational overhead of repeated / parallel gdalwarp calls
 */
class gdalwarp_client {
    struct gdalcubes_transform_info{  // adapted from https://github.com/OSGeo/gdal/blob/master/gdal/alg/gdaltransformer.cpp, struct GDALGenImgProjTransformInfo
        double adfSrcGeoTransform[6];
        double adfSrcInvGeoTransform[6];

        void *pReprojectArg;
        GDALTransformerFunc pReproject;

        double adfDstGeoTransform[6];
        double adfDstInvGeoTransform[6];
    };

    struct gdalcubes_reprojection_info {  // adapted from https://github.com/OSGeo/gdal/blob/master/gdal/alg/gdaltransformer.cpp, struct GDALReprojectionTransformInfo
        OGRCoordinateTransformation *poForwardTransform = nullptr;
        OGRCoordinateTransformation *poReverseTransform = nullptr;
    } ;

   public:
    /**
     * Cache for reprojection transformations given a pair of source and destination coordinate reference systems
     */
    class gdal_transformation_cache {
       public:
        static gdal_transformation_cache *instance() {
            static gdal_transformation_cache instance;
            return &instance;
        }

        gdalcubes_reprojection_info *get(std::string srs_in_str, std::string srs_out_str);

       private:
        gdal_transformation_cache(const gdal_transformation_cache &) = delete;
        gdal_transformation_cache(gdal_transformation_cache &&) = delete;
        gdal_transformation_cache &operator=(const gdal_transformation_cache &) = delete;
        gdal_transformation_cache &operator=(gdal_transformation_cache &&) = delete;
        gdal_transformation_cache() {}
        ~gdal_transformation_cache();

        std::map<std::pair<std::string, std::string>, gdalcubes_reprojection_info *> _cache;
        std::mutex _mutex;
    };

   public:
    /**
     * Warp source GDAL dataset to a target grid
     * @param in source GDAL dataset, will be closed at the end of this function
     * @param s_srs spatial reference system of source image, given as string understandable for OGRSpatialReference::SetFromUserInput()
     * @param t_srs target spatial reference system, given as string understandable for OGRSpatialReference::SetFromUserInput()
     * @param te_left left (minimum x) coordinate of the target grid, given in the target SRS
     * @param te_right right (maximum x) coordinate of the target grid, given in the target SRS
     * @param te_top left (maximum y) coordinate of the target grid, given in the target SRS
     * @param te_bottom left (minimum y) coordinate of the target grid, given in the target SRS
     * @param ts_x number of pixels of the target grid in x direction
     * @param ts_y number of pixels of the target grid in y direction
     * @param resampling  resampling method, given as a string (see https://gdal.org/programs/gdalwarp.html#cmdoption-gdalwarp-r for possible options)
     * @param srcnodata vector with no data values of the source dataset per band
     * @return A new in-memory GDALDataset object
     */
    static GDALDataset *warp(GDALDataset *in, std::string s_srs, std::string t_srs, double te_left, double te_right, double te_top, double te_bottom, uint32_t ts_x, uint32_t ts_y, std::string resampling, std::vector<double> srcnodata);

    static gdalcubes_transform_info *create_transform(GDALDataset *in, GDALDataset *out, std::string srs_in_str, std::string srs_out_str);
    static void destroy_transform(gdalcubes_transform_info *transform);

    // implements GDALTransformerFunc signature
    static int transform(void *pTransformerArg,
                         int bDstToSrc, int nPointCount,
                         double *x, double *y, double *z = nullptr, int *panSuccess = nullptr);

    static gdalcubes_reprojection_info *create_reprojection(std::string srs_in_str, std::string srs_out_str);
    static void destroy_reprojection(gdalcubes_reprojection_info *reprojection);

    // implements GDALTransformerFunc signature
    static int reproject(void *pTransformerArg,
                         int bDstToSrc, int nPointCount,
                         double *x, double *y, double *z = nullptr, int *panSuccess = nullptr);
};

}  // namespace gdalcubes
#endif  // WARP_H
