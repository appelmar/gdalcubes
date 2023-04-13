/*
    MIT License

    Copyright (c) 2019 Marius Appel <marius.appel@uni-muenster.de>

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

#include "../external/catch.hpp"
#include "../gdalcubes.h"
using namespace gdalcubes;
/**
 * THIS TEST CASE DEPENDS ON LOCAL FILES AND IS COMMENTED OUT.
 */

//TEST_CASE("Example test 1 with sample data", "[example_1]") {
//    config::instance()->gdalcubes_init();
//
//    std::string in = "HDF4_EOS:EOS_GRID:\"../../test/MOD13A2.A2015193.h23v03.006.2015304013141.hdf\":MODIS_Grid_16DAY_1km_VI:1 km 16 days NDVI";
//
//    cube_view v;
//    v.left() = 5559752.598333000205457;
//    v.right() = 6671703.118;
//    v.bottom() = 5559752.5983330002055;
//    v.top() = 6671703.117999999783933;
//    v.proj() = "+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs";
//    v.ny() = 1200;
//    v.nx() = 1200;
//
//    REQUIRE(v.ny() == 1200);
//    REQUIRE(v.nx() == 1200);
//    REQUIRE(v.dx() == 926.625433055833014);
//    REQUIRE(v.dy() == 926.625433055833014);
//
//    v.t0() = datetime::from_string("2015-07-12");
//    v.t1() = datetime::from_string("2015-07-12");
//    v.dt(duration::from_string("P1D"));
//
//    REQUIRE(v.nt() == 1);
//
//    v.aggregation_method() = aggregation::aggregation_type::AGG_MIN;
//    v.resampling_method() = resampling::resampling_type::RSMPL_NEAR;
//
//    auto ic = image_collection::create(collection_format("../../test/MOD13A2.json"), {in});
//
//    auto cube = image_collection_cube::create(ic, v);
//
//    cube->write_netcdf_file("test_example_1_write_netcdf_out.nc");
//
//    cube->write_gtiff_directory("test_example_1_write_gtiff_dir_out");
//    reduce_cube::create(cube, "min")->write_gdal_image("test_example_1_write_gdal.tif");
//
//    config::instance()->gdalcubes_cleanup();
//}
