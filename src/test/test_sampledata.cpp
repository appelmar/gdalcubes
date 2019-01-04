
#include "../external/catch.hpp"
#include "../gdalcubes.h"

TEST_CASE("Example test 1 with sample data", "[example_1]") {
    config::instance()->gdalcubes_init();

    std::string in = "HDF4_EOS:EOS_GRID:\"../../test/MOD13A2.A2015193.h23v03.006.2015304013141.hdf\":MODIS_Grid_16DAY_1km_VI:1 km 16 days NDVI";

    cube_view v;
    v.left() = 5559752.598333000205457;
    v.right() = 6671703.118;
    v.bottom() = 5559752.5983330002055;
    v.top() = 6671703.117999999783933;
    v.proj() = "+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs";
    v.ny() = 1200;
    v.nx() = 1200;

    REQUIRE(v.ny() == 1200);
    REQUIRE(v.nx() == 1200);
    REQUIRE(v.dx() == 926.625433055833014);
    REQUIRE(v.dy() == 926.625433055833014);

    v.t0() = datetime::from_string("2015-07-12");
    v.t1() = datetime::from_string("2015-07-12");
    v.dt() = duration::from_string("P1D");

    REQUIRE(v.nt() == 1);

    v.aggregation_method() = aggregation::aggregation_type::AGG_MIN;
    v.resampling_method() = resampling::resampling_type::RSMPL_NEAR;

    auto ic = image_collection::create(collection_format("../../test/MOD13A2.json"), {in});

    auto cube = image_collection_cube::create(ic, v);

    cube->write_netcdf_file("test_example_1_write_netcdf_out.nc");

    cube->write_gtiff_directory("test_example_1_write_gtiff_dir_out");
    reduce_cube::create(cube, "min")->write_gdal_image("test_example_1_write_gdal.tif");

    config::instance()->gdalcubes_cleanup();
}
