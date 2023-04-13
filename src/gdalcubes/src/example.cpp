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

/**
 * THIS FILE IS JUST FOR TESTING WITH SOME LOCAL FILES INCLUDING ABSOLUTE PATHS.
 * IT WILL NOT RUN ON YOUR MACHINE:
 */

#include <fstream>
#include <iostream>

#include "cube_factory.h"
#include "gdalcubes.h"
#include "image_collection_ops.h"

using namespace gdalcubes;

std::vector<std::string> string_list_from_text_file(std::string filename) {
    std::vector<std::string> out;

    std::string line;
    std::ifstream infile(filename);
    while (std::getline(infile, line))
        out.push_back(line);
    return out;
}

int main(int argc, char* argv[]) {
    config::instance()->gdalcubes_init();
    config::instance()->set_error_handler(error_handler::error_handler_debug);
    config::instance()->set_default_progress_bar(std::make_shared<progress_simple_stdout_with_time>());
    //config::instance()->set_default_progress_bar(std::make_shared<progress_none>());
    config::instance()->set_default_chunk_processor(std::make_shared<chunk_processor_multithread>(1));

    //    auto prsts = collection_format::list_presets();
    //    for (auto it = prsts.begin(); it != prsts.end(); ++it) {
    //        std::cout << it->first << "    " << it->second << std::endl;
    //    }

    try {
        timer t0;

        /**************************************************************************/
        // Test create image collection
        //        {
        //                        collection_format f("Sentinel2_L2A");
        //                        auto ic = image_collection::create(f, image_collection::unroll_archives(string_list_from_text_file("/home/marius/eodata/Sentinel2/file_list.txt")), false);
        //                        ic->write("test.db");
        //                        std::cout << ic->to_string() << std::endl;
        //                        std::dynamic_pointer_cast<cube_view>(image_collection_cube::create(ic)->st_reference())->write_json("view_default.json");
        //        }
        /**************************************************************************/

        /**************************************************************************/
        // Test addo
        //        {
        //            //            auto ic = std::make_shared<image_collection>("test.db");
        //            //            image_collection_ops::create_overviews(ic);
        //        }
        /**************************************************************************/

        //cube_view v = cube_view::read_json("view.json");

        /**************************************************************************/
        // test fill_time
        //        {
        //            auto c = image_collection_cube::create("test.db", v);
        //            auto cb = select_bands_cube::create(c, std::vector<std::string>{"B04", "B08"});
        //            cb->write_netcdf_file("test_fill_time_in.nc");
        //            auto cf = fill_time_cube::create(cb, "linear");
        //            cf->write_netcdf_file("test_fill_time_linear.nc");
        //        }
        /**************************************************************************/

        /**************************************************************************/
        // test reduction
        //        {
        //            auto c = image_collection_cube::create("test.db", v);
        //            auto cb = select_bands_cube::create(c, std::vector<std::string>{"B04", "B08"});
        //            cb->write_netcdf_file("band_select.nc");
        //            auto cr = reduce_cube::create(cb, "median");
        //            cr->write_gdal_image("test_A.tif");
        //
        //            c = image_collection_cube::create("test.db", v);
        //            cr = reduce_cube::create(c, "max");
        //            cb = select_bands_cube::create(cr, std::vector<std::string>{"B04_max", "B08_max"});
        //            reduce_cube::create(cb, "max")->write_gdal_image("test_B.tif");
        //        }
        /**************************************************************************/

        /**************************************************************************/
        // test packed export
        {
            //            auto c = image_collection_cube::create("test.db", v);

            //            auto cb = select_bands_cube::create(c, std::vector<std::string>{"B04", "B08"});
            //            cb->write_tif_collection("/home/marius/Desktop/test_pack1",
            //                                     "", true, true, std::map<std::string, std::string>(), "NEAREST", packed_export::make_uint8(1, 0));
        }
        /**************************************************************************/

        /**************************************************************************/
        // test masking
        //        {
        //                    cube_view w;
        //                    w.left() = 300000.000;
        //                    w.top() = 5800020.000;
        //                    w.bottom() = 5690220.000;
        //                    w.right() = 409800.000;
        //                    w.srs() = "EPSG:32632";
        //                    w.nx() = 500;
        //                    w.ny() = 500;
        //                    w.dt(duration::from_string("P1D"));
        //                    w.t0() = datetime::from_string("2018-06-14");
        //                    w.t1() = datetime::from_string("2018-06-14");
        //
        //                    auto c = image_collection_cube::create("test.db", w);
        //                    std::shared_ptr<image_mask> mask = std::make_shared<value_mask>(std::unordered_set<double>{8, 9});
        //                    c->set_mask("SCL", mask);
        //                    auto cb = select_bands_cube::create(c, std::vector<std::string>{"SCL", "B08"});
        //
        //                    cb->write_netcdf_file("mask.nc");
        //  }

        //auto v = cube_view::read_json("view.json");

        /**************************************************************************/

        /**************************************************************************/
        // test filter predicate over time
        //        {
        //            auto c = image_collection_cube::create("test.db", v);
        //            auto cb = select_bands_cube::create(c, std::vector<std::string>{"B04", "B08"});
        //            auto cf = filter_pixel_cube::create(cb, "B04 < 1000");
        //            auto cr = reduce_time_cube::create(cf, {{"max", "B04"}, {"count", "B04"}});
        //            cr->write_gdal_image("test_filter_predicate_2.tif");
        //        }
        /**************************************************************************/

        //        /**************************************************************************/
        //        // test reduction over space
        //        {
        //            auto c = image_collection_cube::create("test.db", v);
        //            auto cb = select_bands_cube::create(c, std::vector<std::string>{"B04", "B08"});
        //            auto cr = reduce_space_cube::create(cb, {{"count", "B04"}, {"mean", "B04"}});
        //            //auto cr = reduce_time_cube::create(cb, {{"min", "B04"}, {"max", "B04"}, {"median", "B04"}, {"var", "B04"}, {"which_min", "B04"}, {"which_max", "B04"}});
        //            cr->write_netcdf_file("test_reduce_space.nc");
        //        }
        //        /**************************************************************************/

        //        /**************************************************************************/
        //        // test window time over space
        //        {
        //            auto c = image_collection_cube::create("test.db", v);
        //            c->st_reference()->nt(5);
        //            auto cb = select_bands_cube::create(c, std::vector<std::string>{"B04"});
        //            auto cw = window_time_cube::create(cb, {{"mean", "B04"}}, 1, 1);
        //            //   cw->write_netcdf_file("test_window_time_reduce.nc");
        //
        //            auto cw2 = window_time_cube::create(cb, {-1.0, 2, -1.0}, 1, 1);
        //            cw2->write_netcdf_file("test_window_time_kernel.nc");
        //        }
        //      /**************************************************************************/

        /**************************************************************************/
        // Test apply_pixel
        //                {
        //                                        auto c = image_collection_cube::create("test.db", v);
        //                                        //auto capply_err = apply_pixel_cube::create(select_bands_cube::create(c, std::vector<std::string>({"B04", "B08"})), {"(B08 - B04)/(B08 + B04 -c Bsss)"});
        //                                        //auto capply = apply_pixel_cube::create(select_bands_cube::create(c, std::vector<std::string>({"B04", "B08"})), {"(B08 - B04)/(B08 + B04)"});
        //                                        //auto capply = apply_pixel_cube::create(select_bands_cube::create(c, std::vector<std::string>({"B02", "B03", "B04"})), {"sqrt((B02+B03+B04)^2)"});
        //                                        // auto capply = apply_pixel_cube::create(select_bands_cube::create(c, std::vector<std::string>({"B02", "B03", "B04"})), {"B02/B03"});
        //
        //                                        auto capply = apply_pixel_cube::create(c, {"(B08 - B04)/(B08 + B04)"});
        //
        //                                        auto cr = reduce_cube::create(capply, "median");
        //                                        // cr->write_gdal_image("test_apply_reduce.tif");
        //                                        cr->write_netcdf_file("test_apply_reduce.nc");
        //                }

        // Test apply_pixel
        //        {
        //                    auto c = dummy_cube::create(v, 1, 1.0);
        //                    auto capply = apply_pixel_cube::create(c, {"day(t0)"});
        //                    auto cr = reduce_cube::create(capply, "median");
        //                    cr->write_netcdf_file("test_apply_reduce.nc");
        //        }

        // Test query_points
        //        {
        //            auto c = image_collection_cube::create("test.db", v);
        //            //c->set_chunk_size(1000, 1000, 1000);  // single chunk
        //            auto cb = select_bands_cube::create(c, std::vector<std::string>{"B04"});
        //
        //            auto ca = apply_pixel_cube::create(cb, {"left", "top"}, {"x", "y"}, true);
        //            //cb->write_netcdf_file("cube.nc");
        //            auto cr = reduce_time_cube::create(ca, {{"median", "B04"}, {"median", "x"}, {"median", "y"}});
        //            cr->write_netcdf_file("cube_red_xy_1.nc");
        //
        //            std::vector<double> ux = {653469.0, 693953.1, 734437.2, 774921.3, 815405.4, 855889.6, 896373.7, 936857.8, 977341.9, 1017826.0};
        //            std::vector<double> uy = {6670781, 6692220, 6713658, 6735097, 6756536, 6777974, 6799413, 6820852, 6842290, 6863729};
        //            std::vector<std::string> ut = {"2018-06-04", "2018-06-12", "2018-06-20", "2018-06-30", "2018-07-10"};
        //
        //            std::vector<double> x;
        //            std::vector<double> y;
        //            std::vector<std::string> t;
        //
        //            for (uint16_t ix = 0; ix < ux.size(); ++ix) {
        //                for (uint16_t iy = 0; iy < uy.size(); ++iy) {
        //                    for (uint16_t it = 0; it < ut.size(); ++it) {
        //                        x.push_back(ux[ix]);
        //                        y.push_back(uy[iy]);
        //                        t.push_back(ut[it]);
        //                    }
        //                }
        //            }
        //            std::string srs = "EPSG:3857";
        //
        //            auto res = vector_queries::query_points(cr, x, y, t, srs);
        //
        //            //std::cout << std::fixed << std::setprecision(std::numeric_limits<double>::digits10 + 1);
        //            for (uint32_t ip = 0; ip < x.size(); ++ip) {
        //                std::cout << ip + 1 << "\t";
        //                std::cout << x[ip] << "\t";
        //                std::cout << y[ip] << "\t";
        //                std::cout << t[ip] << "\t";
        //                for (uint32_t ib = 0; ib < res.size(); ++ib) {
        //                    std::cout << res[ib][ip] << "\t";
        //                }
        //                std::cout << x[ip] - res[1][ip] << "\t";
        //                std::cout << y[ip] - res[2][ip] << "\t";
        //                std::cout << std::endl;
        //            }
        //        }

        // test zonal statistics
        //                {
        //                    cube_view w;
        //                    w.left() = -180;
        //                    w.top() = 50;
        //                    w.bottom() = -50;
        //                    w.right() = 180;
        //                    w.srs() = "EPSG:4326";
        //                    w.dx(0.2);
        //                    w.dy(0.2);
        //                    w.t0() = datetime::from_string("2019-01-01");
        //                    w.t1() = datetime::from_string("2019-01-01");
        //                    w.dt(duration::from_string("P1D"));
        //                    w.resampling_method() = resampling::resampling_type::RSMPL_AVERAGE;
        //
        //                    //            auto ch = _helper_cube::create(w);
        //                    //            ch->set_chunk_size(1,100,100);
        //                    //            ch->write_netcdf_file("/home/marius/Desktop/gdalcubes_model.nc");
        //
        //                    auto c = dummy_cube::create(w, 1, 1.0);
        //
        //                    //vector_queries::zonal_statistics(c, "/home/marius/sciebo/global_grid_5deg.gpkg", {{"count", "band1"}}, "/tmp/zonal_stats", true);
        //
        //                    auto c1 = dummy_cube::create(w, 1, 1.0);
        //                    auto c2 = apply_pixel_cube::create(c1, {"left", "top"}, {"left", "top"}, false);
        //
        //                    //vector_queries::zonal_statistics(c2,"/home/marius/sciebo/test_features.gpkg",{{"min","left"},{"max","left"},{"mean","left"},{"min","top"},{"max","top"},{"mean","top"}}, "/tmp/zonal_stats_coords_");
        //
        //                    // Real world data (CHIRPS)
        //                    collection_format f("/home/marius/github/collection_formats/formats/CHIRPS_v2_0_daily_p05_tif.json");
        //                    std::vector<std::string> files;
        //                    filesystem::iterate_directory("/home/marius/eodata/CHIRPS/", [&files](const std::string& s) {
        //                        if (s.find(".tif") != s.npos) {
        //                            files.push_back(s);
        //                        }
        //                    });
        //                    auto ic = image_collection::create(f, files, false);
        //                    //ic->write("CHIRPS.db");
        //
        //                    w.t0() = datetime::from_string("2018-01-01");
        //                    w.t1() = datetime::from_string("2018-01-04");
        //
        //                    auto chirps_cube = image_collection_cube::create("CHIRPS.db", w);
        //                    chirps_cube->set_chunk_size(10, 256, 256);
        //
        //                    // chirps_cube->write_netcdf_file("/home/marius/sciebo/chirps.nc");
        //
        //                    vector_queries::zonal_statistics(chirps_cube, "/home/marius/sciebo/world_polygons.gpkg", {{"mean", "precipitation"}}, "/tmp/zonal_stats_chirps2.gpkg", true);
        //
        //                }

        // test collection format for spacetime GDAL datasets
        {
            //            std::vector<std::string> files = {"/home/marius/Desktop/pr_day_HadGEM2-ES_historical_r1i1p1_EWEMBI_19810101-19901231.nc4",
            //                                              "/home/marius/Desktop/pr_day_HadGEM2-ES_historical_r1i1p1_EWEMBI_19910101-20001231.nc4",
            //                                              "/home/marius/Desktop/pr_day_HadGEM2-ES_historical_r1i1p1_EWEMBI_20010101-20051231.nc4"};
            //
            //            collection_format f("/home/marius/sciebo/pr_day_HadGEM2-ES_historical_r1i1p1.json");
            //            auto ic = image_collection::create(f, files, false);
            //            ic->write("/home/marius/Desktop/test.db");

            //            std::ifstream i("/tmp/cube.json");
            //            std::stringstream buf;
            //            buf << i.rdbuf();
            //            std::string err;
            //            json11::Json j = json11::Json::parse(buf.str(),err);
            //            auto cube = cube_factory::instance()->create_from_json(j);
            //            cube->write_netcdf_file("/tmp/cube.nc");
        }

        /**************************************************************************/
        // test streaming
        //        {
        //            // test streaming
        //            auto c = image_collection_cube::create("test.db", v);
        //            auto sc = stream_cube::create(c, "Rscript --vanilla stream_example.R", "stdout");
        //            auto cr = reduce_cube::create(sc, "median");
        //            cr->write_gdal_image("test_stream.tif");
        //        }

        //
        //        /******************************************/
        //        // Test CHIRPS dataset and update_view()
        //        {
        //            chdir("/home/marius/Desktop/CHIRPS/");
        //            config::instance()->set_default_chunk_processor(std::make_shared<chunk_processor_multithread>(1));
        //            auto x = image_collection_cube::create("/home/marius/Desktop/CHIRPS/CHIRPS.db", "/home/marius/Desktop/CHIRPS/view_debug.json");
        //            std::cout << x->view()->write_json_string() << std::endl;
        //            auto xmax = reduce_cube::create(x, "max");
        //            std::shared_ptr<cube_st_reference> vv = x->st_reference();
        //            vv->nt(1);
        //            vv->nx() = 100;
        //            vv->ny() = 100;
        //            xmax->update_st_reference(vv);
        //            xmax->write_gdal_image("test_max.tif");
        //        }
        //        /******************************************/

        //
        //        /******************************************/
        // Test NetCDF export
        {
            //            chdir("/home/marius/Desktop/MODIS/MOD13A3.A2018");
            //            auto cc = image_collection_cube::create("MOD13A3.db");
            //            cc->view()->aggregation_method() = aggregation::aggregation_type::AGG_MEDIAN;
            //            cc->write_netcdf_file("full.nc");
        }

        //setenv("GDAL_DISABLE_READDIR_ON_OPEN", "TRUE", 1);
        //        setenv("CPL_DEBUG", "ON", 1);
        //        setenv("CPL_LOG_ERRORS", "ON", 1);
        //        setenv("CPL_LOG", "/tmp/gdal.log", 1);

        // cube_factory::instance()->create_from_json_file("/tmp/cube.json")->write_netcdf_file("/tmp/cube.nc");

        //        int i=1;
        //        filesystem::iterate_directory("/tmp/RtmpRgZMzS/gdalcubes_debug", [&i](const std::string& f)->void {
        //            if (filesystem::extension(f) == "json") {
        //                std::cout << i << ": " << f << "...";
        //                cube_factory::instance()->create_from_json_file(f)->write_netcdf_file("/tmp/cube" + std::to_string(i) + ".nc");
        //                std::cout << std::endl;
        //                ++i;
        //            }
        //
        //        });

        //        std::cout << datetime::from_string("2020-01-01T04:56:22").to_string() << std::endl;
        //        std::cout << datetime::from_string("2020-01-01T04:56:22Z").to_string() << std::endl;
        //        std::cout << datetime::from_string("2020-01-01T04:56:22+11").to_string() << std::endl;
        //        std::cout << datetime::from_string("2020-01-01T04:56:22-11").to_string() << std::endl;
        //        std::cout << datetime::from_string("2020-01-01T04:56:22+11:00").to_string() << std::endl;
        //        std::cout << datetime::from_string("2020-01-01T04:56:22-11:00").to_string() << std::endl;
        //        std::cout << datetime::from_string("2020-01-01T04:56:22+1100").to_string() << std::endl;
        //        std::cout << datetime::from_string("2020-01-01T04:56:22-1100").to_string()  << std::endl;

       // auto c = cube_factory::instance()->create_from_json_file("/tmp/cube.json");
       // auto c_agg = aggregate_time_cube::create(c, "P1M", "mean");

        //c_agg->write_netcdf_file("/tmp/cube.nc");


//        cube_view r;
//        r.srs("EPSG:3857");
//        r.left(-6180000);
//        r.right(-6080000);
//        r.top(-450000);
//        r.bottom(-550000);
//        r.nx(100);
//        r.ny(100);
//        r.t0(datetime::from_string("2014-01-01"));
//        r.t1(datetime::from_string("2014-12-31"));
//        r.dt(duration::from_string("P1D"));
//
//        auto c = dummy_cube::create(r, 1, 1.0);
//        auto ca = crop_cube::create(c, 0, 10, 0, 10, 0, 11);
//        std::cout <<  ca->st_reference()->nx() << "," << ca->st_reference()->ny() << "," <<   ca->st_reference()->nt() <<  std::endl; // 29
//
//        auto ca1 = crop_cube::create(c, r.left(), r.right(), r.bottom(), r.top(), "2014-01-01", "2014-01-20", "in");
//        std::cout <<  ca1->st_reference()->nx() << "," << ca1->st_reference()->ny() << "," <<   ca1->st_reference()->nt() <<  std::endl; // 29
//
//
//        // check snapping
//        ca1 = crop_cube::create(c, r.left() + 0.01, r.right() - 0.01, r.bottom() + 0.01, r.top() - 0.01, "2014-01-01", "2014-01-20", "in");
//        std::cout <<  ca1->st_reference()->nx() << "," << ca1->st_reference()->ny() << "," <<   ca1->st_reference()->nt() <<  std::endl; // 29
//        std::cout <<  std::dynamic_pointer_cast<cube_stref_regular>(ca1->st_reference())->t0().to_string()<< "," << std::dynamic_pointer_cast<cube_stref_regular>(ca1->st_reference())->t1().to_string() <<  std::endl;
//
//        ca1 = crop_cube::create(c, r.left() + 0.01, r.right() - 0.01, r.bottom() + 0.01, r.top() - 0.01, "2014-01-01", "2014-01-20", "near");
//        std::cout <<  ca1->st_reference()->nx() << "," << ca1->st_reference()->ny() << "," <<   ca1->st_reference()->nt() <<  std::endl; // 29
//        std::cout <<  std::dynamic_pointer_cast<cube_stref_regular>(ca1->st_reference())->t0().to_string()<< "," << std::dynamic_pointer_cast<cube_stref_regular>(ca1->st_reference())->t1().to_string() <<  std::endl;
//
//        ca1 = crop_cube::create(c, r.left() + 0.01, r.right() - 0.01, r.bottom() + 0.01, r.top() - 0.01, "2014-01-01", "2014-01-20", "out");
//        std::cout <<  ca1->st_reference()->nx() << "," << ca1->st_reference()->ny() << "," <<   ca1->st_reference()->nt() <<  std::endl; // 29
//        std::cout <<  std::dynamic_pointer_cast<cube_stref_regular>(ca1->st_reference())->t0().to_string()<< "," << std::dynamic_pointer_cast<cube_stref_regular>(ca1->st_reference())->t1().to_string() <<  std::endl;
//
//
//
//
//        auto c1 = aggregate_time_cube::create(c, 3);
//        ca1 = crop_cube::create(c1, r.left(), r.right(), r.bottom(), r.top(), "2014-01-02", "2014-01-20", "in");
//        std::cout <<  ca1->st_reference()->nx() << "," << ca1->st_reference()->ny() << "," <<   ca1->st_reference()->nt() <<  std::endl; // 29
//        std::cout <<  std::dynamic_pointer_cast<cube_stref_regular>(ca1->st_reference())->t0().to_string()<< "," << std::dynamic_pointer_cast<cube_stref_regular>(ca1->st_reference())->t1().to_string() <<  std::endl;
//
//        ca1 = crop_cube::create(c1, r.left(), r.right(), r.bottom(), r.top(), "2014-01-02", "2014-01-20", "out");
//        std::cout <<  ca1->st_reference()->nx() << "," << ca1->st_reference()->ny() << "," <<   ca1->st_reference()->nt() <<  std::endl; // 29
//        std::cout <<  std::dynamic_pointer_cast<cube_stref_regular>(ca1->st_reference())->t0().to_string()<< "," << std::dynamic_pointer_cast<cube_stref_regular>(ca1->st_reference())->t1().to_string() <<  std::endl;
//
//        ca1 = crop_cube::create(c1, r.left(), r.right(), r.bottom(), r.top(), "2014-01-02", "2014-01-20", "near");
//        std::cout <<  ca1->st_reference()->nx() << "," << ca1->st_reference()->ny() << "," <<   ca1->st_reference()->nt() <<  std::endl; // 29
//        std::cout <<  std::dynamic_pointer_cast<cube_stref_regular>(ca1->st_reference())->t0().to_string()<< "," << std::dynamic_pointer_cast<cube_stref_regular>(ca1->st_reference())->t1().to_string() <<  std::endl;
//
//
//        auto c2 = select_time_cube::create(c, {"2014-01-02", "2014-01-06", "2014-01-20", "2014-01-21", "2014-01-28", "2014-02-01"});
//        ca1 = crop_cube::create(c2, r.left(), r.right(), r.bottom(), r.top(), "2014-01-05", "2014-01-27", "near");
//        std::cout <<  ca1->st_reference()->nx() << "," << ca1->st_reference()->ny() << "," <<   ca1->st_reference()->nt() <<  std::endl; // 29
//        std::cout <<  std::dynamic_pointer_cast<cube_stref_regular>(ca1->st_reference())->t0().to_string()<< "," << std::dynamic_pointer_cast<cube_stref_regular>(ca1->st_reference())->t1().to_string() <<  std::endl;
//
//        ca1 = crop_cube::create(c2, r.left(), r.right(), r.bottom(), r.top(), "2014-01-05", "2014-01-27", "out");
//        std::cout <<  ca1->st_reference()->nx() << "," << ca1->st_reference()->ny() << "," <<   ca1->st_reference()->nt() <<  std::endl; // 29
//        std::cout <<  std::dynamic_pointer_cast<cube_stref_regular>(ca1->st_reference())->t0().to_string()<< "," << std::dynamic_pointer_cast<cube_stref_regular>(ca1->st_reference())->t1().to_string() <<  std::endl;
//
//        ca1 = crop_cube::create(c2, r.left(), r.right(), r.bottom(), r.top(), "2014-01-05", "2014-01-27", "in");
//        std::cout <<  ca1->st_reference()->nx() << "," << ca1->st_reference()->ny() << "," <<   ca1->st_reference()->nt() <<  std::endl; // 29
//        std::cout <<  std::dynamic_pointer_cast<cube_stref_regular>(ca1->st_reference())->t0().to_string()<< "," << std::dynamic_pointer_cast<cube_stref_regular>(ca1->st_reference())->t1().to_string() <<  std::endl;


        // TODO: test for labeled time axis


//        cube_stref_regular r;
//        r.srs("EPSG:3857");
//        r.left(-6180000);
//        r.right(-6080000);
//        r.top(-450000);
//        r.bottom(-550000);
//        r.dx(1000);
//        r.dy(1000);
//        r.t0(datetime::from_string("2014-01-01"));
//        r.t1(datetime::from_string("2014-12-31"));
//        r.dt(duration::from_string("P10D"));
//
//
//        std::cout <<  r.index_at_datetime(datetime::from_string("2014-10-18")) << std::endl; // 29
//        std::cout <<  r.index_at_datetime(datetime::from_string("2014-10-19")) << std::endl; // 29
//        std::cout <<  r.index_at_datetime(datetime::from_string("2014-10-20")) << std::endl; // 29
//        std::cout <<  r.index_at_datetime(datetime::from_string("2014-10-21")) << std::endl; // 29
//        std::cout <<  r.index_at_datetime(datetime::from_string("2014-10-22")) << std::endl; // 29
//        std::cout <<  r.index_at_datetime(datetime::from_string("2014-10-23")) << std::endl; // 29
//        std::cout <<  r.index_at_datetime(datetime::from_string("2014-10-24")) << std::endl; // 29
//        std::cout <<  r.index_at_datetime(datetime::from_string("2014-10-25")) << std::endl; // 29
//        std::cout <<  r.index_at_datetime(datetime::from_string("2014-10-26")) << std::endl; // 29
//        std::cout <<  r.index_at_datetime(datetime::from_string("2014-10-27")) << std::endl; // 29
//        std::cout <<  r.index_at_datetime(datetime::from_string("2014-10-28")) << std::endl; // 30?
//
//        std::cout <<  r.index_at_datetime(datetime::from_string("2014-11-01")) - 1 << std::endl; // 30?


        //cube_factory::instance()->create_from_json_file("/tmp/cube.json")->write_png_collection("/tmp/testpng_error","ttt", {"B03","B04","B05"}, {0,0,0}, {1500, 1500, 1500},1,{},true);
        //cube_factory::instance()->create_from_json_file("/tmp/cube.json")->write_netcdf_file("/tmp/test.nc");

        //cube_factory::instance()->create_from_json_file("/tmp/cube.json")->write_png_collection("/tmp/testpng_rgb_alpha","ttt", {"B04","B03","B02"}, {0,0,0}, {1500, 1500, 1500},1,{},true);
        //cube_factory::instance()->create_from_json_file("/tmp/cube.json")->write_png_collection("/tmp/testpng_rgb_nacol","ttt", {"B04","B03","B02"}, {0,0,0}, {1500, 1500, 1500},0.5,{255,0,0},false);
        //cube_factory::instance()->create_from_json_file("/tmp/cube.json")->write_png_collection("/tmp/testpng_gray_nacol","ttt", {"B04"}, {100}, {2000},0.7,{255,0,0},false);
        //cube_factory::instance()->create_from_json_file("/tmp/cube.json")->write_png_collection("/tmp/testpng_gray_nacolblack","ttt", {"B04"}, {100}, {2000},0.7,{0},false);
        //cube_factory::instance()->create_from_json_file("/tmp/cube.json")->write_png_collection("/tmp/testpng_gray_alpha","ttt", {"B04"}, {100}, {2000},0.7,{255,0,0},true);
        //cube_factory::instance()->create_from_json_file("/tmp/cube.json")->write_png_collection("/tmp/testpng_rgb_default","ttt");

        //cube_factory::instance()->create_from_json_file("/tmp/cube.json")->write_png_collection("/tmp/testpng","ttt", {"B03"}, {0}, {1500},1,{255},false);
        //cube_factory::instance()->create_from_json_file("/tmp/cube.json")->write_png_collection("/tmp/testpng","ttt", {"B03"}, {0}, {1500},1,{},true);

        //test_multiprocess::write_chunks_netcdf(cube_factory::instance()->create_from_json_file("/tmp/cube.json"),"/tmp", "test");
        //config::instance()->set_gdal_use_overviews(false);
        //        auto c = cube_factory::instance()->create_from_json_file("/tmp/cube.json");
        //        c->write_netcdf_file("/tmp/xxx.nc");
        //c->write_tif_collection("/tmp/TESTTIF", "xxx");

        /******************************************/

//        std::vector<std::string> files{"/media/marius/Samsung_T5/eodata/GPM/IMERG_3B_DAY_GIS_V06A/3B-DAY-GIS.MS.MRG.3IMERG.20150101-S000000-E235959.0000.V06A/3B-DAY-GIS.MS.MRG.3IMERG.20150101-S000000-E235959.0000.V06A.total.accum.tif", "/media/marius/Samsung_T5/eodata/GPM/IMERG_3B_DAY_GIS_V06A/3B-DAY-GIS.MS.MRG.3IMERG.20150102-S000000-E235959.0030.V06A/3B-DAY-GIS.MS.MRG.3IMERG.20150102-S000000-E235959.0030.V06A.total.accum.tif", "/media/marius/Samsung_T5/eodata/GPM/IMERG_3B_DAY_GIS_V06A/3B-DAY-GIS.MS.MRG.3IMERG.20150103-S000000-E235959.0060.V06A/3B-DAY-GIS.MS.MRG.3IMERG.20150103-S000000-E235959.0060.V06A.total.accum.tif", "/media/marius/Samsung_T5/eodata/GPM/IMERG_3B_DAY_GIS_V06A/3B-DAY-GIS.MS.MRG.3IMERG.20150104-S000000-E235959.0090.V06A/3B-DAY-GIS.MS.MRG.3IMERG.20150104-S000000-E235959.0090.V06A.total.accum.tif", "/media/marius/Samsung_T5/eodata/GPM/IMERG_3B_DAY_GIS_V06A/3B-DAY-GIS.MS.MRG.3IMERG.20150105-S000000-E235959.0120.V06A/3B-DAY-GIS.MS.MRG.3IMERG.20150105-S000000-E235959.0120.V06A.total.accum.tif", "/media/marius/Samsung_T5/eodata/GPM/IMERG_3B_DAY_GIS_V06A/3B-DAY-GIS.MS.MRG.3IMERG.20150106-S000000-E235959.0150.V06A/3B-DAY-GIS.MS.MRG.3IMERG.20150106-S000000-E235959.0150.V06A.total.accum.tif", "/media/marius/Samsung_T5/eodata/GPM/IMERG_3B_DAY_GIS_V06A/3B-DAY-GIS.MS.MRG.3IMERG.20150107-S000000-E235959.0180.V06A/3B-DAY-GIS.MS.MRG.3IMERG.20150107-S000000-E235959.0180.V06A.total.accum.tif", "/media/marius/Samsung_T5/eodata/GPM/IMERG_3B_DAY_GIS_V06A/3B-DAY-GIS.MS.MRG.3IMERG.20150108-S000000-E235959.0210.V06A/3B-DAY-GIS.MS.MRG.3IMERG.20150108-S000000-E235959.0210.V06A.total.accum.tif", "/media/marius/Samsung_T5/eodata/GPM/IMERG_3B_DAY_GIS_V06A/3B-DAY-GIS.MS.MRG.3IMERG.20150109-S000000-E235959.0240.V06A/3B-DAY-GIS.MS.MRG.3IMERG.20150109-S000000-E235959.0240.V06A.total.accum.tif", "/media/marius/Samsung_T5/eodata/GPM/IMERG_3B_DAY_GIS_V06A/3B-DAY-GIS.MS.MRG.3IMERG.20150110-S000000-E235959.0270.V06A/3B-DAY-GIS.MS.MRG.3IMERG.20150110-S000000-E235959.0270.V06A.total.accum.tif", "/media/marius/Samsung_T5/eodata/GPM/IMERG_3B_DAY_GIS_V06A/3B-DAY-GIS.MS.MRG.3IMERG.20150111-S000000-E235959.0300.V06A/3B-DAY-GIS.MS.MRG.3IMERG.20150111-S000000-E235959.0300.V06A.total.accum.tif", "/media/marius/Samsung_T5/eodata/GPM/IMERG_3B_DAY_GIS_V06A/3B-DAY-GIS.MS.MRG.3IMERG.20150112-S000000-E235959.0330.V06A/3B-DAY-GIS.MS.MRG.3IMERG.20150112-S000000-E235959.0330.V06A.total.accum.tif", "/media/marius/Samsung_T5/eodata/GPM/IMERG_3B_DAY_GIS_V06A/3B-DAY-GIS.MS.MRG.3IMERG.20150113-S000000-E235959.0360.V06A/3B-DAY-GIS.MS.MRG.3IMERG.20150113-S000000-E235959.0360.V06A.total.accum.tif", "/media/marius/Samsung_T5/eodata/GPM/IMERG_3B_DAY_GIS_V06A/3B-DAY-GIS.MS.MRG.3IMERG.20150114-S000000-E235959.0390.V06A/3B-DAY-GIS.MS.MRG.3IMERG.20150114-S000000-E235959.0390.V06A.total.accum.tif", "/media/marius/Samsung_T5/eodata/GPM/IMERG_3B_DAY_GIS_V06A/3B-DAY-GIS.MS.MRG.3IMERG.20150115-S000000-E235959.0420.V06A/3B-DAY-GIS.MS.MRG.3IMERG.20150115-S000000-E235959.0420.V06A.total.accum.tif"};
//        std::vector<std::string> datetime{"2015-01-01", "2015-01-02", "2015-01-03", "2015-01-04", "2015-01-05", "2015-01-06", "2015-01-07", "2015-01-08", "2015-01-09", "2015-01-10", "2015-01-11", "2015-01-12", "2015-01-13", "2015-01-14", "2015-01-15"};
//        std::vector<std::string> band_names{"TOTAL_ACCUM"};
//
//        auto c1 = simple_cube::create(files, datetime, {}, band_names, 0.5,  0.5);
//        auto cr = reduce_time_cube::create(c1, {{"median", "TOTAL_ACCUM"}});
//        cr->write_netcdf_file("/tmp/cube1.nc", 6);

        //                    auto c2 = apply_pixel_cube::create(c1, {"left", "top"}, {"left", "top"}, false);





        /* test extract_geom */
        {
//
//            cube_view r;
//            r.srs("EPSG:3857");
//            r.left(-6180000);
//            r.right(-6080000);
//            r.top(-450000);
//            r.bottom(-550000);
//            r.nx(100);
//            r.ny(100);
//            r.t0(datetime::from_string("2014-01-01"));
//            r.t1(datetime::from_string("2014-12-31"));
//            r.dt(duration::from_string("P1D"));
//
//            auto d= dummy_cube::create(r, 2,1.2);
//            d->set_chunk_size(10,100,100);
//            //d->write_netcdf_file("/tmp/test_extract.nc");
//
//            auto ex1 = extract_geom::create(d, "/home/marius/Desktop/test_extract_cube/points1000_time.gpkg", "date");


//            auto c = cube_factory::instance()->create_from_json_file("/home/marius/Desktop/cube.json");
//            auto ex1 = extract_geom::create(c, "/home/marius/Desktop/test.gpkg", "time");
//
//            for (uint32_t i=0; i<ex1->count_chunks(); ++i) {
//                std::cout <<  std::endl <<  "CHUNK ID " << i << std::endl;
//                std::cout << "---------------------------------------------------------"<< std::endl;
//                std::shared_ptr<chunk_data> dat = ex1->read_chunk(i);
//                if (!dat->empty()) {
//                    uint32_t ncol = dat->size()[0];
//                    uint32_t nrow = dat->size()[1];
//
//                    for (uint32_t row = 0; row < nrow; ++row) {
//                        for (uint32_t col = 0; col < ncol; ++col) {
//                            std::cout << ((double*)dat->buf())[col * nrow + row] << " ";
//                        }
//                        std::cout << std::endl;
//                    }
//                }
//            }


        }


        std::vector<std::string> items;
        filesystem::iterate_directory("/home/marius/github", [&items](const std::string& f) {
            items.push_back(f);
        });
        for (auto it = items.begin(); it != items.end(); ++it) {
            std::cout << *it << std::endl;
        }
        std::cout << items.size() << " files" << std::endl;





    } catch (std::string e) {
        std::cout << e << std::endl;
    }

    config::instance()->gdalcubes_cleanup();
}