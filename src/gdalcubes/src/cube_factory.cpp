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

#include "cube_factory.h"

#include <fstream>

#include "aggregate_time.h"
#include "aggregate_space.h"
#include "apply_pixel.h"
#include "crop.h"
#include "dummy.h"
#include "external/json11/json11.hpp"
#include "extract_geom.h"
#include "filesystem.h"
#include "fill_time.h"
#include "filter_geom.h"
#include "filter_pixel.h"
#include "image_collection_cube.h"
#include "join_bands.h"
#include "ncdf_cube.h"
#include "reduce_time.h"
#include "reduce_space.h"
#include "rename_bands.h"
#include "select_bands.h"
#include "select_time.h"
#include "simple_cube.h"
#include "slice_time.h"
#include "slice_space.h"
#include "stream.h"
#include "stream_apply_pixel.h"
#include "stream_apply_time.h"
#include "stream_reduce_space.h"
#include "stream_reduce_time.h"
#include "window_space.h"
#include "window_time.h"

namespace gdalcubes {

cube_factory* cube_factory::_instance = 0;

std::shared_ptr<cube> cube_factory::create_from_json_file(std::string path) {
    std::ifstream i(path);
    std::stringstream buf;
    buf << i.rdbuf();
    std::string err;
    json11::Json j = json11::Json::parse(buf.str(), err);
    return (create_from_json(j));
}

std::shared_ptr<cube> cube_factory::create_from_json(json11::Json j) {
    if (j["cube_type"].is_null()) {
        throw std::string("ERROR in cube_factory::create_from_json(): invalid object, missing cube_type key.");
    }

    std::string cube_type = j["cube_type"].string_value();

    return (cube_generators[cube_type](j));  //recursive creation
}

void cube_factory::register_cube_type(std::string type_name,
                                      std::function<std::shared_ptr<cube>(json11::Json&)> generator) {
    cube_generators.insert(std::make_pair(type_name, generator));
}

void cube_factory::register_default() {
    /* register data cube types */
    cube_generators.insert(std::make_pair<std::string, std::function<std::shared_ptr<cube>(json11::Json&)>>(
        "reduce_time", [](json11::Json& j) {
            std::vector<std::pair<std::string, std::string>> band_reducers;
            for (uint16_t i = 0; i < j["reducer_bands"].array_items().size(); ++i) {
                band_reducers.push_back(std::make_pair(j["reducer_bands"][i][0].string_value(), j["reducer_bands"][i][1].string_value()));
            }
            if (!j["names"].is_null()) {
              std::vector<std::string> names;
              for (uint16_t i = 0; i < j["names"].array_items().size(); ++i) {
                names.push_back(j["names"][i].string_value());
              }
              auto x = reduce_time_cube::create(instance()->create_from_json(j["in_cube"]), band_reducers, names);
              return x;
            }
            else {
              auto x = reduce_time_cube::create(instance()->create_from_json(j["in_cube"]), band_reducers);
              return x;
            }
            
        }));
    cube_generators.insert(std::make_pair<std::string, std::function<std::shared_ptr<cube>(json11::Json&)>>(
        "reduce_space", [](json11::Json& j) {
            std::vector<std::pair<std::string, std::string>> band_reducers;
            for (uint16_t i = 0; i < j["reducer_bands"].array_items().size(); ++i) {
                band_reducers.push_back(std::make_pair(j["reducer_bands"][i][0].string_value(), j["reducer_bands"][i][1].string_value()));
            }
            if (!j["names"].is_null()) {
              std::vector<std::string> names;
              for (uint16_t i = 0; i < j["names"].array_items().size(); ++i) {
                names.push_back(j["names"][i].string_value());
              }
              auto x = reduce_space_cube::create(instance()->create_from_json(j["in_cube"]), band_reducers, names);
              return x;
            }
            else {
              auto x = reduce_space_cube::create(instance()->create_from_json(j["in_cube"]), band_reducers);
              return x;
            }
        }));

    cube_generators.insert(std::make_pair<std::string, std::function<std::shared_ptr<cube>(json11::Json&)>>(
        "window_time", [](json11::Json& j) {
            if (!j["kernel"].is_null()) {
                std::vector<double> kernel;
                for (uint16_t i = 0; i < j["kernel"].array_items().size(); ++i) {
                    kernel.push_back(j["kernel"][i].number_value());
                }
                return window_time_cube::create(instance()->create_from_json(j["in_cube"]), kernel,
                                                j["win_size_l"].int_value(), j["win_size_r"].int_value());
            } else {
                std::vector<std::pair<std::string, std::string>> band_reducers;
                for (uint16_t i = 0; i < j["reducer_bands"].array_items().size(); ++i) {
                    band_reducers.push_back(std::make_pair(j["reducer_bands"][i][0].string_value(), j["reducer_bands"][i][1].string_value()));
                }
                return window_time_cube::create(instance()->create_from_json(j["in_cube"]), band_reducers,
                                                j["win_size_l"].int_value(), j["win_size_r"].int_value());
            }
        }));
    cube_generators.insert(std::make_pair<std::string, std::function<std::shared_ptr<cube>(json11::Json&)>>(
        "window_space", [](json11::Json& j) {
            if (!j["kernel"].is_null()) {
                std::vector<double> kernel;
                for (uint16_t i = 0; i < j["kernel"].array_items().size(); ++i) {
                    kernel.push_back(j["kernel"][i].number_value());
                }
                return window_space_cube::create(instance()->create_from_json(j["in_cube"]), kernel,
                                                j["win_size_y"].int_value(), j["win_size_x"].int_value(), 
                                                j["keep_bands"].bool_value(), j["pad_str"].string_value(), j["pad_fill"].number_value());
            }
            else {
                std::vector<std::pair<std::string, std::string>> band_reducers;
                for (uint16_t i = 0; i < j["reducer_bands"].array_items().size(); ++i) {
                    band_reducers.push_back(std::make_pair(j["reducer_bands"][i][0].string_value(), j["reducer_bands"][i][1].string_value()));
                }
                return window_space_cube::create(instance()->create_from_json(j["in_cube"]), band_reducers,
                                                j["win_size_y"].int_value(), j["win_size_x"].int_value(),
                                                j["keep_bands"].bool_value(), j["pad_str"].string_value(), j["pad_fill"].number_value());
            }
        }));
    cube_generators.insert(std::make_pair<std::string, std::function<std::shared_ptr<cube>(json11::Json&)>>(
        "select_bands", [](json11::Json& j) {
            std::vector<std::string> bands;
            for (uint16_t i = 0; i < j["bands"].array_items().size(); ++i) {
                bands.push_back(j["bands"][i].string_value());
            }
            auto x = select_bands_cube::create(instance()->create_from_json(j["in_cube"]), bands);
            return x;
        }));

    cube_generators.insert(std::make_pair<std::string, std::function<std::shared_ptr<cube>(json11::Json&)>>(
        "rename_bands", [](json11::Json& j) {
            std::map<std::string, std::string> band_names;
            for (auto it = j["band_names"].object_items().begin(); it != j["band_names"].object_items().end(); ++it) {
                band_names[it->first] = it->second.string_value();
            }
            auto x = rename_bands_cube::create(instance()->create_from_json(j["in_cube"]), band_names);
            return x;
        }));

    cube_generators.insert(std::make_pair<std::string, std::function<std::shared_ptr<cube>(json11::Json&)>>(
        "filter_pixel", [](json11::Json& j) {
            auto x = filter_pixel_cube::create(instance()->create_from_json(j["in_cube"]), j["predicate"].string_value());
            return x;
        }));
    cube_generators.insert(std::make_pair<std::string, std::function<std::shared_ptr<cube>(json11::Json&)>>(
        "filter_geom", [](json11::Json& j) {
            auto x = filter_geom_cube::create(instance()->create_from_json(j["in_cube"]), j["wkt"].string_value(), j["srs"].string_value());
            return x;
        }));

    cube_generators.insert(std::make_pair<std::string, std::function<std::shared_ptr<cube>(json11::Json&)>>(
        "fill_time", [](json11::Json& j) {
            auto x = fill_time_cube::create(instance()->create_from_json(j["in_cube"]), j["method"].string_value());
            return x;
        }));

    cube_generators.insert(std::make_pair<std::string, std::function<std::shared_ptr<cube>(json11::Json&)>>(
        "aggregate_time", [](json11::Json& j) {
            auto x = aggregate_time_cube::create(instance()->create_from_json(j["in_cube"]), j["dt"].string_value(), j["func"].string_value());
            return x;
        }));
    cube_generators.insert(std::make_pair<std::string, std::function<std::shared_ptr<cube>(json11::Json&)>>(
        "aggregate_space", [](json11::Json& j) {
            auto x = aggregate_space_cube::create(instance()->create_from_json(j["in_cube"]), j["dx"].number_value(), j["dy"].number_value(), j["func"].string_value());
            return x;
        }));

    cube_generators.insert(std::make_pair<std::string, std::function<std::shared_ptr<cube>(json11::Json&)>>(
        "simple_cube", [](json11::Json& j) {
            std::vector<std::string> files;
            for (uint16_t i = 0; i < j["files"].array_items().size(); ++i) {
                files.push_back(j["files"][i].string_value());
            }

            std::vector<std::string> datetime;
            for (uint16_t i = 0; i < j["datetime"].array_items().size(); ++i) {
                datetime.push_back(j["datetime"][i].string_value());
            }

            std::vector<std::string> bands;
            if (!j["bands"].is_null()) {
                for (uint16_t i = 0; i < j["bands"].array_items().size(); ++i) {
                    bands.push_back(j["bands"][i].string_value());
                }
            }

            std::vector<std::string> band_names;
            if (!j["band_names"].is_null()) {
                for (uint16_t i = 0; i < j["band_names"].array_items().size(); ++i) {
                    band_names.push_back(j["band_names"][i].string_value());
                }
            }
            double dx = -1;
            double dy = -1;
            if (!j["dx"].is_null()) {
                dx = j["dx"].number_value();
            }
            if (!j["dy"].is_null()) {
                dy = j["dy"].number_value();
            }
            auto x = simple_cube::create(files, datetime, bands, band_names, dx, dy);
            if (!j["strict"].is_null()) {
                x->set_strict(j["strict"].bool_value());
            }
            if (!j["chunk_size"].is_null()) {
                if (j["chunk_size"].array_items().size() == 3) {
                    x->set_chunk_size(j["chunk_size"][0].int_value(), j["chunk_size"][1].int_value(), j["chunk_size"][2].int_value());
                }
            }

            return x;
        }));

    cube_generators.insert(std::make_pair<std::string, std::function<std::shared_ptr<cube>(json11::Json&)>>(
        "apply_pixel", [](json11::Json& j) {
            std::vector<std::string> expr;
            for (uint16_t i = 0; i < j["expr"].array_items().size(); ++i) {
                expr.push_back(j["expr"][i].string_value());
            }

            if (!j["band_names"].is_null()) {
                std::vector<std::string> bandnames;
                for (uint16_t i = 0; i < j["band_names"].array_items().size(); ++i) {
                    bandnames.push_back(j["band_names"][i].string_value());
                }
                auto x = apply_pixel_cube::create(instance()->create_from_json(j["in_cube"]), expr, bandnames, j["keep_bands"].bool_value());
                return x;
            } else {
                auto x = apply_pixel_cube::create(instance()->create_from_json(j["in_cube"]), expr, {}, j["keep_bands"].bool_value());
                return x;
            }
        }));
    cube_generators.insert(std::make_pair<std::string, std::function<std::shared_ptr<cube>(json11::Json&)>>(
        "join_bands", [](json11::Json& j) {
            std::vector<std::shared_ptr<cube>> cubes;
            for (uint16_t i = 0; i < j["in_cubes"].array_items().size(); ++i) {
                cubes.push_back(instance()->create_from_json(j["in_cubes"][i]));
            }
            std::vector<std::string> prefixes;
            for (uint16_t i = 0; i < j["prefixes"].array_items().size(); ++i) {
                prefixes.push_back(j["prefixes"][i].string_value());
            }
            auto x = join_bands_cube::create(cubes, prefixes);
            return x;
        }));

    cube_generators.insert(std::make_pair<std::string, std::function<std::shared_ptr<cube>(json11::Json&)>>(
        "stream", [](json11::Json& j) {
            auto x = stream_cube::create(instance()->create_from_json(j["in_cube"]), j["command"].string_value());
            return x;
        }));

    cube_generators.insert(std::make_pair<std::string, std::function<std::shared_ptr<cube>(json11::Json&)>>(
        "select_time", [](json11::Json& j) {
            std::vector<std::string> t;
            for (uint16_t i = 0; i < j["t"].array_items().size(); ++i) {
                t.push_back(j["t"][i].string_value());
            }
            auto x = select_time_cube::create(instance()->create_from_json(j["in_cube"]), t);
            return x;
        }));

    cube_generators.insert(std::make_pair<std::string, std::function<std::shared_ptr<cube>(json11::Json&)>>(
        "slice_time", [](json11::Json& j) {
            auto x = slice_time_cube::create(instance()->create_from_json(j["in_cube"]), j["t"].int_value());
            return x;
        }));

    cube_generators.insert(std::make_pair<std::string, std::function<std::shared_ptr<cube>(json11::Json&)>>(
        "slice_space", [](json11::Json& j) {
            auto x = slice_space_cube::create(instance()->create_from_json(j["in_cube"]), j["ix"].int_value(), j["iy"].int_value());
            return x;
        }));

    cube_generators.insert(std::make_pair<std::string, std::function<std::shared_ptr<cube>(json11::Json&)>>(
        "crop", [](json11::Json& j) {
            auto x = crop_cube::create(instance()->create_from_json(j["in_cube"]),
                                              j["ix_min"].int_value(), j["ix_max"].int_value(),
                                              j["iy_min"].int_value(), j["iy_max"].int_value(),
                                              j["it_min"].int_value(), j["it_max"].int_value());
            return x;
        }));

    cube_generators.insert(std::make_pair<std::string, std::function<std::shared_ptr<cube>(json11::Json&)>>(
        "extract", [](json11::Json& j) {
            auto x = extract_geom::create(instance()->create_from_json(j["in_cube"]),
                                          j["ogr_dataset"].string_value(), j["time_column"].string_value(),
                                          j["ogr_layer"].string_value());
            return x;
        }));

    cube_generators.insert(std::make_pair<std::string, std::function<std::shared_ptr<cube>(json11::Json&)>>(
        "image_collection", [](json11::Json& j) {
            if (!filesystem::exists(j["file"].string_value())) {
                throw std::string("ERROR in cube_generators[\"image_collection\"](): image collection file does not exist.");
            }
            cube_view v = cube_view::read_json_string(j["view"].dump());
            auto x = image_collection_cube::create(j["file"].string_value(), v);
            x->set_chunk_size(j["chunk_size"][0].int_value(), j["chunk_size"][1].int_value(), j["chunk_size"][2].int_value());
            if (!j["strict"].is_null()) {
                x->set_strict(j["strict"].bool_value());
            }
            if (!j["mask"].is_null()) {
                if (j["mask"]["mask_type"].is_null()) {
                    GCBS_WARN("ERROR in cube_generators[\"image_collection\"](): missing mask type, mask will be ignored");
                } else {
                    std::string mask_type = j["mask"]["mask_type"].string_value();
                    std::vector<uint8_t> bits;
                    for (uint16_t i = 0; i < j["mask"]["bits"].array_items().size(); ++i) {
                        bits.push_back(j["mask"]["bits"][i].int_value());
                    }
                    if (mask_type == "value_mask") {
                        std::unordered_set<double> vals;
                        for (uint16_t i = 0; i < j["mask"]["values"].array_items().size(); ++i) {
                            vals.insert(j["mask"]["values"][i].number_value());
                        }
                        x->set_mask(j["mask_band"].string_value(), std::make_shared<value_mask>(vals, j["mask"]["invert"].bool_value(), bits));
                    } else if (mask_type == "range_mask") {
                        x->set_mask(j["mask_band"].string_value(), std::make_shared<range_mask>(j["mask"]["min"].number_value(), j["mask"]["max"].number_value(), j["mask"]["invert"].bool_value(), bits));
                    } else {
                        GCBS_WARN("ERROR in cube_generators[\"image_collection\"](): invalid mask type, mask will be ignored");
                    }
                }
            }
            return x;
        }));

    cube_generators.insert(std::make_pair<std::string, std::function<std::shared_ptr<cube>(json11::Json&)>>(
        "dummy", [](json11::Json& j) {
            cube_view v = cube_view::read_json_string(j["view"].dump());
            auto x = dummy_cube::create(v, j["nbands"].int_value(), j["fill"].number_value());
            x->set_chunk_size(j["chunk_size"][0].int_value(), j["chunk_size"][1].int_value(), j["chunk_size"][2].int_value());
            return x;
        }));
    cube_generators.insert(std::make_pair<std::string, std::function<std::shared_ptr<cube>(json11::Json&)>>(
        "empty", [](json11::Json& j) {
            cube_view v = cube_view::read_json_string(j["view"].dump());
            auto x = empty_cube::create(v, j["nbands"].int_value());
            x->set_chunk_size(j["chunk_size"][0].int_value(), j["chunk_size"][1].int_value(), j["chunk_size"][2].int_value());
            return x;
        }));

    cube_generators.insert(std::make_pair<std::string, std::function<std::shared_ptr<cube>(json11::Json&)>>(
        "stream_reduce_time", [](json11::Json& j) {
            std::vector<std::string> names;
            for (uint16_t i = 0; i < j["names"].array_items().size(); ++i) {
                names.push_back(j["names"][i].string_value());
            }
            auto x = stream_reduce_time_cube::create(instance()->create_from_json(j["in_cube"]), j["cmd"].string_value(), j["nbands"].int_value(), names);
            return x;
        }));

    cube_generators.insert(std::make_pair<std::string, std::function<std::shared_ptr<cube>(json11::Json&)>>(
        "stream_reduce_space", [](json11::Json& j) {
            std::vector<std::string> names;
            for (uint16_t i = 0; i < j["names"].array_items().size(); ++i) {
                names.push_back(j["names"][i].string_value());
            }
            auto x = stream_reduce_space_cube::create(instance()->create_from_json(j["in_cube"]), j["cmd"].string_value(), j["nbands"].int_value(), names);
            return x;
        }));

    cube_generators.insert(std::make_pair<std::string, std::function<std::shared_ptr<cube>(json11::Json&)>>(
        "stream_apply_pixel_cube", [](json11::Json& j) {  // FIXME
            std::vector<std::string> names;
            for (uint16_t i = 0; i < j["names"].array_items().size(); ++i) {
                names.push_back(j["names"][i].string_value());
            }
            auto x = stream_apply_pixel_cube::create(instance()->create_from_json(j["in_cube"]), j["cmd"].string_value(), j["nbands"].int_value(), names, j["keep_bands"].bool_value());
            return x;
        }));
    cube_generators.insert(std::make_pair<std::string, std::function<std::shared_ptr<cube>(json11::Json&)>>(
        "stream_apply_time_cube", [](json11::Json& j) {  // FIXME
            std::vector<std::string> names;
            for (uint16_t i = 0; i < j["names"].array_items().size(); ++i) {
                names.push_back(j["names"][i].string_value());
            }
            auto x = stream_apply_time_cube::create(instance()->create_from_json(j["in_cube"]), j["cmd"].string_value(), j["nbands"].int_value(), names, j["keep_bands"].bool_value());
            return x;
        }));

    cube_generators.insert(std::make_pair<std::string, std::function<std::shared_ptr<cube>(json11::Json&)>>(
        "ncdf", [](json11::Json& j) {
            bool auto_unpack = j["auto_unpack"].bool_value();
            auto x = ncdf_cube::create(j["file"].string_value(), auto_unpack);
            if (!j["chunk_size"].is_null()) {
                x->set_chunk_size(j["chunk_size"][0].int_value(), j["chunk_size"][1].int_value(), j["chunk_size"][2].int_value());
            }
            if (!j["band_selection"].is_null()) {
                std::vector<std::string> bands;
                for (uint32_t i = 0; i < j["band_selection"].array_items().size(); ++i) {
                    bands.push_back(j["band_selection"][i].string_value());
                }
                x->select_bands(bands);
            }
            return x;
        }));
}

}  // namespace gdalcubes
