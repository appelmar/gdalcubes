
/*
   Copyright 2018 Marius Appel <marius.appel@uni-muenster.de>

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.
*/

/**
 * This file contains the main entry for the command line client.
 */

#include <boost/algorithm/string.hpp>
#include <boost/program_options.hpp>
#include "build_info.h"
#include "filesystem.h"
#include "image_collection.h"
#include "image_collection_cube.h"
#include "reduce.h"
#include "stream.h"
#include "swarm.h"
#include "utils.h"

std::vector<std::string> string_list_from_text_file(std::string filename) {
    std::vector<std::string> out;

    std::string line;
    std::ifstream infile(filename);
    while (std::getline(infile, line))
        out.push_back(line);
    return out;
}

void print_usage(std::string command = "") {
    if (command == "create_collection") {
        std::cout << "Usage: gdalcubes create_collection [options] IN DEST" << std::endl;
        std::cout << std::endl;
        std::cout << "Create a new GDAL image collection (an SQLite database) from a list of GDAL Datasets (files, URLs, or other descriptors for GDALOpen()) and "
                     "a collection format definition. IN can be either a directory or a simple text file where each line is interpreted as a potential GDALDataset reference. If IN "
                     "is a directory, all containing files will be considered as potential GDALDatasets first and checked if they match the collection format afterwards. "
                  << std::endl;
        std::cout << std::endl;
        std::cout << "Options:" << std::endl;
        std::cout << "  -f, --format                  Path of the collection format description JSON file, this option is required" << std::endl;
        std::cout << "  -R, --recursive               If IN is a directory, do a recursive file listing" << std::endl;
        std::cout << "    , --noarchives              If given, do not scan within zip, tar, gz, tar.gz archive files" << std::endl;
        std::cout << "  -s, --strict                  Cancel if a single GDALDataset cannot be added to the collection. If not given, ignore failing datasets in the output collection" << std::endl;
        std::cout << "  -d, --debug                   Print debug messages" << std::endl;
        std::cout << std::endl;
    } else if (command == "info") {
        std::cout << "Usage: gdalcubes info SOURCE" << std::endl;
        std::cout << std::endl;
        std::cout << "Print information about a specified GDAL image collection (SOURCE)." << std::endl;
        std::cout << std::endl;
    } else if (command == "reduce") {
        std::cout << "Usage: gdalcubes reduce [options] SOURCE DEST" << std::endl;
        std::cout << std::endl;
        std::cout << "Reduce a given image collection (SOURCE) over time with the specified method and produce a single output image (DEST). The specified reducer "
                     "will be applied over all bands of the input collection. To convert the image collection to a cube, a data view "
                     "JSON file must be specified as -v or --view option. Depending on the collection's size and the location of their data "
                     "the reduction might be time-consuming."
                  << std::endl;
        std::cout << std::endl;
        std::cout << "Options:" << std::endl;
        std::cout << "  -v, --view               Filename of the JSON data view description, this option is required" << std::endl;
        std::cout << "  -r, --reducer            Reduction method, currently 'mean', 'median', 'min', 'max', or 'count', defaults to 'mean'" << std::endl;
        std::cout << "      --gdal-of            GDAL output format, defaults to GTiff" << std::endl;
        std::cout << "      --gdal-co            GDAL create options as 'KEY=VALUE' strings, can be passed multiple times" << std::endl;
        std::cout << "  -t, --threads            Number of threads used for parallel chunk processing, defaults to 1" << std::endl;
        std::cout << "      --swarm              Filename of a simple text file where each line points to a gdalcubes server API endpoint" << std::endl;
        std::cout << "  -d, --debug              Print debug messages" << std::endl;
        std::cout << std::endl;
    } else if (command == "stream") {
        std::cout << "Usage: gdalcubes stream [options] SOURCE DEST" << std::endl;
        std::cout << std::endl;
        std::cout << "Streams chunks of a cube (SOURCE) to stdin of an external program call such as R and stores output as DEST. THIS IS HIGHLY EXPERIMENTAL AND UNDER DEVELOPMENT."
                  << std::endl;
        std::cout << std::endl;
        std::cout << "Options:" << std::endl;
        std::cout << "      --exec               External program call, this option is required" << std::endl;
        std::cout << "  -v, --view               Filename of the JSON data view description, this option is required" << std::endl;
        std::cout << "  -c, --chunking           Chunk sizes in t,y,x dimensions as integer numbers, e.g. -c 16 256 256" << std::endl;
        std::cout << "  -r, --reducer            Reduction method, currently 'mean', 'median', 'min', or 'max', if not given, no reduction is performed on the result chunks." << std::endl;
        std::cout << "      --gdal-of            GDAL output format for optional reduction, defaults to GTiff, only relevant if -r is given" << std::endl;
        std::cout << "      --gdal-co            GDAL create options as 'KEY=VALUE' strings for optional redutction, can be passed multiple times, only relevant if -r is given" << std::endl;
        std::cout << "  -t, --threads            Number of threads used for parallel chunk processing, defaults to 1" << std::endl;
        std::cout << "      --swarm              Filename of a simple text file where each line points to a gdalcubes server API endpoint" << std::endl;
        std::cout << "  -d, --debug              Print debug messages" << std::endl;

        std::cout << std::endl;

    } else {
        std::cout << "Usage: gdalcubes command [arguments]" << std::endl;
        std::cout << "   or: gdalcubes [--help | --version]" << std::endl;
        std::cout << std::endl;
        std::cout << "Commands:" << std::endl;
        std::cout << "  info                     Print metadata of a GDAL image collection file " << std::endl;
        std::cout << "  create_collection        Create a new image collection from GDAL datasets" << std::endl;
        std::cout << "  reduce                   Reduce a GDAL cube over time to a single GDAL image" << std::endl;
        std::cout << "  stream                   Stream chunks of a GDAL cube to stdin of other programs" << std::endl;
        std::cout << std::endl;
        std::cout << "Please use 'gdalcubes command --help' for further information about command-specific arguments." << std::endl;
    }
}

int main(int argc, char* argv[]) {
    config::instance()->gdalcubes_init();
    config::instance()->set_default_progress_bar(std::make_shared<progress_simple_stdout_with_time>());

    namespace po = boost::program_options;
    // see https://stackoverflow.com/questions/15541498/how-to-implement-subcommands-using-boost-program-options

    po::options_description global_args("Global arguments");
    global_args.add_options()("help,h", "")("version", "")("debug,d", "")("command", po::value<std::string>(), "")("subargs", po::value<std::vector<std::string>>(), "");

    po::positional_options_description pos;
    pos.add("command", 1).add("subargs", -1);

    po::variables_map vm;

    try {
        po::parsed_options parsed = po::command_line_parser(argc, argv).options(global_args).positional(pos).allow_unregistered().run();
        po::store(parsed, vm);
        if (vm.count("version")) {
            version_info v = config::instance()->get_version_info();
            std::cout << "gdalcubes " << v.VERSION_MAJOR << "." << v.VERSION_MINOR << "." << v.VERSION_PATCH << " (" << v.GIT_COMMIT << ") built on " << v.BUILD_DATE << " " << v.BUILD_TIME << std::endl;
            return 0;
        }
        if (vm.count("help") && !vm.count("command")) {
            print_usage();
            return 0;
        }

        // if command and --help is given, print command-specific usage
        if (vm.count("help") && vm.count("command")) {
            print_usage(vm["command"].as<std::string>());
            return 0;
        }

        if (vm.count("debug")) {
            config::instance()->set_error_handler(error_handler::error_handler_debug);
        }

        std::string cmd = vm["command"].as<std::string>();
        if (cmd == "create_collection") {
            po::options_description cc_desc("create_collection arguments");
            cc_desc.add_options()("recursive,R", "Scan provided directory recursively")("format,f",
                                                                                        po::value<std::string>(), "")(
                "strict,s", "")("noarchives", "")("input", po::value<std::string>(), "")("output",
                                                                                         po::value<std::string>(),
                                                                                         "");

            po::positional_options_description cc_pos;
            cc_pos.add("input", 1).add("output", 1);

            try {
                std::vector<std::string> opts = po::collect_unrecognized(parsed.options, po::include_positional);
                opts.erase(opts.begin());
                po::store(po::command_line_parser(opts).options(cc_desc).positional(cc_pos).run(), vm);
            } catch (...) {
                std::cout << "ERROR in gdalcubes create_collection: invalid arguments." << std::endl;
                std::cout << cc_desc << std::endl;
                return 1;
            }

            bool scan_archives = true;
            bool recursive = false;
            bool strict = false;
            if (vm.count("recursive")) {
                recursive = true;
            }
            if (vm.count("strict")) {
                strict = true;
            }
            if (vm.count("noarchives")) {
                scan_archives = false;
            }

            std::string input = vm["input"].as<std::string>();
            std::string output = vm["output"].as<std::string>();
            std::string format = vm["format"].as<std::string>();

            std::vector<std::string> in;

            if (filesystem::is_directory(input)) {
                if (recursive) {
                    filesystem::iterate_directory_recursive(input, [&in](const std::string& p) {
                        if (filesystem::is_regular_file(p)) {
                            in.push_back(filesystem::make_absolute(p));  // TODO make absolute
                        }
                    });

                } else {
                    filesystem::iterate_directory(input, [&in](const std::string& p) {
                        if (filesystem::is_regular_file(p)) {
                            in.push_back(filesystem::make_absolute(p));  // TODO make absolute
                        }
                    });
                }
            } else if (filesystem::is_regular_file(input)) {
                in = string_list_from_text_file(input);
            } else {
                throw std::string("ERROR in gdalcubes create_collection: Invalid input, provide a text file or directory.");
            }

            if (scan_archives) {
                in = image_collection::unroll_archives(in);
            }

            collection_format f(format);
            auto ic = image_collection::create(f, in, strict);
            ic->write(output);
            std::cout << ic->to_string() << std::endl;

        } else if (cmd == "info") {
            po::options_description info_desc("info arguments");
            info_desc.add_options()("input", po::value<std::string>(), "Filename of the image collection.");

            po::positional_options_description info_pos;
            info_pos.add("input", 1);

            try {
                std::vector<std::string> opts = po::collect_unrecognized(parsed.options, po::include_positional);
                opts.erase(opts.begin());
                po::store(po::command_line_parser(opts).options(info_desc).positional(info_pos).run(), vm);
            } catch (...) {
                std::cout << "ERROR in gdalcubes info: invalid arguments." << std::endl;
                std::cout << info_desc << std::endl;
            }

            std::string input = vm["input"].as<std::string>();

            image_collection ic = image_collection(input);

            std::vector<image_collection::bands_row> bands = ic.get_bands();

            bounds_st e = ic.extent();
            std::cout << ic.to_string() << std::endl;
            std::cout << "DIMENSIONS: " << std::endl;
            std::cout << "  BANDS:       ";
            for (uint16_t i = 0; i < bands.size(); ++i) {
                std::cout << "(" << bands[i].name << ")";
                if (i == bands.size() - 1) {
                    std::cout << std::endl;
                } else {
                    std::cout << " ";
                }
            }

            std::cout << "  DATETIME:    "
                      << "(" << e.t0.to_string() << " - " << e.t1.to_string() << ")" << std::endl;
            std::cout << "  Y / LAT:     "
                      << "(" << e.s.bottom << " - " << e.s.top << ")" << std::endl;
            std::cout << "  X / LON:     "
                      << "(" << e.s.left << " - " << e.s.right << ")" << std::endl;

        } else if (cmd == "reduce") {
            po::options_description reduce_desc("reduce arguments");
            reduce_desc.add_options()("view,v", po::value<std::string>(), "Path to the JSON data view description");
            reduce_desc.add_options()("reducer,r", po::value<std::string>()->default_value("mean"), "Reduction method, currently mean, median, min, and max are implemented.");
            reduce_desc.add_options()("gdal-of", po::value<std::string>()->default_value("GTiff"), "GDAL output format, defaults to GTiff");
            reduce_desc.add_options()("gdal-co", po::value<std::vector<std::string>>(), "GDAL create options");
            reduce_desc.add_options()("input", po::value<std::string>(), "Filename of the input image collection.");
            reduce_desc.add_options()("output", po::value<std::string>(), "Filename of the output image.");
            reduce_desc.add_options()("threads,t", po::value<uint16_t>()->default_value(1), " Number of threads used for parallel chunk processing, defaults to 1");
            reduce_desc.add_options()("swarm", po::value<std::string>(), "Simple text file defining a gdalcubes_server swarm");

            po::positional_options_description reduce_pos;
            reduce_pos.add("input", 1);
            reduce_pos.add("output", 1);

            try {
                std::vector<std::string> opts = po::collect_unrecognized(parsed.options, po::include_positional);
                opts.erase(opts.begin());
                po::store(po::command_line_parser(opts).options(reduce_desc).positional(reduce_pos).run(), vm);
            } catch (...) {
                std::cout << "ERROR in gdalcubes reduce: invalid arguments." << std::endl;
                print_usage("reduce");
                return 1;
            }

            std::string input = vm["input"].as<std::string>();
            std::string output = vm["output"].as<std::string>();

            std::vector<std::string> create_options;
            if (vm.count("gdal-co") > 0) {
                create_options = vm["gdal-co"].as<std::vector<std::string>>();
            }
            std::string reducer = vm["reducer"].as<std::string>();
            std::string outformat = vm["gdal-of"].as<std::string>();
            std::string json_view_path = vm["view"].as<std::string>();

            uint16_t nthreads = vm["threads"].as<uint16_t>();

            if (vm.count("swarm")) {
                auto p = gdalcubes_swarm::from_txtfile(vm["swarm"].as<std::string>());
                p->set_threads(nthreads);
                config::instance()->set_default_chunk_processor(p);
            } else {
                if (nthreads > 1) {
                    config::instance()->set_default_chunk_processor(std::dynamic_pointer_cast<chunk_processor>(std::make_shared<chunk_processor_multithread>(nthreads)));
                }
            }

            std::shared_ptr<image_collection> ic = std::make_shared<image_collection>(input);
            auto c_in = image_collection_cube::create(ic, json_view_path);
            auto c_reduce = reduce_cube::create(c_in, reducer);
            c_reduce->write_gdal_image(output, outformat, create_options);

        } else if (cmd == "stream") {
            po::options_description stream_desc("stream arguments");
            stream_desc.add_options()("exec", po::value<std::string>(), "External program call for each chunk");
            stream_desc.add_options()("view,v", po::value<std::string>(), "Path to the JSON data view description");
            stream_desc.add_options()("chunking,c", po::value<std::string>(), "Chunk sizes in order t,y,x");
            stream_desc.add_options()("reducer,r", po::value<std::string>(), "Reduction method, currently mean, median, min, and max are implemented.");
            stream_desc.add_options()("gdal-of", po::value<std::string>(), "GDAL output format");
            stream_desc.add_options()("gdal-co", po::value<std::vector<std::string>>(), "GDAL create options");
            stream_desc.add_options()("input", po::value<std::string>(), "Filename of the input image collection.");
            stream_desc.add_options()("output", po::value<std::string>(), "Output file / directory.");
            stream_desc.add_options()("threads,t", po::value<uint16_t>()->default_value(1), " Number of threads used for parallel chunk processing, defaults to 1");

            po::positional_options_description stream_pos;
            stream_pos.add("input", 1);
            stream_pos.add("output", 1);

            try {
                std::vector<std::string> opts = po::collect_unrecognized(parsed.options, po::include_positional);
                opts.erase(opts.begin());
                po::store(po::command_line_parser(opts).options(stream_desc).positional(stream_pos).run(), vm);
            } catch (...) {
                print_usage("stream");
                return 1;
            }

            std::string input = vm["input"].as<std::string>();
            std::string output = vm["output"].as<std::string>();
            std::string exec = vm["exec"].as<std::string>();
            std::string json_view_path = vm["view"].as<std::string>();
            uint16_t nthreads = vm["threads"].as<uint16_t>();
            if (vm.count("swarm")) {
                auto p = gdalcubes_swarm::from_txtfile(vm["swarm"].as<std::string>());
                p->set_threads(nthreads);
                config::instance()->set_default_chunk_processor(p);
            } else {
                if (nthreads > 1) {
                    config::instance()->set_default_chunk_processor(std::dynamic_pointer_cast<chunk_processor>(std::make_shared<chunk_processor_multithread>(nthreads)));
                }
            }

            std::shared_ptr<image_collection> ic = std::make_shared<image_collection>(input);

            auto c_in = image_collection_cube::create(ic, json_view_path);

            std::vector<uint32_t> chunk_sizes;
            std::string chunkstr = vm["chunking"].as<std::string>();
            std::vector<std::string> chunkstr_tokens;
            // split chunkstr
            boost::split(chunkstr_tokens, chunkstr, boost::is_any_of(" ,;"));
            if (chunkstr_tokens.size() == 3) {
                chunk_sizes.push_back(std::stoi(chunkstr_tokens[0]));
                chunk_sizes.push_back(std::stoi(chunkstr_tokens[1]));
                chunk_sizes.push_back(std::stoi(chunkstr_tokens[2]));
            } else if (chunkstr_tokens.size() == 1) {
                if (chunkstr == "auto_temporal") {
                    chunk_sizes.push_back(c_in->size()[1]);
                    uint32_t csspatial = (uint32_t)std::ceil(std::sqrt((double)(1024 * 1024 * 8) / (double)(sizeof(double) * chunk_sizes[0] * c_in->bands().count())));  // default 8 MB
                    chunk_sizes.push_back(csspatial);
                    chunk_sizes.push_back(csspatial);
                    std::cout << "Using chunk size (t,y,x)=(" << chunk_sizes[0] << "," << chunk_sizes[1] << "," << chunk_sizes[2] << ")" << std::endl;
                } else if (chunkstr == "auto_spatial") {
                    chunk_sizes.push_back(1);
                    chunk_sizes.push_back(c_in->size()[2]);
                    chunk_sizes.push_back(c_in->size()[3]);

                    uint32_t cstemporal = (uint32_t)std::ceil((double)(1024 * 1024 * 8) / (double)(sizeof(double) * chunk_sizes[1] * chunk_sizes[2] * c_in->bands().count()));  // default 8 MB
                    chunk_sizes[0] = cstemporal;
                    std::cout << "Using chunk size (t,y,x)=(" << chunk_sizes[0] << "," << chunk_sizes[1] << "," << chunk_sizes[2] << ")" << std::endl;
                } else {
                    throw std::string("ERROR in gdalcubes stream: expected exactly three numbers, 'auto_temporal', or 'auto_spatial' as chunk size.");
                }
            } else {
                throw std::string("ERROR in gdalcubes stream: expected exactly three numbers, 'auto_temporal', or 'auto_spatial' as chunk size.");
            }

            c_in->set_chunk_size(chunk_sizes[0], chunk_sizes[1], chunk_sizes[2]);

            auto c_stream = stream_cube::create(c_in, exec);

            // TODO: do something even if no reducer is given

            if (vm.count("reducer")) {
                std::string reducer = vm["reducer"].as<std::string>();
                std::vector<std::string> create_options;
                if (vm.count("gdal-co") > 0) {
                    create_options = vm["gdal-co"].as<std::vector<std::string>>();
                }
                std::string outformat = "GTiff";
                if (vm.count("gdal-of") > 0) {
                    outformat = vm["gdal-of"].as<std::string>();
                }

                auto c_reduce = reduce_cube::create(c_stream, reducer);
                c_reduce->write_gdal_image(output, outformat, create_options);
            }

        } else {
            std::cout << "ERROR in gdalcubes: unrecognized command." << std::endl;
            print_usage();
            return 1;
        }
    } catch (std::string s) {
        std::cout << s << std::endl;
        config::instance()->gdalcubes_cleanup();
        return 1;
    } catch (std::exception& e) {
        std::cout << "ERROR in gdalcubes: unexpected exception" << std::endl;
        std::cout << "what():" << e.what() << std::endl;
        config::instance()->gdalcubes_cleanup();
        return 1;
    } catch (...) {
        std::cout << "ERROR in gdalcubes: unexpected exception" << std::endl;
        config::instance()->gdalcubes_cleanup();
        return 1;
    }
    config::instance()->gdalcubes_cleanup();
}
