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

/**
 * This file contains the main entry for the command line client.
 */

#include <boost/algorithm/string.hpp>
#include <boost/program_options.hpp>
#include <fstream>

#include "build_info.h"
#include "cube_factory.h"
#include "filesystem.h"
#include "image_collection.h"
#include "image_collection_cube.h"
#include "image_collection_ops.h"
#include "stream.h"
#ifndef GDALCUBES_NO_SWARM
#include "swarm.h"
#endif
#include "utils.h"

using namespace gdalcubes;

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
    } else if (command == "exec") {
        std::cout << "Usage: gdalcubes exec [options] SOURCE DEST" << std::endl;
        std::cout << std::endl;
        std::cout << "Evaluate a JSON-serialized SOURCE data cube and store the result as a NetCDF file (DEST)."
                  << std::endl;
        std::cout << std::endl;
        std::cout << "Options:" << std::endl;
        std::cout << "    , --deflate            Deflate compression level for output NetCDF file (0=no compression, 9=max compression), defaults to 1" << std::endl;
        std::cout << "  -t, --threads            Number of threads used for parallel chunk processing, defaults to 1" << std::endl;
        std::cout << "  -c, --chunk              Compute only one specific chunk, specified by its integer identifier" << std::endl;
#ifndef GDALCUBES_NO_SWARM
        std::cout << "      --swarm              Filename of a simple text file where each line points to a gdalcubes server API endpoint" << std::endl;
#endif
        std::cout << "  -d, --debug              Print debug messages" << std::endl;
        std::cout << std::endl;
    } else if (command == "addo") {
        std::cout << "Usage: gdalcubes addo [options] SOURCE" << std::endl;
        std::cout << std::endl;
        std::cout << "Build overview images for an existing image collection (SOURCE)." << std::endl;
        std::cout << std::endl;
        std::cout << "Options:" << std::endl;
        std::cout << "  -t, --threads            Number of threads used for parallel processing, defaults to 1" << std::endl;
        std::cout << "  -l, --levels             Overview levels, defaults to 2 4 8 16 32" << std::endl;
        std::cout << "  -r, --resampling         Resampling algorithm, one of \"AVERAGE\", \"AVERAGE_MAGPHASE\", \"BILINEAR\", \"CUBIC\", \"CUBICSPLINE\", \"GAUSS\", \"LANCZOS\", \"MODE\", \"NEAREST\", or \"NONE\"." << std::endl;
        std::cout << "  -d, --debug              Print debug messages" << std::endl;
        std::cout << std::endl;
    } else if (command == "translate_gtiff") {
        std::cout << "Usage: gdalcubes translate_gtiff [options] SOURCE DEST" << std::endl;
        std::cout << std::endl;
        std::cout << "Translate all images in a collection (SOURCE) to cloud-optimized GeoTiffs under the DEST directory." << std::endl;
        std::cout << std::endl;
        std::cout << "Options:" << std::endl;
        std::cout << "  -t, --threads            Number of threads used for parallel processing, defaults to 1" << std::endl;
        std::cout << "  -f, --force              Force translation to GTiff even for already existing files" << std::endl;
        std::cout << "  -d, --debug              Print debug messages" << std::endl;
        std::cout << std::endl;
    } else {
        std::cout << "Usage: gdalcubes command [arguments]" << std::endl;
        std::cout << "   or: gdalcubes [--help | --version]" << std::endl;
        std::cout << std::endl;
        std::cout << "Commands:" << std::endl;
        std::cout << "  info                     Print metadata of a GDAL image collection file " << std::endl;
        std::cout << "  create_collection        Create a new image collection from GDAL datasets" << std::endl;
        std::cout << "  exec                     Evaluate a data cube and store the result as a NetCDF file" << std::endl;
        std::cout << "  addo                     Build overview images for an existing image collection" << std::endl;
        std::cout << "  translate_gtiff          Translate all images in a collection to (tiled) GeoTiff files" << std::endl;
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
                            in.push_back(filesystem::make_absolute(p));
                        }
                    });

                } else {
                    filesystem::iterate_directory(input, [&in](const std::string& p) {
                        if (filesystem::is_regular_file(p)) {
                            in.push_back(filesystem::make_absolute(p));
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

            std::vector<image_collection::bands_row> bands = ic.get_available_bands();

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

        } else if (cmd == "exec") {
            po::options_description exec_desc("exec arguments");
            exec_desc.add_options()("input", po::value<std::string>(), "");
            exec_desc.add_options()("output", po::value<std::string>(), "");
            exec_desc.add_options()("chunk,c", po::value<uint32_t>(), "");
            exec_desc.add_options()("threads,t", po::value<uint16_t>()->default_value(1), "");
#ifndef GDALCUBES_NO_SWARM
            exec_desc.add_options()("swarm", po::value<std::string>(), "");
#endif
            exec_desc.add_options()("deflate", po::value<uint8_t>()->default_value(1), "");

            po::positional_options_description exec_pos;
            exec_pos.add("input", 1);
            exec_pos.add("output", 1);

            try {
                std::vector<std::string> opts = po::collect_unrecognized(parsed.options, po::include_positional);
                opts.erase(opts.begin());
                po::store(po::command_line_parser(opts).options(exec_desc).positional(exec_pos).run(), vm);
            } catch (...) {
                std::cout << "ERROR in gdalcubes exec: invalid arguments." << std::endl;
                print_usage("exec");
                return 1;
            }

            std::string input = vm["input"].as<std::string>();
            std::string output = vm["output"].as<std::string>();

            uint16_t nthreads = vm["threads"].as<uint16_t>();
            uint8_t deflate = vm["deflate"].as<uint8_t>();

#ifndef GDALCUBES_NO_SWARM
            if (vm.count("swarm")) {
                auto p = gdalcubes_swarm::from_txtfile(vm["swarm"].as<std::string>());
                p->set_threads(nthreads);
                config::instance()->set_default_chunk_processor(p);
            } else {
#endif
                if (nthreads > 1) {
                    config::instance()->set_default_chunk_processor(std::dynamic_pointer_cast<chunk_processor>(std::make_shared<chunk_processor_multithread>(nthreads)));
                }
#ifndef GDALCUBES_NO_SWARM
            }
#endif
            std::shared_ptr<cube> c = cube_factory::instance()->create_from_json_file(input);
            if (vm.count("chunk")) {
                c->write_single_chunk_netcdf(vm["chunk"].as<chunkid_t>(), output, deflate);
            } else {
                c->write_netcdf_file(output, deflate);
            }

        } else if (cmd == "addo") {
            po::options_description addo_desc("addo arguments");
            addo_desc.add_options()("input", po::value<std::string>(), "");
            addo_desc.add_options()("threads,t", po::value<uint16_t>()->default_value(1), "");
            addo_desc.add_options()("resampling,r", po::value<std::string>()->default_value("NEAREST"), "");
            addo_desc.add_options()("levels,l", po::value<std::vector<int>>()->multitoken()->default_value(std::vector<int>{2, 4, 8, 16, 32}, "2 4 8 16 32"), "");

            po::positional_options_description addo_pos;
            addo_pos.add("input", 1);

            try {
                std::vector<std::string> opts = po::collect_unrecognized(parsed.options, po::include_positional);
                opts.erase(opts.begin());
                po::store(po::command_line_parser(opts).options(addo_desc).positional(addo_pos).run(), vm);
            } catch (...) {
                std::cout << "ERROR in gdalcubes addo: invalid arguments." << std::endl;
                std::cout << addo_desc << std::endl;
            }

            std::string input = vm["input"].as<std::string>();

            std::shared_ptr<image_collection> ic = std::make_shared<image_collection>(input);
            uint16_t nthreads = vm["threads"].as<uint16_t>();
            std::string resampling = vm["resampling"].as<std::string>();
            std::vector<int> levels = vm["levels"].as<std::vector<int>>();

            image_collection_ops::create_overviews(ic, levels, resampling, nthreads);

        } else if (cmd == "translate_gtiff") {
            po::options_description gtiff_description("exec arguments");
            gtiff_description.add_options()("input", po::value<std::string>(), "");
            gtiff_description.add_options()("output", po::value<std::string>(), "");
            gtiff_description.add_options()("force,f", "");
            gtiff_description.add_options()("threads,t", po::value<uint16_t>()->default_value(1), "");

            po::positional_options_description cog_pos;
            cog_pos.add("input", 1);
            cog_pos.add("output", 1);

            try {
                std::vector<std::string> opts = po::collect_unrecognized(parsed.options, po::include_positional);
                opts.erase(opts.begin());
                po::store(po::command_line_parser(opts).options(gtiff_description).positional(cog_pos).run(), vm);
            } catch (...) {
                std::cout << "ERROR in gdalcubes translate_cog: invalid arguments." << std::endl;
                print_usage("exec");
                return 1;
            }

            std::string input = vm["input"].as<std::string>();
            std::string output = vm["output"].as<std::string>();

            uint16_t nthreads = vm["threads"].as<uint16_t>();
            bool force = false;
            if (vm.count("force")) {
                force = true;
            }
            std::shared_ptr<image_collection> ic = std::make_shared<image_collection>(input);

            image_collection_ops::translate_gtiff(ic, output, nthreads, force);

        }

        else {
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
