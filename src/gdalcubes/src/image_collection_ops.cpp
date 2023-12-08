/*
   Copyright 2019 Marius Appel <marius.appel@hs-bochum.de>

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

#include "image_collection_ops.h"

#include <gdal_utils.h>
#include <sqlite3.h>

#include <thread>
#include <unordered_set>

#include "cube.h"

namespace gdalcubes {

void image_collection_ops::translate_gtiff(std::shared_ptr<gdalcubes::image_collection> in, std::string out_dir, uint16_t nthreads, bool overwrite, std::vector<std::string> creation_options) {
    if (!filesystem::exists(out_dir)) {
        filesystem::mkdir_recursive(out_dir);
    }

    if (!filesystem::is_directory(out_dir)) {
        throw std::string("ERROR in image_collection_ops::translate_cog(): output is not a directory.");
    }

    std::vector<std::thread> thrds;

    std::shared_ptr<progress> prg = config::instance()->get_default_progress_bar()->get();
    prg->set(0);  // explicitly set to zero to show progress bar immediately

    in->write(filesystem::join(out_dir, filesystem::filename(in->get_filename())));

    std::mutex mutex;
    std::vector<image_collection::gdalrefs_row> gdalrefs = in->get_gdalrefs();

    for (uint16_t it = 0; it < nthreads; ++it) {
        thrds.push_back(std::thread([it, nthreads, &out_dir, &gdalrefs, &prg, in, &mutex, overwrite, &creation_options]() {
            for (uint32_t i = it; i < gdalrefs.size(); i += nthreads) {
                prg->increment((double)1 / (double)gdalrefs.size());
                std::string descr = gdalrefs[i].descriptor;

                CPLStringList translate_args;

                translate_args.AddString("-of");
                translate_args.AddString("GTiff");
                for (auto it = creation_options.begin(); it != creation_options.end(); ++it) {
                    translate_args.AddString("-co");
                    translate_args.AddString(it->c_str());
                }

                translate_args.AddString("-b");
                translate_args.AddString(std::to_string(gdalrefs[i].band_num).c_str());  // band_num is 1 based

                GDALTranslateOptions* trans_options = GDALTranslateOptionsNew(translate_args.List(), NULL);
                if (trans_options == NULL) {
                    GCBS_WARN("Cannot create gdal_translate options.");
                    continue;
                }
                GDALDataset* dataset = (GDALDataset*)GDALOpen(descr.c_str(), GA_ReadOnly);
                if (!dataset) {
                    GCBS_WARN("Cannot open GDAL dataset '" + descr + "'.");
                    GDALTranslateOptionsFree(trans_options);
                    continue;
                }
                std::string outimgdir = filesystem::join(out_dir, std::to_string(gdalrefs[i].image_id));
                if (!filesystem::exists(outimgdir)) {
                    filesystem::mkdir(outimgdir);
                }
                std::string outfile = filesystem::join(outimgdir, std::to_string(gdalrefs[i].band_id) + ".tif");
                if (filesystem::exists(outfile) && !overwrite) {
                    GCBS_DEBUG(outfile + " already exists; set overwrite=true to force recreation of existing files.");
                } else {
                    GDALDatasetH out = GDALTranslate(outfile.c_str(), (GDALDatasetH)dataset, trans_options, NULL);
                    if (!out) {
                        GCBS_WARN("Cannot translate GDAL dataset '" + descr + "'.");
                        GDALClose((GDALDatasetH)dataset);
                        GDALTranslateOptionsFree(trans_options);
                    }
                    GDALClose(out);
                    GDALTranslateOptionsFree(trans_options);
                }
                GDALClose((GDALDatasetH)dataset);

                // Run SQL update anyway to fix broken links etc. if needed
                std::string sql = "UPDATE gdalrefs SET descriptor='" + outfile + "', band_num=1 " + "WHERE image_id=" + std::to_string(gdalrefs[i].image_id) + " AND band_id=" + std::to_string(gdalrefs[i].band_id) + ";";

                mutex.lock();
                if (sqlite3_exec(in->get_db_handle(), sql.c_str(), NULL, NULL, NULL) != SQLITE_OK) {
                    GCBS_WARN("Skipping image " + std::to_string(gdalrefs[i].image_id) + " due to failed band table update");
                }
                mutex.unlock();
            }
        }));
    }
    for (uint16_t it = 0; it < nthreads; ++it) {
        thrds[it].join();
    }
    prg->finalize();
}

void image_collection_ops::translate_cog(std::shared_ptr<gdalcubes::image_collection> in, std::string out_dir, uint16_t nthreads, bool overwrite, std::vector<std::string> creation_options) {
    if (!filesystem::exists(out_dir)) {
        filesystem::mkdir_recursive(out_dir);
    }

    if (!filesystem::is_directory(out_dir)) {
        throw std::string("ERROR in image_collection_ops::translate_cog(): output is not a directory.");
    }

    std::vector<std::thread> thrds;

    std::shared_ptr<progress> prg = config::instance()->get_default_progress_bar()->get();
    prg->set(0);  // explicitly set to zero to show progress bar immediately

    in->write(filesystem::join(out_dir, filesystem::filename(in->get_filename())));

    std::mutex mutex;
    std::vector<image_collection::gdalrefs_row> gdalrefs = in->get_gdalrefs();

    bool has_COG = GetGDALDriverManager()->GetDriverByName("COG") != NULL;
    if (!has_COG) {
        throw std::string("Direct translation to COG requires GDAL >= 3.1, please combine translate_gtiff and create_overviews instead");
    }

    for (uint16_t it = 0; it < nthreads; ++it) {
        thrds.push_back(std::thread([it, nthreads, &out_dir, &gdalrefs, &prg, in, &mutex, overwrite, &creation_options]() {
            for (uint32_t i = it; i < gdalrefs.size(); i += nthreads) {
                prg->increment((double)1 / (double)gdalrefs.size());
                std::string descr = gdalrefs[i].descriptor;

                CPLStringList translate_args;
                translate_args.AddString("-of");
                translate_args.AddString("COG");

                for (auto it = creation_options.begin(); it != creation_options.end(); ++it) {
                    translate_args.AddString("-co");
                    translate_args.AddString(it->c_str());
                }

                translate_args.AddString("-b");
                translate_args.AddString(std::to_string(gdalrefs[i].band_num).c_str());  // band_num is 1 based

                GDALTranslateOptions* trans_options = GDALTranslateOptionsNew(translate_args.List(), NULL);
                if (trans_options == NULL) {
                    GCBS_WARN("Cannot create gdal_translate options.");
                    continue;
                }
                GDALDataset* dataset = (GDALDataset*)GDALOpen(descr.c_str(), GA_ReadOnly);
                if (!dataset) {
                    GCBS_WARN("Cannot open GDAL dataset '" + descr + "'.");
                    GDALTranslateOptionsFree(trans_options);
                    continue;
                }
                std::string outimgdir = filesystem::join(out_dir, std::to_string(gdalrefs[i].image_id));
                if (!filesystem::exists(outimgdir)) {
                    filesystem::mkdir(outimgdir);
                }

                std::string outfile = filesystem::join(outimgdir, std::to_string(gdalrefs[i].band_id) + ".tif");
                if (filesystem::exists(outfile) && !overwrite) {
                    GCBS_DEBUG(outfile + " already exists; set overwrite=true to force recreation of existing files.");
                } else {
                    GDALDatasetH out = GDALTranslate(outfile.c_str(), (GDALDatasetH)dataset, trans_options, NULL);
                    if (!out) {
                        GCBS_WARN("Cannot translate GDAL dataset '" + descr + "'.");
                        GDALClose((GDALDatasetH)dataset);
                        GDALTranslateOptionsFree(trans_options);
                    }
                    GDALClose(out);
                    GDALTranslateOptionsFree(trans_options);
                }
                GDALClose((GDALDatasetH)dataset);

                // Run SQL update anyway to fix broken links etc. if needed
                std::string sql = "UPDATE gdalrefs SET descriptor='" + outfile + "', band_num=1 " + "WHERE image_id=" + std::to_string(gdalrefs[i].image_id) + " AND band_id=" + std::to_string(gdalrefs[i].band_id) + ";";

                mutex.lock();
                if (sqlite3_exec(in->get_db_handle(), sql.c_str(), NULL, NULL, NULL) != SQLITE_OK) {
                    GCBS_WARN("Skipping image " + std::to_string(gdalrefs[i].image_id) + " due to failed band table update");
                }
                mutex.unlock();
            }
        }));
    }
    for (uint16_t it = 0; it < nthreads; ++it) {
        thrds[it].join();
    }
    prg->finalize();
}

void image_collection_ops::create_overviews(std::shared_ptr<image_collection> in, std::vector<int> levels, std::string resampling, uint16_t nthreads) {
    std::vector<image_collection::gdalrefs_row> gdalrefs = in->get_gdalrefs();

    std::unordered_set<std::string> done;
    std::mutex m;

    std::vector<std::thread> thrds;

    std::shared_ptr<progress> prg = config::instance()->get_default_progress_bar()->get();
    prg->set(0);  // explicitly set to zero to show progress bar immediately

    for (uint16_t it = 0; it < nthreads; ++it) {
        thrds.push_back(std::thread([it, nthreads, &done, &m, &gdalrefs, &resampling, &levels, &prg]() {
            for (uint32_t i = it; i < gdalrefs.size(); i += nthreads) {
                prg->increment((double)1 / (double)gdalrefs.size());
                std::string descr = gdalrefs[i].descriptor;
                m.lock();
                if (done.count(descr) > 0) {
                    m.unlock();
                    continue;
                }
                done.insert(descr);
                m.unlock();

                GDALDataset* dataset = (GDALDataset*)GDALOpen(descr.c_str(), GA_Update);
                if (!dataset) {
                    dataset = (GDALDataset*)GDALOpen(descr.c_str(), GA_ReadOnly);
                    if (!dataset) {
                        GCBS_WARN("Cannot open GDAL dataset '" + descr + "'.");
                        continue;
                    }
                }
                if (dataset->BuildOverviews(resampling.c_str(), levels.size(), levels.data(), 0, nullptr, NULL, nullptr) == CE_Failure) {
                    GCBS_WARN("Cannot build overviews for dataset '" + descr + "'.");
                }
                GDALClose((GDALDatasetH)dataset);
            }
        }));
    }
    for (uint16_t it = 0; it < nthreads; ++it) {
        thrds[it].join();
    }
    prg->finalize();
}

}  // namespace gdalcubes
