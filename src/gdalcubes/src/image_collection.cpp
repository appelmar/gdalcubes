/*
 MIT License

 Copyright (c) 2019 Marius Appel <marius.appel@hs-bochum.de>

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

#include "image_collection.h"

#include <gdalwarper.h>
#include <sqlite3.h>

#include <boost/regex.hpp>
#include <set>
#include <unordered_set>

#include "config.h"
#include "external/date.h"
#include "filesystem.h"
#include "utils.h"

namespace gdalcubes {

image_collection::image_collection() : _format(), _filename(""), _db(nullptr) {
    if (sqlite3_open_v2("", &_db, SQLITE_OPEN_READWRITE | SQLITE_OPEN_CREATE | SQLITE_OPEN_FULLMUTEX, NULL) != SQLITE_OK) {
        std::string msg = "ERROR in image_collection::create(): cannot create temporary image collection file.";
        throw msg;
    }

    // Enable foreign key constraints
    sqlite3_db_config(_db, SQLITE_DBCONFIG_ENABLE_FKEY, 1, NULL);

    // Create tables

    // key value metadata for collection
    std::string sql_schema_collection_md = "CREATE TABLE collection_md(key TEXT PRIMARY KEY, value TEXT);";
    if (sqlite3_exec(_db, sql_schema_collection_md.c_str(), NULL, NULL, NULL) != SQLITE_OK) {
        throw std::string("ERROR in image_collection::create(): cannot create image collection schema (i).");
    }

    std::string sql_insert_gdalcubes_version = "INSERT INTO collection_md(key, value) VALUES('GDALCUBES_VERSION','" + std::to_string(GDALCUBES_VERSION_MAJOR) + "." + std::to_string(GDALCUBES_VERSION_MINOR) + "." + std::to_string(GDALCUBES_VERSION_PATCH) + "');";
    if (sqlite3_exec(_db, sql_insert_gdalcubes_version.c_str(), NULL, NULL, NULL) != SQLITE_OK) {
        throw std::string("ERROR in image_collection::create(): cannot insert gdalcubes version to database.");
    }

    // Create bands table
    std::string sql_schema_bands = "CREATE TABLE bands (id INTEGER PRIMARY KEY, name TEXT, type VARCHAR(16), offset NUMERIC DEFAULT 0.0, scale NUMERIC DEFAULT 1.0, unit VARCHAR(16) DEFAULT '', nodata VARCHAR(16) DEFAULT '');";
    if (sqlite3_exec(_db, sql_schema_bands.c_str(), NULL, NULL, NULL) != SQLITE_OK) {
        throw std::string("ERROR in image_collection::create(): cannot create image collection schema (ii).");
    }

    // Create band metadata table
    std::string sql_schema_band_md = "CREATE TABLE band_md(band_id INTEGER, key TEXT, value TEXT, PRIMARY KEY (band_id, key), FOREIGN KEY (band_id) REFERENCES bands(id) ON DELETE CASCADE);";
    if (sqlite3_exec(_db, sql_schema_band_md.c_str(), NULL, NULL, NULL) != SQLITE_OK) {
        throw std::string("ERROR in image_collection::create(): cannot create image collection schema (iii).");
    }

    // Create image table
    std::string sql_schema_images = "CREATE TABLE images (id INTEGER PRIMARY KEY, name TEXT, left NUMERIC, top NUMERIC, bottom NUMERIC, right NUMERIC, datetime TEXT, proj TEXT, UNIQUE(name));CREATE INDEX idx_image_names ON images(name);";
    if (sqlite3_exec(_db, sql_schema_images.c_str(), NULL, NULL, NULL) != SQLITE_OK) {
        throw std::string("ERROR in image_collection::create(): cannot create image collection schema (iv).");
    }

    // Create image metadata table
    std::string sql_schema_image_md = "CREATE TABLE image_md(image_id INTEGER, key TEXT, value TEXT, PRIMARY KEY (image_id, key), FOREIGN KEY (image_id) REFERENCES images(id) ON DELETE CASCADE);";
    if (sqlite3_exec(_db, sql_schema_image_md.c_str(), NULL, NULL, NULL) != SQLITE_OK) {
        throw std::string("ERROR in image_collection::create(): cannot create image collection schema (v).");
    }

    // create gdal references table
    std::string sql_schema_gdalrefs =
        "CREATE TABLE gdalrefs (image_id INTEGER, band_id INTEGER, descriptor TEXT, band_num INTEGER, FOREIGN KEY (image_id) REFERENCES images(id) ON DELETE CASCADE, PRIMARY KEY (image_id, band_id), FOREIGN KEY (band_id) REFERENCES bands(id) ON DELETE CASCADE);"
        "CREATE INDEX idx_gdalrefs_bandid ON gdalrefs(band_id);"
        "CREATE INDEX idx_gdalrefs_imageid ON gdalrefs(image_id);";
    if (sqlite3_exec(_db, sql_schema_gdalrefs.c_str(), NULL, NULL, NULL) != SQLITE_OK) {
        throw std::string("ERROR in collection_format::apply(): cannot create image collection schema (vi).");
    }
}

image_collection::image_collection(collection_format format) : image_collection() {
    _format = format;
    if (format.json()["bands"].is_null()) {
        throw std::string("ERROR in image_collection::create(): image collection format does not contain any bands.");
    }

    if (_format.json()["bands"].object_items().size() == 0) {
        throw std::string("ERROR in image_collection::create(): image collection format does not contain any bands.");
    }

    std::string sql_insert_format = "INSERT INTO collection_md(key, value) VALUES('collection_format','" + sqlite_escape_singlequotes(_format.json().dump()) + "');";
    if (sqlite3_exec(_db, sql_insert_format.c_str(), NULL, NULL, NULL) != SQLITE_OK) {
        throw std::string("ERROR in image_collection::create(): cannot insert collection format to database.");
    }

    uint16_t band_id = 0;
    for (auto it = _format.json()["bands"].object_items().begin(); it != _format.json()["bands"].object_items().end(); ++it) {
        std::string sql_insert_band;
        sql_insert_band = "INSERT INTO bands(id, name";
        if (!it->second["nodata"].is_null()) sql_insert_band += ",nodata";
        if (!it->second["offset"].is_null()) sql_insert_band += ",offset";
        if (!it->second["scale"].is_null()) sql_insert_band += ",scale";
        if (!it->second["unit"].is_null()) sql_insert_band += ",unit";
        sql_insert_band += ") VALUES(" + std::to_string(band_id) + ",'" + sqlite_escape_singlequotes(it->first) + "'";
        if (!it->second["nodata"].is_null()) sql_insert_band += ",'" + std::to_string(it->second["nodata"].number_value()) + "'";
        if (!it->second["offset"].is_null()) sql_insert_band += "," + std::to_string(it->second["offset"].number_value()) + "";
        if (!it->second["scale"].is_null()) sql_insert_band += "," + std::to_string(it->second["scale"].number_value()) + "";
        if (!it->second["unit"].is_null()) sql_insert_band += ",'" + sqlite_escape_singlequotes(it->second["unit"].string_value()) + "'";
        sql_insert_band += ");";

        ++band_id;
        if (sqlite3_exec(_db, sql_insert_band.c_str(), NULL, NULL, NULL) != SQLITE_OK) {
            throw std::string("ERROR in collection_format::apply(): cannot insert bands to collection database.");
        }
    }
}

image_collection::image_collection(std::string filename) : _format(), _filename(filename), _db(nullptr) {
    // TODO: IMPLEMENT VERSIONING OF COLLECTION FORMATS AND CHECK COMPATIBILITY HERE
    if (!filesystem::exists(filename)) {
        throw std::string("ERROR in image_collection::image_collection(): input collection '" + filename + "' does not exist.");
    }
    if (sqlite3_open_v2(filename.c_str(), &_db, SQLITE_OPEN_READWRITE | SQLITE_OPEN_CREATE | SQLITE_OPEN_FULLMUTEX, NULL) != SQLITE_OK) {
        std::string msg = "ERROR in image_collection::image_collection(): cannot open existing image collection file.";
        throw msg;
    }
    // Enable foreign key constraints
    sqlite3_db_config(_db, SQLITE_DBCONFIG_ENABLE_FKEY, 1, NULL);

    // load format from database
    std::string sql_select_format = "SELECT value FROM \"collection_md\" WHERE key='collection_format';";
    sqlite3_stmt* stmt;
    sqlite3_prepare_v2(_db, sql_select_format.c_str(), -1, &stmt, NULL);
    if (!stmt) {
        GCBS_DEBUG("Failed to extract collection format from existing image collection file");
    }
    if (sqlite3_step(stmt) == SQLITE_DONE) {
        std::string msg = "No collection format from existing image collection file.";
        GCBS_DEBUG("No collection format from existing image collection file");
    } else {
        _format.load_string(sqlite_as_string(stmt, 0));
    }
    sqlite3_finalize(stmt);
}

image_collection::~image_collection() {
    if (_db) {
        sqlite3_close(_db);
        _db = nullptr;
    }
}

std::shared_ptr<image_collection> image_collection::create(collection_format format, std::vector<std::string> descriptors, bool strict) {
    std::shared_ptr<image_collection> o = std::make_shared<image_collection>(format);
    o->add_with_collection_format(descriptors, strict);
    return o;
}

std::shared_ptr<image_collection> image_collection::create(std::vector<std::string> descriptors,
                                                           std::vector<std::string> date_time,
                                                           std::vector<std::string> band_names,
                                                           bool use_subdatasets, bool one_band_per_file) {
    std::shared_ptr<image_collection> o = std::make_shared<image_collection>();
    if (one_band_per_file) {
        o->add_with_datetime_bands(descriptors, date_time, band_names, use_subdatasets);
    }
    else {
        o->add_with_datetime(descriptors, date_time, band_names, use_subdatasets);
    }

    return o;
}

std::shared_ptr<image_collection> image_collection::create() {
    std::shared_ptr<image_collection> o = std::make_shared<image_collection>();
    return o;
}

struct image_band {
    GDALDataType type;
    std::string unit;
    double scale;
    double offset;
    std::string nodata;
};

void image_collection::add_with_datetime(std::vector<std::string> descriptors, std::vector<std::string> date_time,
                                         std::vector<std::string> band_names, bool use_subdatasets) {
    if (!_format.is_null()) {
        GCBS_WARN("Image collection has nonempty format; trying to apply the format to provided datasets");
        add_with_collection_format(descriptors);
        return;
    }

    if (descriptors.size() != date_time.size()) {
        GCBS_ERROR("The number of provided datasets must be identical to the number of provided date/time strings");
        throw std::string("The number of provided datasets must be identical to the number of provided date/time strings");
    }

    if (use_subdatasets) {
        std::vector<std::string> subdatasets;
        for (auto it = descriptors.begin(); it != descriptors.end(); ++it) {
            GDALDataset* dataset = (GDALDataset*)GDALOpen((*it).c_str(), GA_ReadOnly);
            if (!dataset) {
                GCBS_WARN("GDAL failed to open " + *it);
                continue;
            }

            // Is there a SUBDATASETS metadata domain?
            char** md_domains = dataset->GetMetadataDomainList();
            if (md_domains != NULL) {
                if (CSLFindString(md_domains, "SUBDATASETS") != -1) {
                    // if yes, list all metadata keys ending with _NAME
                    char** md_sd = dataset->GetMetadata("SUBDATASETS");
                    if (md_sd != NULL) {
                        for (uint16_t imd = 0; imd < CSLCount(md_sd); ++imd) {
                            std::string s(md_sd[imd]);
                            size_t ii = s.find("_NAME=");
                            if (ii != std::string::npos) {
                                // found
                                subdatasets.push_back(s.substr(ii + 6));
                            }
                        }
                        // Don't call CSLDestroy(md_sd);
                    }
                }
                CSLDestroy(md_domains);
            }
            GDALClose((GDALDatasetH)dataset);
        }
        descriptors = subdatasets;  // TODO: how to handle input datasets if they do not have any subdatasets?
    }

    std::string sql_select_bands = "SELECT id, name, type, offset, scale, unit FROM bands";

    std::vector<image_band> bands;
    std::vector<std::string> band_names_db;
    std::vector<uint32_t> band_ids;

    sqlite3_stmt* stmt;
    sqlite3_prepare_v2(_db, sql_select_bands.c_str(), -1, &stmt, NULL);
    while (sqlite3_step(stmt) == SQLITE_ROW) {
        image_band b;
        b.type = utils::gdal_type_from_string(sqlite_as_string(stmt, 2));
        b.offset = sqlite3_column_double(stmt, 3);
        b.scale = sqlite3_column_double(stmt, 4);
        b.unit = sqlite_as_string(stmt, 5);
        bands.push_back(b);
        band_ids.push_back(sqlite3_column_int(stmt, 0));
        band_names_db.push_back(sqlite_as_string(stmt, 1));
    }
    sqlite3_finalize(stmt);
    bool collection_contains_bands = !bands.empty();

    if (collection_contains_bands) {
        if (!band_names.empty()) {
            if (band_names.size() != band_names_db.size()) {
                GCBS_ERROR("Provided band names are in conflict with existing bands in the collection");
                throw std::string("Provided band names are in conflict with existing bands in the collection");
            }

            // compare names
            for (uint16_t ib = 0; ib < band_names.size(); ++ib) {
                if (band_names[ib] != band_names_db[ib]) {
                    GCBS_ERROR("Provided band names are in conflict with existing bands in the collection");
                    throw std::string("Provided band names are in conflict with existing bands in the collection");
                }
            }
        } else {
            band_names = band_names_db;
        }
        // now, we can be sure that if collection already has bands, band_names and band_names_db are identical
    }

    std::shared_ptr<progress> p = config::instance()->get_default_progress_bar()->get();
    p->set(0);  // explicitly set to zero to show progress bar immediately
    for (uint32_t i = 0; i < descriptors.size(); ++i) {
        GDALDataset* dataset = (GDALDataset*)GDALOpen(descriptors[i].c_str(), GA_ReadOnly);
        if (!dataset) {
            GCBS_WARN("GDAL failed to open '" + descriptors[i] + "'; dataset will be skipped");
            continue;
        }

        double affine_in[6] = {0, 0, 1, 0, 0, 1};
        bounds_2d<double> bbox;
        std::string srs_str;
        if (dataset->GetGeoTransform(affine_in) != CE_None) {
            GCBS_WARN("GDAL failed to fetch geotransform parameters for '" + descriptors[i] + "'; dataset will be skipped");
            GDALClose(dataset);
            continue;
        }

        bbox.left = affine_in[0];
        bbox.right = affine_in[0] + affine_in[1] * dataset->GetRasterXSize() + affine_in[2] * dataset->GetRasterYSize();
        bbox.top = affine_in[3];
        bbox.bottom = affine_in[3] + affine_in[4] * dataset->GetRasterXSize() + affine_in[5] * dataset->GetRasterYSize();
        OGRSpatialReference srs_in;

        srs_in.SetFromUserInput(dataset->GetProjectionRef());
        if ( srs_in.GetAuthorityName(NULL) != NULL &&  srs_in.GetAuthorityCode(NULL) != NULL) {
            srs_str = std::string( srs_in.GetAuthorityName(NULL)) + ":" + std::string(srs_in.GetAuthorityCode(NULL));
        }
        else {
            char *tmp;
            srs_in.exportToWkt(&tmp);
            srs_str = std::string(tmp);
            CPLFree(tmp);
        }

        bbox.transform(srs_str, "EPSG:4326");

        if (!collection_contains_bands) {
            // first dataset, bands table is empty
            if (!band_names.empty()) {
                if (dataset->GetRasterCount() != (int)band_names.size()) {
                    std::string msg = "Got " + std::to_string(band_names.size()) + " names but image '" + descriptors[i] + "' has " + std::to_string(dataset->GetRasterCount()) + " bands; please make sure that numbers of names and bands of all datasets are compatible";
                    GCBS_ERROR(msg);
                    GDALClose(dataset);
                    throw msg;
                }
            } else {
                for (uint16_t ib = 0; ib < dataset->GetRasterCount(); ++ib) {
                    band_names.push_back("band" + std::to_string(ib + 1));
                }
            }
            // now, we can be sure that band_names is not empty

            for (uint16_t ib = 0; ib < dataset->GetRasterCount(); ++ib) {
                image_band b;
                b.type = dataset->GetRasterBand(ib + 1)->GetRasterDataType();
                b.offset = dataset->GetRasterBand(ib + 1)->GetOffset();
                b.scale = dataset->GetRasterBand(ib + 1)->GetScale();
                b.unit = dataset->GetRasterBand(ib + 1)->GetUnitType();
                b.nodata = "";
                int hasnodata = 0;
                double nd = dataset->GetRasterBand(ib + 1)->GetNoDataValue(&hasnodata);
                if (hasnodata)
                    b.nodata = std::to_string(nd);
                bands.push_back(b);

                std::string sql_band_insert = "INSERT INTO bands(name, type, offset, scale, unit, nodata) VALUES ('" +
                                              sqlite_escape_singlequotes(band_names[ib]) + "', '" + utils::string_from_gdal_type(b.type) + "'," + std::to_string(b.offset) + "," +
                                              std::to_string(b.scale) + ",'" + sqlite_escape_singlequotes(b.unit) + "','" + b.nodata + "')";
                if (sqlite3_exec(_db, sql_band_insert.c_str(), NULL, NULL, NULL) != SQLITE_OK) {
                    GCBS_ERROR("Failed to insert band into image collection database");
                    GDALClose(dataset);
                    throw std::string("Failed to insert band into image collection database");
                }
                band_ids.push_back(sqlite3_last_insert_rowid(_db));
                collection_contains_bands = true;
            }
        } else {
            // check consistency with other images
            bool is_compatible = true;
            if (dataset->GetRasterCount() != (int)bands.size()) {
                is_compatible = false;
            }
            for (uint16_t ib = 0; ib < dataset->GetRasterCount(); ++ib) {
                is_compatible = (bands[ib].type == dataset->GetRasterBand(ib + 1)->GetRasterDataType()) &&
                                (bands[ib].offset == dataset->GetRasterBand(ib + 1)->GetOffset()) &&
                                (bands[ib].scale == dataset->GetRasterBand(ib + 1)->GetScale()) &&
                                (bands[ib].unit == dataset->GetRasterBand(ib + 1)->GetUnitType());
            }
            if (!is_compatible) {
                GCBS_WARN("Bands of image '" + descriptors[i] + "' are not identical to bands of other images in the collection; dataset will be skipped");
                GDALClose(dataset);
                continue;
            }
        }

        datetime d = datetime::from_string(date_time[i]);

        // add image to database
        std::string image_name = descriptors[i];
        sqlite3_exec(_db, "BEGIN TRANSACTION;", NULL, NULL, NULL);  // what if this fails?!

        // TODO: change from srs_str to WKT
        std::string sql_insert_image = "INSERT INTO images(name, datetime, left, top, bottom, right, proj) VALUES('" + sqlite_escape_singlequotes(image_name) + "','" +
                                       d.to_string() + "'," +
                                       std::to_string(bbox.left) + "," + std::to_string(bbox.top) + "," + std::to_string(bbox.bottom) + "," + std::to_string(bbox.right) + ",'" + sqlite_escape_singlequotes(srs_str) + "')";
        if (sqlite3_exec(_db, sql_insert_image.c_str(), NULL, NULL, NULL) != SQLITE_OK) {
            GCBS_WARN("Failed to add image '" + descriptors[i] + "'; dataset will be skipped");
            GDALClose(dataset);
            sqlite3_exec(_db, "ROLLBACK TRANSACTION;", NULL, NULL, NULL);  // what if this fails?!
            continue;
        }
        uint32_t image_id = sqlite3_last_insert_rowid(_db);  // take care of race conditions if things run parallel at some point

        // add gdalrefs (one for each band) to database
        for (uint16_t ib = 0; ib < bands.size(); ++ib) {
            std::string sql_insert_gdalref = "INSERT INTO gdalrefs(descriptor, image_id, band_id, band_num) VALUES('" + sqlite_escape_singlequotes(descriptors[i]) + "'," + std::to_string(image_id) + "," + std::to_string(band_ids[ib]) + "," + std::to_string(ib + 1) + ");";
            if (sqlite3_exec(_db, sql_insert_gdalref.c_str(), NULL, NULL, NULL) != SQLITE_OK) {
                GCBS_WARN("Failed to add image '" + descriptors[i] + "'; dataset will be skipped");
                sqlite3_exec(_db, "ROLLBACK TRANSACTION;", NULL, NULL, NULL);  // what if this fails?!
                break;
            }
        }
        sqlite3_exec(_db, "COMMIT TRANSACTION;", NULL, NULL, NULL);  // what if this fails?!
        GDALClose(dataset);
        p->increment(double(1) / double(descriptors.size()));
    }

    p->set(1);
    p->finalize();
}


void image_collection::add_with_datetime_bands(std::vector<std::string> descriptors, std::vector<std::string> date_time,
                                         std::vector<std::string> band_names, bool use_subdatasets) {
    if (!_format.is_null()) {
        GCBS_WARN("Image collection has nonempty format; trying to apply the format to provided datasets");
        add_with_collection_format(descriptors);
        return;
    }

    if (descriptors.size() != date_time.size()) {
        GCBS_ERROR("The number of provided datasets must be identical to the number of provided date/time strings");
        throw std::string("The number of provided datasets must be identical to the number of provided date/time strings");
    }
    if (descriptors.size() != band_names.size()) {
        GCBS_ERROR("The number of provided datasets must be identical to the number of provided band names");
        throw std::string("The number of provided datasets must be identical to the number of provided band names");
    }

    if (use_subdatasets) {
        std::vector<std::string> subdatasets;
        for (auto it = descriptors.begin(); it != descriptors.end(); ++it) {
            GDALDataset* dataset = (GDALDataset*)GDALOpen((*it).c_str(), GA_ReadOnly);
            if (!dataset) {
                GCBS_WARN("GDAL failed to open " + *it);
                continue;
            }

            // Is there a SUBDATASETS metadata domain?
            char** md_domains = dataset->GetMetadataDomainList();
            if (md_domains != NULL) {
                if (CSLFindString(md_domains, "SUBDATASETS") != -1) {
                    // if yes, list all metadata keys ending with _NAME
                    char** md_sd = dataset->GetMetadata("SUBDATASETS");
                    if (md_sd != NULL) {
                        for (uint16_t imd = 0; imd < CSLCount(md_sd); ++imd) {
                            std::string s(md_sd[imd]);
                            size_t ii = s.find("_NAME=");
                            if (ii != std::string::npos) {
                                // found
                                subdatasets.push_back(s.substr(ii + 6));
                            }
                        }
                        // Don't call CSLDestroy(md_sd);
                    }
                }
                CSLDestroy(md_domains);
            }
            GDALClose((GDALDatasetH)dataset);
        }
        descriptors = subdatasets;  // TODO: how to handle input datasets if they do not have any subdatasets?
    }

    std::string sql_select_bands = "SELECT id, name, type, offset, scale, unit FROM bands";

    std::map<std::string, image_band> bands;
    std::map<std::string, uint32_t> band_ids;

    sqlite3_stmt* stmt;
    sqlite3_prepare_v2(_db, sql_select_bands.c_str(), -1, &stmt, NULL);
    while (sqlite3_step(stmt) == SQLITE_ROW) {
        image_band b;
        b.type = utils::gdal_type_from_string(sqlite_as_string(stmt, 2));
        b.offset = sqlite3_column_double(stmt, 3);
        b.scale = sqlite3_column_double(stmt, 4);
        b.unit = sqlite_as_string(stmt, 5);
        std::string name = sqlite_as_string(stmt, 1);
        uint32_t id = sqlite3_column_int(stmt, 0);
        bands.insert(std::make_pair(name, b));
        band_ids.insert(std::make_pair(name, id));
    }
    sqlite3_finalize(stmt);

    std::shared_ptr<progress> p = config::instance()->get_default_progress_bar()->get();
    p->set(0);  // explicitly set to zero to show progress bar immediately
    for (uint32_t i = 0; i < descriptors.size(); ++i) {
        GDALDataset* dataset = (GDALDataset*)GDALOpen(descriptors[i].c_str(), GA_ReadOnly);
        if (!dataset) {
            GCBS_WARN("GDAL failed to open '" + descriptors[i] + "'; dataset will be skipped");
            continue;
        }

        double affine_in[6] = {0, 0, 1, 0, 0, 1};
        bounds_2d<double> bbox;
        std::string srs_str;
        if (dataset->GetGeoTransform(affine_in) != CE_None) {
            GCBS_WARN("GDAL failed to fetch geotransform parameters for '" + descriptors[i] + "'; dataset will be skipped");
            GDALClose(dataset);
            continue;
        }

        bbox.left = affine_in[0];
        bbox.right = affine_in[0] + affine_in[1] * dataset->GetRasterXSize() + affine_in[2] * dataset->GetRasterYSize();
        bbox.top = affine_in[3];
        bbox.bottom = affine_in[3] + affine_in[4] * dataset->GetRasterXSize() + affine_in[5] * dataset->GetRasterYSize();
        OGRSpatialReference srs_in;

        srs_in.SetFromUserInput(dataset->GetProjectionRef());
        if ( srs_in.GetAuthorityName(NULL) != NULL &&  srs_in.GetAuthorityCode(NULL) != NULL) {
            srs_str = std::string( srs_in.GetAuthorityName(NULL)) + ":" + std::string(srs_in.GetAuthorityCode(NULL));
        }
        else {
            char *tmp;
            srs_in.exportToWkt(&tmp);
            srs_str = std::string(tmp);
            CPLFree(tmp);
        }

        bbox.transform(srs_str, "EPSG:4326");

        if (bands.find(band_names[i]) == bands.end()) {
            if (dataset->GetRasterCount() > 1) {
                GCBS_WARN("Dataset '" + descriptors[i] + " has > 1 bands, only band 1 will be considered");
            }
            image_band b;
            b.type = dataset->GetRasterBand(1)->GetRasterDataType();
            b.offset = dataset->GetRasterBand(1)->GetOffset();
            b.scale = dataset->GetRasterBand(1)->GetScale();
            b.unit = dataset->GetRasterBand(1)->GetUnitType();
            b.nodata = "";
            int hasnodata = 0;
            double nd = dataset->GetRasterBand(1)->GetNoDataValue(&hasnodata);
            if (hasnodata)
                b.nodata = std::to_string(nd);
            std::string name = band_names[i];
            bands.insert(std::make_pair(name, b));
            std::string sql_band_insert = "INSERT INTO bands(name, type, offset, scale, unit, nodata) VALUES ('" +
                                          sqlite_escape_singlequotes(name) + "', '" + utils::string_from_gdal_type(b.type) + "'," + std::to_string(b.offset) + "," +
                                          std::to_string(b.scale) + ",'" + sqlite_escape_singlequotes(b.unit) + "','" + b.nodata + "')";
            if (sqlite3_exec(_db, sql_band_insert.c_str(), NULL, NULL, NULL) != SQLITE_OK) {
                GCBS_ERROR("Failed to insert band into image collection database");
                GDALClose(dataset);
                throw std::string("Failed to insert band into image collection database");
            }
            band_ids.insert(std::make_pair(name, sqlite3_last_insert_rowid(_db)));
        }
        else { // band already exists
            bool is_compatible = true;
            if (dataset->GetRasterCount() > 1) {
                GCBS_WARN("Dataset '" + descriptors[i] + " has > 1 bands, only band 1 will be considered");
            }
            std::string name = band_names[i];
            is_compatible = (bands[name].type == dataset->GetRasterBand(1)->GetRasterDataType()) &&
                            (bands[name].offset == dataset->GetRasterBand(1)->GetOffset()) &&
                            (bands[name].scale == dataset->GetRasterBand(1)->GetScale()) &&
                            (bands[name].unit == dataset->GetRasterBand(1)->GetUnitType());

            if (!is_compatible) {
                GCBS_WARN("Band " + name + " of image '" + descriptors[i] + "' is not identical to the same band of other images in the collection; dataset will be skipped");
                GDALClose(dataset);
                continue;
            }
        }

        datetime d = datetime::from_string(date_time[i]);

        sqlite3_exec(_db, "BEGIN TRANSACTION;", NULL, NULL, NULL);

        // TODO: change from srs_str to WKT
        // Note: This results in a separate rows in the images table even if the files contain different bands of the same image
        std::string image_name = filesystem::stem(descriptors[i]);
        std::string sql_insert_image = "INSERT INTO images(name, datetime, left, top, bottom, right, proj) VALUES('" + sqlite_escape_singlequotes(image_name) + "','" +
                                       d.to_string() + "'," +
                                       std::to_string(bbox.left) + "," + std::to_string(bbox.top) + "," + std::to_string(bbox.bottom) + "," + std::to_string(bbox.right) + ",'" + srs_str + "')";
        if (sqlite3_exec(_db, sql_insert_image.c_str(), NULL, NULL, NULL) != SQLITE_OK) {
            GCBS_WARN("Failed to add image '" + descriptors[i] + "'; dataset will be skipped");
            GDALClose(dataset);
            sqlite3_exec(_db, "ROLLBACK TRANSACTION;", NULL, NULL, NULL);  // what if this fails?!
            continue;
        }
        uint32_t image_id = sqlite3_last_insert_rowid(_db);

        // add to gdalrefs table
        std::string sql_insert_gdalref = "INSERT INTO gdalrefs(descriptor, image_id, band_id, band_num) VALUES('" + sqlite_escape_singlequotes(descriptors[i]) + "'," + std::to_string(image_id) + "," + std::to_string(band_ids[band_names[i]]) + "," + std::to_string(1) + ");";
        if (sqlite3_exec(_db, sql_insert_gdalref.c_str(), NULL, NULL, NULL) != SQLITE_OK) {
            GCBS_WARN("Failed to add '" + descriptors[i] + "'; dataset will be skipped");
            sqlite3_exec(_db, "ROLLBACK TRANSACTION;", NULL, NULL, NULL);
            continue;
        }
        sqlite3_exec(_db, "COMMIT TRANSACTION;", NULL, NULL, NULL);

        GDALClose(dataset);
        p->increment(double(1) / double(descriptors.size()));
    }

    p->set(1);
    p->finalize();
}

void image_collection::add_with_collection_format(std::vector<std::string> descriptors, bool strict) {
    std::vector<boost::regex> regex_band_pattern;

    if (_format.is_null()) {
        GCBS_ERROR("Failed to add GDAL dataset to image collection due to missing collection format entry");
        throw std::string("Failed to add GDAL dataset to image collection due to missing collection format entry");
    }

    /* TODO: The following will fail if other applications create image collections and assign ids to bands differently.
     * A better solution would be to load band ids, names, and nums from the database bands table directly
     */

    std::vector<std::string> band_name;
    std::vector<uint16_t> band_num;
    std::vector<uint16_t> band_ids;
    std::vector<bool> band_complete;

    uint16_t band_id = 0;
    for (auto it = _format.json()["bands"].object_items().begin(); it != _format.json()["bands"].object_items().end(); ++it) {
        band_name.push_back(it->first);
        regex_band_pattern.push_back(boost::regex(it->second["pattern"].string_value()));
        if (!it->second["band"].is_null()) {
            band_num.push_back(it->second["band"].int_value());
        } else {
            band_num.push_back(1);
        }
        band_ids.push_back(band_id);
        band_complete.push_back(false);
        ++band_id;
    }

    std::string global_pattern = "";
    if (!_format.json()["pattern"].is_null()) global_pattern = _format.json()["pattern"].string_value();
    boost::regex regex_global_pattern(global_pattern);

    if (_format.json()["images"]["pattern"].is_null())
        throw std::string("ERROR in image_collection::add(): image collection format does not contain a composition rule for images.");
    boost::regex regex_images(_format.json()["images"]["pattern"].string_value());

    // @TODO: Make datetime optional, e.g., for DEMs

    std::string datetime_format = "%Y-%m-%d";
    if (_format.json()["datetime"].is_null() || _format.json()["datetime"]["pattern"].is_null()) {
        throw std::string("ERROR in image_collection::add(): image collection format does not contain a rule to derive date/time.");
    }

    boost::regex regex_datetime(_format.json()["datetime"]["pattern"].string_value());
    if (!_format.json()["datetime"]["format"].is_null()) {
        datetime_format = _format.json()["datetime"]["format"].string_value();
    }

    bool time_as_bands = false;  // time is stored as bands in the datasets
    duration band_time_delta;
    if (!_format.json()["datetime"]["bands"].is_null()) {
        // check if any band has band_num != 1
        for (uint32_t ib = 0; ib < band_num.size(); ++ib) {
            if (band_num[ib] != 1) {
                GCBS_ERROR("Collection format has conflicting declaration of bands; one dataset may either accept multiple bands, or multiple points in time, but not both");
                throw std::string("Collection format has conflicting declaration of bands; one dataset may either accept multiple bands, or multiple points in time, but not both");
            }
        }

        if (!_format.json()["datetime"]["bands"]["dt"].is_null()) {
            time_as_bands = true;
            band_time_delta = duration::from_string(_format.json()["datetime"]["bands"]["dt"].string_value());
        } else {
            GCBS_WARN("Collection format seems to have time information in dataset bands, bot does not define how to relate bands to time; time extraction might be incorrect");
        }
    }

    bool use_subdatasets = false;
    if (!_format.json()["subdatasets"].is_null()) {
        use_subdatasets = _format.json()["subdatasets"].bool_value();
    }

    if (use_subdatasets) {
        std::vector<std::string> subdatasets;
        for (auto it = descriptors.begin(); it != descriptors.end(); ++it) {
            GDALDataset* dataset = (GDALDataset*)GDALOpen((*it).c_str(), GA_ReadOnly);
            if (!dataset) {
                if (strict) throw std::string("ERROR in image_collection::add(): GDAL cannot open '" + *it + "'.");
                GCBS_WARN("GDAL failed to open " + *it);
                continue;
            }

            // Is there a SUBDATASETS metadata domain?
            char** md_domains = dataset->GetMetadataDomainList();
            if (md_domains != NULL) {
                if (CSLFindString(md_domains, "SUBDATASETS") != -1) {
                    // if yes, list all metadata keys ending with _NAME
                    char** md_sd = dataset->GetMetadata("SUBDATASETS");
                    if (md_sd != NULL) {
                        for (uint16_t imd = 0; imd < CSLCount(md_sd); ++imd) {
                            std::string s(md_sd[imd]);
                            size_t ii = s.find("_NAME=");
                            if (ii != std::string::npos) {
                                // found
                                subdatasets.push_back(s.substr(ii + 6));
                            }
                        }
                        // Don't call CSLDestroy(md_sd);
                    }
                }
                CSLDestroy(md_domains);
            }
            GDALClose((GDALDatasetH)dataset);
        }
        descriptors = subdatasets;  // TODO: how to handle input datasets if they do not have any subdatasets?
    }

    std::string global_srs_str = "";
    OGRSpatialReference global_srs;
    if (!_format.json()["srs"].is_null()) {
        global_srs_str = _format.json()["srs"].string_value();
        if (global_srs.SetFromUserInput(global_srs_str.c_str()) != OGRERR_NONE) {
            GCBS_WARN("Cannot read global SRS definition in collection format, trying to extract from individual datasets.");
            global_srs_str = "";
        }
    }

    uint32_t counter = -1;
    std::shared_ptr<progress> p = config::instance()->get_default_progress_bar()->get();
    p->set(0);  // explicitly set to zero to show progress bar immediately
    for (auto it = descriptors.begin(); it != descriptors.end(); ++it) {
        ++counter;

        p->set((double)counter / (double)descriptors.size());
        if (!global_pattern.empty()) {  // prevent unnecessary GDALOpen calls
            if (!boost::regex_match(*it, regex_global_pattern)) {
                // std::cout << "ignoring " << *it << std::endl;
                GCBS_DEBUG("Dataset " + *it + " doesn't match the global collection pattern and will be ignored");
                continue;
            }
        }

        // Read GDAL metadata
        GDALDataset* dataset = (GDALDataset*)GDALOpen((*it).c_str(), GA_ReadOnly);
        if (!dataset) {
            if (strict) throw std::string("ERROR in image_collection::add(): GDAL cannot open '" + *it + "'.");
            GCBS_WARN("GDAL failed to open " + *it);
            continue;
        }
        // if check = false, the following is not really needed if image is already in the database due to another file.
        double affine_in[6] = {0, 0, 1, 0, 0, 1};
        bounds_2d<double> bbox;
        std::string srs_str;
        if (dataset->GetGeoTransform(affine_in) != CE_None) {
            // No affine transformation, maybe GCPs?
            if (dataset->GetGCPCount() > 0) {
                dataset->GetGCPProjection();  // TODO replace with GetGCPSpatialRef if available (GDAL > 2.5)

                // First try, find GCPs for corner pixels
                double xmin = std::numeric_limits<double>::infinity();
                double ymin = std::numeric_limits<double>::infinity();
                double xmax = -std::numeric_limits<double>::infinity();
                double ymax = -std::numeric_limits<double>::infinity();
                bool x1 = false, x2 = false, x3 = false, x4 = false;

                for (int32_t igcp = 0; igcp < dataset->GetGCPCount(); ++igcp) {
                    GDAL_GCP gcp = dataset->GetGCPs()[igcp];
                    if (gcp.dfGCPLine == 0 && gcp.dfGCPPixel == 0) {
                        x1 = true;
                        if (gcp.dfGCPX < xmin) xmin = gcp.dfGCPX;
                        if (gcp.dfGCPX > xmax) xmax = gcp.dfGCPX;
                        if (gcp.dfGCPY < ymin) ymin = gcp.dfGCPY;
                        if (gcp.dfGCPY > ymax) ymax = gcp.dfGCPY;
                    } else if (gcp.dfGCPLine == dataset->GetRasterYSize() - 1 && gcp.dfGCPPixel == 0) {
                        x2 = true;
                        if (gcp.dfGCPX < xmin) xmin = gcp.dfGCPX;
                        if (gcp.dfGCPX > xmax) xmax = gcp.dfGCPX;
                        if (gcp.dfGCPY < ymin) ymin = gcp.dfGCPY;
                        if (gcp.dfGCPY > ymax) ymax = gcp.dfGCPY;
                    } else if (gcp.dfGCPLine == 0 && gcp.dfGCPPixel == dataset->GetRasterXSize() - 1) {
                        x3 = true;
                        if (gcp.dfGCPX < xmin) xmin = gcp.dfGCPX;
                        if (gcp.dfGCPX > xmax) xmax = gcp.dfGCPX;
                        if (gcp.dfGCPY < ymin) ymin = gcp.dfGCPY;
                        if (gcp.dfGCPY > ymax) ymax = gcp.dfGCPY;
                    } else if (gcp.dfGCPLine == dataset->GetRasterYSize() - 1 && gcp.dfGCPPixel == dataset->GetRasterXSize() - 1) {
                        x4 = true;
                        if (gcp.dfGCPX < xmin) xmin = gcp.dfGCPX;
                        if (gcp.dfGCPX > xmax) xmax = gcp.dfGCPX;
                        if (gcp.dfGCPY < ymin) ymin = gcp.dfGCPY;
                        if (gcp.dfGCPY > ymax) ymax = gcp.dfGCPY;
                    }
                }

                if (x1 && x2 && x3 && x4) {
                    // use extent from corner GCPS
                    bbox.left = xmin;
                    bbox.right = xmax;
                    bbox.top = ymax;
                    bbox.bottom = ymin;
                    OGRSpatialReference srs_in;
                    srs_in.SetFromUserInput(dataset->GetGCPProjection());

                    if ( srs_in.GetAuthorityName(NULL) != NULL &&  srs_in.GetAuthorityCode(NULL) != NULL) {
                        srs_str = std::string( srs_in.GetAuthorityName(NULL)) + ":" + std::string(srs_in.GetAuthorityCode(NULL));
                    }
                    else {
                        char *tmp;
                        srs_in.exportToWkt(&tmp);
                        srs_str = std::string(tmp);
                        CPLFree(tmp);
                    }
                    bbox.transform(srs_str, "EPSG:4326");
                } else {
                    //approximate extent based on gdalwarp
                    double approx_geo_transform[6];
                    int nx = 0, ny = 0;
                    double extent[4];
                    OGRSpatialReference srs_in;
                    srs_in.SetFromUserInput(dataset->GetGCPProjection());
                    if ( srs_in.GetAuthorityName(NULL) != NULL &&  srs_in.GetAuthorityCode(NULL) != NULL) {
                        srs_str = std::string( srs_in.GetAuthorityName(NULL)) + ":" + std::string(srs_in.GetAuthorityCode(NULL));
                    }
                    else {
                        char *tmp;
                        srs_in.exportToWkt(&tmp);
                        srs_str = std::string(tmp);
                        CPLFree(tmp);
                    }

                    CPLStringList transform_args;
                    transform_args.AddString(("SRC_SRS=" + std::string(srs_str)).c_str());
                    transform_args.AddString("DST_SRS=EPSG:4326");
                    transform_args.AddString("GCPS_OK=TRUE");

                    // TODO: add further transformation options if needed
                    void* transform = GDALCreateGenImgProjTransformer2(dataset, NULL, transform_args.List());
                    if (GDALSuggestedWarpOutput2(dataset,
                                                 GDALGenImgProjTransform, transform,
                                                 approx_geo_transform, &nx, &ny, extent, 0) != CE_None) {
                        if (strict) throw std::string("ERROR in image_collection::add(): GDAL cannot derive extent for '" + *it + "'.");
                        GCBS_WARN("Failed to derive spatial extent from " + *it);
                    }

                    // TODO: error handling
                    bbox.left = extent[0];
                    bbox.right = extent[2];
                    bbox.top = extent[3];
                    bbox.bottom = extent[1];
                }

            }
            else {
                char **slist = dataset->GetMetadata("GEOLOCATION");
                if (slist != NULL) {
                    std::string x_dataset = std::string(CSLFetchNameValue(slist, "X_DATASET"));
                    std::string y_dataset = std::string(CSLFetchNameValue(slist, "Y_DATASET"));
                    std::string geoloc_srs_str = std::string(CSLFetchNameValue(slist, "SRS"));
                    int16_t x_band = std::stoi(std::string(CSLFetchNameValue(slist, "X_BAND")));
                    int16_t y_band = std::stoi(std::string(CSLFetchNameValue(slist, "Y_BAND")));

                    double adfMinMax[2];
                    int bGotMin, bGotMax;

                    GDALDataset* gd_x = (GDALDataset*)GDALOpen(x_dataset.c_str(), GA_ReadOnly);
                    if (!gd_x) {
                        GDALClose((GDALDatasetH)dataset);
                        if (strict) throw std::string("ERROR in image_collection::add(): GDAL cannot open '" + x_dataset + "'.");
                        GCBS_WARN("GDAL failed to open " + x_dataset);
                        continue;
                    }
                    GDALRasterBand* b_x = gd_x->GetRasterBand(x_band);
                    adfMinMax[0] = b_x->GetMinimum( &bGotMin );
                    adfMinMax[1] = b_x->GetMaximum( &bGotMax );
                    if(!(bGotMin && bGotMax))
                        GDALComputeRasterMinMax((GDALRasterBandH)b_x, FALSE, adfMinMax);
                    bbox.left = adfMinMax[0];
                    bbox.right = adfMinMax[1];
                    GDALClose(gd_x);

                    GDALDataset* gd_y = (GDALDataset*)GDALOpen(y_dataset.c_str(), GA_ReadOnly);
                    if (!gd_y) {
                        GDALClose((GDALDatasetH)dataset);
                        if (strict) throw std::string("ERROR in image_collection::add(): GDAL cannot open '" + y_dataset + "'.");
                        GCBS_WARN("GDAL failed to open " + y_dataset);
                        continue;
                    }
                    GDALRasterBand* b_y = gd_y->GetRasterBand(y_band);
                    adfMinMax[0] = b_y->GetMinimum( &bGotMin );
                    adfMinMax[1] = b_y->GetMaximum( &bGotMax );
                    if(!(bGotMin && bGotMax))
                        GDALComputeRasterMinMax((GDALRasterBandH)b_y, FALSE, adfMinMax);
                    bbox.bottom = adfMinMax[0];
                    bbox.top = adfMinMax[1];
                    GDALClose(gd_y);

                    bbox.transform(geoloc_srs_str, "EPSG:4326");
                }
                else { // No extent??? 
                    GDALClose((GDALDatasetH)dataset);
                    if (strict) throw std::string("ERROR in image_collection::add(): GDAL cannot derive spatial extent for '" + *it + "'.");
                    GCBS_WARN("Failed to derive spatial extent from " + *it);
                    continue;
                }
            
             }  
        } else {
            bbox.left = affine_in[0];
            bbox.right = affine_in[0] + affine_in[1] * dataset->GetRasterXSize() + affine_in[2] * dataset->GetRasterYSize();
            bbox.top = affine_in[3];
            bbox.bottom = affine_in[3] + affine_in[4] * dataset->GetRasterXSize() + affine_in[5] * dataset->GetRasterYSize();
            OGRSpatialReference srs_in;

            if (global_srs_str.empty()) {  // if no global SRS is given
                srs_in.SetFromUserInput(dataset->GetProjectionRef());
            } else {
                if (dataset->GetProjectionRef() != NULL && !std::string(dataset->GetProjectionRef()).empty()) {
                    srs_in.SetFromUserInput(dataset->GetProjectionRef());
                    if (!srs_in.IsSame(&global_srs)) {
                        GCBS_WARN("SRS of dataset '" + (*it) + "' is different from global SRS and will be overwritten.");
                    }
                }
                srs_in = global_srs;
            }

            if ( srs_in.GetAuthorityName(NULL) != NULL &&  srs_in.GetAuthorityCode(NULL) != NULL) {
                srs_str = std::string( srs_in.GetAuthorityName(NULL)) + ":" + std::string(srs_in.GetAuthorityCode(NULL));
            }
            else {
                char *tmp;
                srs_in.exportToWkt(&tmp);
                srs_str = std::string(tmp);
                CPLFree(tmp);
            }
            bbox.transform(srs_str, "EPSG:4326");
        }

        // TODO: check consistency for all files of an image?!
        // -> add parameter checks=true / false

        boost::cmatch res_image;
        if (!boost::regex_match(it->c_str(), res_image, regex_images)) {
            if (strict) throw std::string("ERROR in image_collection::add(): image composition rule failed for " + std::string(*it));
            GCBS_WARN("Skipping " + *it + " due to failed image composition rule");
            continue;
        }

        if (!time_as_bands) {
            // Input dataset is a SINGLE image with only one point in time

            std::vector<image_band> bands;
            for (uint16_t i = 0; i < dataset->GetRasterCount(); ++i) {
                image_band b;
                b.type = dataset->GetRasterBand(i + 1)->GetRasterDataType();
                b.offset = dataset->GetRasterBand(i + 1)->GetOffset();
                b.scale = dataset->GetRasterBand(i + 1)->GetScale();
                b.unit = dataset->GetRasterBand(i + 1)->GetUnitType();
                b.nodata = "";
                int hasnodata = 0;
                double nd = dataset->GetRasterBand(i + 1)->GetNoDataValue(&hasnodata);
                if (hasnodata)
                    b.nodata = std::to_string(nd);
                bands.push_back(b);
            }
            if (bands.empty()) {
                if (strict) throw std::string("ERROR in image_collection::add(): " + *it + " doesn't contain any band data and will be ignored");
                GCBS_WARN("Dataset " + *it + " doesn't contain any band data and will be ignored");
                continue;
            }

            uint32_t image_id;
            std::string sql_select_image = "SELECT id FROM images WHERE name='" + res_image[1].str() + "'";
            sqlite3_stmt* stmt;
            sqlite3_prepare_v2(_db, sql_select_image.c_str(), -1, &stmt, NULL);
            if (!stmt) {
                // @TODO do error check here
            }
            if (sqlite3_step(stmt) == SQLITE_DONE) {
                // Empty result --> image has not been added before

                // @TODO: Shall we check that all files óf the same image have the same date / time? Currently we don't.

                // Extract datetime
                boost::cmatch res_datetime;
                if (!boost::regex_match(it->c_str(), res_datetime, regex_datetime)) {  // not sure to continue or throw an exception here...
                    if (strict) throw std::string("ERROR in image_collection::add(): datetime rule failed for " + std::string(*it));
                    GCBS_WARN("Skipping " + *it + " due to failed datetime rule");
                    continue;
                }

                std::stringstream os;
                date::sys_seconds pt;
                pt = datetime::tryparse(datetime_format, res_datetime[1].str());
                os << date::format("%Y-%m-%dT%H:%M:%S", pt);

                // Convert to ISO string including separators (boost::to_iso_string or boost::to_iso_extended_string do not work with SQLite datetime functions)
                std::string sql_insert_image = "INSERT OR IGNORE INTO images(name, datetime, left, top, bottom, right, proj) VALUES('" + res_image[1].str() + "','" +
                                               os.str() + "'," +
                                               std::to_string(bbox.left) + "," + std::to_string(bbox.top) + "," + std::to_string(bbox.bottom) + "," + std::to_string(bbox.right) + ",'" + srs_str + "')";
                if (sqlite3_exec(_db, sql_insert_image.c_str(), NULL, NULL, NULL) != SQLITE_OK) {
                    if (strict) throw std::string("ERROR in image_collection::add(): cannot add image to images table.");
                    GCBS_WARN("Skipping " + *it + " due to failed image table insert");
                    continue;
                }
                image_id = sqlite3_last_insert_rowid(_db);  // take care of race conditions if things run parallel at some point

            } else {
                image_id = sqlite3_column_int(stmt, 0);
                // TODO: if checks, compare l,r,b,t, datetime,srs_str from images table with current GDAL dataset
            }
            sqlite3_finalize(stmt);

            // Insert into gdalrefs table
            for (uint16_t i = 0; i < band_name.size(); ++i) {
                if (boost::regex_match(*it, regex_band_pattern[i])) {
                    // TODO: if checks, check whether bandnum exists in GDALdataset
                    // TODO: if checks, compare band type, offset, scale, unit, etc. with current GDAL dataset

                    if (!band_complete[i]) {
                        std::string sql_band_update = "UPDATE bands SET type='" + utils::string_from_gdal_type(bands[band_num[i] - 1].type) + "'";

                        if (_format.json()["bands"][band_name[i]]["scale"].is_null())
                            sql_band_update += ",scale=" + std::to_string(bands[band_num[i] - 1].scale);
                        if (_format.json()["bands"][band_name[i]]["offset"].is_null())
                            sql_band_update += ",offset=" + std::to_string(bands[band_num[i] - 1].offset);
                        if (_format.json()["bands"][band_name[i]]["unit"].is_null())
                            sql_band_update += ",unit='" + bands[band_num[i] - 1].unit + "'";

                        // TODO: also add no data if not defined in image collection?
                        sql_band_update += " WHERE name='" + band_name[i] + "';";

                        if (sqlite3_exec(_db, sql_band_update.c_str(), NULL, NULL, NULL) != SQLITE_OK) {
                            if (strict) throw std::string("ERROR in image_collection::add(): cannot update band table.");
                            GCBS_WARN("Skipping " + *it + " due to failed band table update");
                            continue;
                        }
                        band_complete[i] = true;
                    }

                    std::string sql_insert_gdalref = "INSERT INTO gdalrefs(descriptor, image_id, band_id, band_num) VALUES('" + sqlite_escape_singlequotes(*it) + "'," + std::to_string(image_id) + "," + std::to_string(band_ids[i]) + "," + std::to_string(band_num[i]) + ");";
                    if (sqlite3_exec(_db, sql_insert_gdalref.c_str(), NULL, NULL, NULL) != SQLITE_OK) {
                        if (strict) throw std::string("ERROR in image_collection::add(): cannot add dataset to gdalrefs table.");
                        GCBS_WARN("Skipping " + *it + "  due to failed gdalrefs insert");
                        break;  // break only works because there is nothing after the loop.
                    }
                }
            }

            // Read image metadata from GDALDataset
            std::unordered_set<std::string> image_md_fields;

            if (!_format.json()["image_md_fields"].is_null()) {
                for (uint16_t imd_fields = 0; imd_fields < _format.json()["image_md_fields"].array_items().size(); ++imd_fields) {
                    image_md_fields.insert(_format.json()["image_md_fields"][imd_fields].string_value());
                }
            }

            if (image_md_fields.size() > 0) {
                char** md_domains = dataset->GetMetadataDomainList();
                for (auto cur_md_key = image_md_fields.begin(); cur_md_key != image_md_fields.end(); ++cur_md_key) {
                    // has domain?
                    std::size_t sep_pos = cur_md_key->find_first_of(":");
                    if (sep_pos != std::string::npos) {
                        // has domain

                        std::string domain = cur_md_key->substr(0, sep_pos);
                        std::string field = cur_md_key->substr(sep_pos + 1, std::string::npos);

                        // does the domain exist?
                        if (CSLFindString(md_domains, domain.c_str()) == -1) {
                            // no
                            continue;
                        } else {
                            // yes
                            const char* value = CSLFetchNameValue(dataset->GetMetadata(domain.c_str()), field.c_str());
                            if (value) {
                                std::string sql_insert_image_md = "INSERT OR IGNORE INTO image_md(image_id, key, value) VALUES(" + std::to_string(image_id) + ",'" + *cur_md_key + "','" + std::string(value) + "');";
                                sqlite3_exec(_db, sql_insert_image_md.c_str(), NULL, NULL, NULL);
                            }
                        }
                    } else {
                        // default domain
                        const char* value = CSLFetchNameValue(dataset->GetMetadata(), cur_md_key->c_str());
                        if (value) {
                            std::string sql_insert_image_md = "INSERT OR IGNORE INTO image_md(image_id, key, value) VALUES(" + std::to_string(image_id) + ",'" + *cur_md_key + "','" + std::string(value) + "');";
                            sqlite3_exec(_db, sql_insert_image_md.c_str(), NULL, NULL, NULL);
                        }
                    }
                }
                CSLDestroy(md_domains);
            }

        } else {
            // Input dataset is multitemporal, bands represent different points in time
            // Add as multiple images to the image collection as

            boost::cmatch res_datetime;
            if (!boost::regex_match(it->c_str(), res_datetime, regex_datetime)) {  // not sure to continue or throw an exception here...
                if (strict) throw std::string("ERROR in image_collection::add(): datetime rule failed for " + std::string(*it));
                GCBS_WARN("Skipping " + *it + " due to failed datetime rule");
                continue;
            }

            date::sys_seconds pt;
            pt = datetime::tryparse(datetime_format, res_datetime[1].str());

            // find the corresponding band of the dataset (there can be only 1 because bands represent time)
            // and update band information in database if needed
            int16_t band_index = -1;
            for (uint16_t i = 0; i < band_name.size(); ++i) {
                if (boost::regex_match(*it, regex_band_pattern[i])) {
                    band_index = i;
                    break;
                }
            }
            if (band_index == -1) {
                continue;
            }

            if (!band_complete[band_index]) {
                image_band b;
                b.type = dataset->GetRasterBand(1)->GetRasterDataType();
                b.offset = dataset->GetRasterBand(1)->GetOffset();
                b.scale = dataset->GetRasterBand(1)->GetScale();
                b.unit = dataset->GetRasterBand(1)->GetUnitType();
                b.nodata = "";
                int hasnodata = 0;
                double nd = dataset->GetRasterBand(1)->GetNoDataValue(&hasnodata);
                if (hasnodata)
                    b.nodata = std::to_string(nd);
                std::string sql_band_update = "UPDATE bands SET type='" + utils::string_from_gdal_type(b.type) + "'";

                if (_format.json()["bands"][band_name[band_index]]["scale"].is_null())
                    sql_band_update += ",scale=" + std::to_string(b.scale);
                if (_format.json()["bands"][band_name[band_index]]["offset"].is_null())
                    sql_band_update += ",offset=" + std::to_string(b.offset);
                if (_format.json()["bands"][band_name[band_index]]["unit"].is_null())
                    sql_band_update += ",unit='" + b.unit + "'";

                // TODO: also add no data if not defined in image collection?
                sql_band_update += " WHERE name='" + band_name[band_index] + "';";

                if (sqlite3_exec(_db, sql_band_update.c_str(), NULL, NULL, NULL) != SQLITE_OK) {
                    if (strict) throw std::string("ERROR in image_collection::add(): cannot update band table.");
                    GCBS_WARN("Skipping " + *it + " due to failed band table update");
                    continue;
                }
                band_complete[band_index] = true;
            }

            // for all time steps (bands in the current dataset)
            for (uint16_t i = 0; i < dataset->GetRasterCount(); ++i) {
                // derive datetime
                datetime t = datetime(pt, band_time_delta.dt_unit) + (band_time_delta * i);

                // add image to collection
                std::string image_name = res_image[1].str() + "_" + t.to_string();
                std::string sql_insert_image = "INSERT OR IGNORE INTO images(name, datetime, left, top, bottom, right, proj) VALUES('" + image_name + "','" +
                                               t.to_string(datetime_unit::SECOND) + "'," +
                                               std::to_string(bbox.left) + "," + std::to_string(bbox.top) + "," + std::to_string(bbox.bottom) + "," + std::to_string(bbox.right) + ",'" + srs_str + "')";
                if (sqlite3_exec(_db, sql_insert_image.c_str(), NULL, NULL, NULL) != SQLITE_OK) {
                    if (strict) throw std::string("ERROR in image_collection::add(): cannot add image to images table.");
                    GCBS_WARN("Skipping " + *it + " due to failed image table insert");
                    continue;
                }

                uint32_t image_id = sqlite3_last_insert_rowid(_db);  // take care of race conditions if things run parallel at some point

                // add gdalref to collection
                std::string sql_insert_gdalref = "INSERT INTO gdalrefs(descriptor, image_id, band_id, band_num) VALUES('" + sqlite_escape_singlequotes(*it) + "'," + std::to_string(image_id) + "," + std::to_string(band_ids[band_index]) + "," + std::to_string(band_num[band_index]) + ");";
                if (sqlite3_exec(_db, sql_insert_gdalref.c_str(), NULL, NULL, NULL) != SQLITE_OK) {
                    if (strict) throw std::string("ERROR in image_collection::add(): cannot add dataset to gdalrefs table.");
                    GCBS_WARN("Skipping " + *it + "  due to failed gdalrefs insert");
                    break;  // break only works because there is nothing after the loop.
                }
            }
        }
        GDALClose((GDALDatasetH)dataset);
    }
    p->set(1);
    p->finalize();
}

void image_collection::add_with_collection_format(std::string descriptor, bool strict) {
    std::vector<std::string> x{descriptor};
    return add_with_collection_format(x, strict);
}

void image_collection::write(const std::string filename) {
    if (_filename.compare(filename) == 0) {
        // nothing to do
        return;
    }

    if (!_db) {
        throw std::string("ERROR in image_collection::write(): database handle is not open");
    }
    _filename = filename;

    // Store DB
    sqlite3_backup* db_backup;
    sqlite3* out_db;

    if (sqlite3_open_v2(_filename.c_str(), &out_db, SQLITE_OPEN_READWRITE | SQLITE_OPEN_CREATE | SQLITE_OPEN_FULLMUTEX, NULL) != SQLITE_OK) {
        throw std::string("ERROR in image_collection::write(): cannot create output database file.");
    }
    db_backup = sqlite3_backup_init(out_db, "main", _db, "main");
    if (!db_backup) {
        throw std::string("ERROR in image_collection::write(): cannot create temporary database backup object.");
    }
    sqlite3_backup_step(db_backup, -1);
    sqlite3_backup_finish(db_backup);

    sqlite3_close(_db);
    _db = out_db;

    // Enable foreign key constraints
    sqlite3_db_config(_db, SQLITE_DBCONFIG_ENABLE_FKEY, 1, NULL);  // this is important!
}

uint16_t image_collection::count_bands() {
    std::string sql = "SELECT COUNT(*) FROM bands;";

    sqlite3_stmt* stmt;
    sqlite3_prepare_v2(_db, sql.c_str(), -1, &stmt, NULL);
    if (!stmt) {
        throw std::string("ERROR in image_collection::count_bands(): cannot read query result");
    }
    sqlite3_step(stmt);
    uint16_t out = sqlite3_column_int(stmt, 0);
    sqlite3_finalize(stmt);
    return out;
}

uint32_t image_collection::count_images() {
    std::string sql = "SELECT COUNT(*) FROM images;";

    sqlite3_stmt* stmt;
    sqlite3_prepare_v2(_db, sql.c_str(), -1, &stmt, NULL);
    if (!stmt) {
        throw std::string("ERROR in image_collection::count_images(): cannot read query result");
    }
    sqlite3_step(stmt);
    uint16_t out = sqlite3_column_int(stmt, 0);
    sqlite3_finalize(stmt);
    return out;
}

uint32_t image_collection::count_gdalrefs() {
    std::string sql = "SELECT COUNT(*) FROM gdalrefs;";

    sqlite3_stmt* stmt;
    sqlite3_prepare_v2(_db, sql.c_str(), -1, &stmt, NULL);
    if (!stmt) {
        throw std::string("ERROR in image_collection::count_gdalrefs(): cannot read query result");
    }
    sqlite3_step(stmt);
    uint16_t out = sqlite3_column_int(stmt, 0);
    sqlite3_finalize(stmt);
    return out;
}

std::string image_collection::to_string() {
    std::stringstream ss;
    ss << "IMAGE COLLECTION '" << (_filename.empty() ? "unnamed" : _filename) << "' has ";
    ss << std::to_string(count_images()) << " images with ";
    ss << std::to_string(count_bands()) << " bands from ";
    ss << std::to_string(count_gdalrefs()) << " GDAL dataset references";

    auto band_info = get_available_bands();
    ss << std::endl
       << "NAME | OFFSET | SCALE | UNIT | NODATA | IMAGE COUNT" << std::endl;
    for (uint16_t i = 0; i < band_info.size(); ++i) {
        ss << band_info[i].name << " | " << band_info[i].offset << " | " << band_info[i].scale << " | " << band_info[i].nodata << " | " << band_info[i].image_count << std::endl;
    }

    return ss.str();
}

void image_collection::filter_bands(std::vector<std::string> bands) {
    // This implementation requires a foreign key constraint for gdalrefs table with cascade delete

    if (bands.empty()) {
        throw std::string("ERROR in image_collection::filter_bands(): no bands selected");
    }
    if (bands.size() == count_bands())
        return;

    std::string bandlist;
    for (uint16_t i = 0; i < bands.size() - 1; ++i) {
        bandlist += "'" + bands[i] + "',";
    }
    bandlist += "'" + bands[bands.size() - 1] + "'";

    std::string sql = "DELETE FROM bands WHERE name NOT IN (" + bandlist + ");";
    if (sqlite3_exec(_db, sql.c_str(), NULL, NULL, NULL) != SQLITE_OK) {
        throw std::string("ERROR in image_collection::filter_bands(): cannot remove bands from collection.");
    }
}

void image_collection::filter_datetime_range(date::sys_seconds start, date::sys_seconds end) {
    // This implementation requires a foreign key constraint for the gdalrefs table with cascade delete

    std::ostringstream os;

    os << date::format("%Y-%m-%dT%H:%M:%S", start);
    std::string start_str = os.str();
    os.clear();
    os << date::format("%Y-%m-%dT%H:%M:%S", end);
    std::string end_str = os.str();

    std::string sql = "DELETE FROM images WHERE datetime(images.datetime) < datetime('" + start_str + "')  OR datetime(images.datetime) > datetime('" + end_str + "');";
    if (sqlite3_exec(_db, sql.c_str(), NULL, NULL, NULL) != SQLITE_OK) {
        throw std::string("ERROR in image_collection::filter_datetime_range(): cannot remove images from collection.");
    }
}

void image_collection::filter_spatial_range(bounds_2d<double> range, std::string proj) {
    // This implementation requires a foreign key constraint for the gdalrefs table with cascade delete

    range.transform(proj, "EPSG:4326");
    std::string sql = "DELETE FROM images WHERE images.right < " + std::to_string(range.left) + " OR images.left > " + std::to_string(range.right) + " OR images.bottom > " + std::to_string(range.top) + " OR images.top < " + std::to_string(range.bottom) + ";";

    if (sqlite3_exec(_db, sql.c_str(), NULL, NULL, NULL) != SQLITE_OK) {
        throw std::string("ERROR in image_collection::filter_spatial_range(): cannot remove images from collection.");
    }
}

uint16_t image_collection::pixel_size_bytes(std::string band) {
    std::string sql = "SELECT type FROM bands";
    if (!band.empty()) sql += " WHERE name='" + band + "'";
    sql += ";";

    sqlite3_stmt* stmt;
    sqlite3_prepare_v2(_db, sql.c_str(), -1, &stmt, NULL);
    if (!stmt) {
        throw std::string("ERROR in image_collection::pixel_size_bytes(): cannot prepare query statement");
    }
    uint16_t out = 0;
    while (sqlite3_step(stmt) == SQLITE_ROW) {
        GDALDataType type = (GDALDataType)sqlite3_column_int(stmt, 0);
        out += GDALGetDataTypeSizeBytes(type);
    }
    sqlite3_finalize(stmt);
    return out;
}

bounds_st image_collection::extent() {
    std::string sql = "SELECT min(left), max(right), min(bottom), max(top), min(datetime), max(datetime) FROM images;";
    sqlite3_stmt* stmt;
    sqlite3_prepare_v2(_db, sql.c_str(), -1, &stmt, NULL);
    if (!stmt) {
        throw std::string("ERROR in image_collection::extent(): cannot prepare query statement");
    }
    bounds_st out;
    if (sqlite3_step(stmt) == SQLITE_ROW) {
        out.s.left = sqlite3_column_double(stmt, 0);
        out.s.right = sqlite3_column_double(stmt, 1);
        out.s.bottom = sqlite3_column_double(stmt, 2);
        out.s.top = sqlite3_column_double(stmt, 3);
        out.t0 = datetime::from_string(sqlite_as_string(stmt, 4));
        out.t1 = datetime::from_string(sqlite_as_string(stmt, 5));
    } else {
        throw std::string("ERROR in image_collection::extent(): cannot fetch query results");
    }
    sqlite3_finalize(stmt);
    return out;
}

std::vector<image_collection::find_range_st_row> image_collection::find_range_st(bounds_st range, std::string srs,
                                                                                 std::vector<std::string> bands, std::vector<std::string> order_by) {
    bounds_2d<double> range_trans = (srs == "EPSG:4326") ? range.s : range.s.transform(srs, "EPSG:4326");
    std::string sql =  // TODO: do we really need image_name ?
        "SELECT gdalrefs.image_id, images.name, gdalrefs.descriptor, images.datetime, bands.name, gdalrefs.band_num, images.proj "
        "FROM images INNER JOIN gdalrefs ON images.id = gdalrefs.image_id INNER JOIN bands ON gdalrefs.band_id = bands.id WHERE "
        "strftime('%Y-%m-%dT%H:%M:%S', images.datetime) >= '" +
        range.t0.to_string(datetime_unit::SECOND) + "' AND strftime('%Y-%m-%dT%H:%M:%S', images.datetime) <= '" + range.t1.to_string(datetime_unit::SECOND) +
        "' AND NOT "
        "(images.right < " +
        std::to_string(range_trans.left) + " OR images.left > " + std::to_string(range_trans.right) + " OR images.bottom > " + std::to_string(range_trans.top) + " OR images.top < " + std::to_string(range_trans.bottom) + ")";

    if (!bands.empty()) {
        std::string bandlist = "";
        for (uint16_t i = 0; i < bands.size() - 1; ++i) {
            bandlist += "'" + bands[i] + "',";
        }
        bandlist += "'" + bands[bands.size() - 1] + "'";
        sql += " AND bands.name IN (" + bandlist + ")";
    }
    if (!order_by.empty()) {
        sql += "ORDER BY ";
        for (uint16_t io = 0; io < order_by.size() - 1; ++io) {
            if (order_by[io] == "gdalrefs.image_id" ||
                order_by[io] == "images.name" ||
                order_by[io] == "gdalrefs.descriptor" ||
                order_by[io] == "images.datetime" ||
                order_by[io] == "bands.name" ||
                order_by[io] == "images.proj" ||
                order_by[io] == "gdalrefs.band_num") {
                sql += order_by[io] + ",";
            } else {
                throw std::string("ERROR in image_collection::find_range_st(): invalid column for sorting");
            }
        }
        uint16_t io = order_by.size() - 1;
        if (order_by[io] == "gdalrefs.image_id" ||
            order_by[io] == "images.name" ||
            order_by[io] == "gdalrefs.descriptor" ||
            order_by[io] == "images.datetime" ||
            order_by[io] == "bands.name" ||
            order_by[io] == "images.proj" ||
            order_by[io] == "gdalrefs.band_num") {
            sql += order_by[io];
        } else {
            throw std::string("ERROR in image_collection::find_range_st(): invalid column for sorting");
        }
    }
    sql += ";";

    sqlite3_stmt* stmt;
    sqlite3_prepare_v2(_db, sql.c_str(), -1, &stmt, NULL);
    if (!stmt) {
        throw std::string("ERROR in image_collection::find_range_st(): cannot prepare query statement");
    }
    std::vector<find_range_st_row> out;
    while (sqlite3_step(stmt) == SQLITE_ROW) {
        find_range_st_row r;
        r.image_id = sqlite3_column_int(stmt, 0);
        r.image_name = sqlite_as_string(stmt, 1);
        r.descriptor = sqlite_as_string(stmt, 2);
        r.datetime = sqlite_as_string(stmt, 3);
        r.band_name = sqlite_as_string(stmt, 4);
        r.band_num = sqlite3_column_int(stmt, 5);
        r.srs = sqlite_as_string(stmt, 6);

        out.push_back(r);
    }
    sqlite3_finalize(stmt);
    return out;
}

std::vector<image_collection::bands_row> image_collection::get_all_bands() {
    std::vector<image_collection::bands_row> out;

    //std::string sql = "SELECT id, name, type, offset,scale, unit, nodata FROM bands ORDER BY name;";  // changing the order my have consequences to data read implementations

    // all bands with image count
    std::string sql = "SELECT id, name, type, offset,scale, unit, nodata , sum(n) FROM (SELECT bands.id, bands.name, bands.type, bands.offset, bands.scale, bands.unit, bands.nodata, count(*) as n FROM bands INNER JOIN gdalrefs ON bands.id = gdalrefs.band_id GROUP BY bands.id UNION SELECT bands.id, bands.name, bands.type, bands.offset, bands.scale, bands.unit, bands.nodata, 0 FROM bands) GROUP BY id ORDER BY name";
    sqlite3_stmt* stmt;
    sqlite3_prepare_v2(_db, sql.c_str(), -1, &stmt, NULL);
    if (!stmt) {
        throw std::string("ERROR in image_collection::get_all_bands(): cannot prepare query statement");
    }

    while (sqlite3_step(stmt) == SQLITE_ROW) {
        image_collection::bands_row row;
        row.id = sqlite3_column_int(stmt, 0);
        row.name = sqlite_as_string(stmt, 1);
        row.type = utils::gdal_type_from_string(sqlite_as_string(stmt, 2));
        row.offset = sqlite3_column_double(stmt, 3);
        row.scale = sqlite3_column_double(stmt, 4);
        row.unit = sqlite_as_string(stmt, 5);
        row.nodata = sqlite_as_string(stmt, 6);
        row.image_count = sqlite3_column_int(stmt, 7);
        out.push_back(row);
    }
    sqlite3_finalize(stmt);
    return out;
}

std::vector<image_collection::bands_row> image_collection::get_available_bands() {
    std::vector<image_collection::bands_row> out;
    std::vector<image_collection::bands_row> all_bands = get_all_bands();
    for (auto it = all_bands.begin(); it != all_bands.end(); ++it) {
        if (it->image_count > 0)
            out.push_back(*it);
    }

    return out;
}

sqlite3* image_collection::get_db_handle() {
    return _db;
}

std::vector<image_collection::gdalrefs_row> image_collection::get_gdalrefs() {
    std::vector<image_collection::gdalrefs_row> out;

    std::string sql = "SELECT image_id, band_id, descriptor, band_num FROM gdalrefs";
    sqlite3_stmt* stmt;
    sqlite3_prepare_v2(_db, sql.c_str(), -1, &stmt, NULL);
    if (!stmt) {
        throw std::string("ERROR in image_collection::get_gdalrefs(): cannot prepare query statement");
    }

    while (sqlite3_step(stmt) == SQLITE_ROW) {
        image_collection::gdalrefs_row row;
        row.image_id = sqlite3_column_int(stmt, 0);
        row.band_id = sqlite3_column_int(stmt, 1);
        row.descriptor = sqlite_as_string(stmt, 2);
        row.band_num = sqlite3_column_int(stmt, 3);
        out.push_back(row);
    }
    sqlite3_finalize(stmt);
    return out;
}

std::vector<image_collection::images_row> image_collection::get_images() {
    std::vector<image_collection::images_row> out;
    std::string sql = "SELECT id, name, left, top, bottom, right, datetime, proj FROM images";
    sqlite3_stmt* stmt;
    sqlite3_prepare_v2(_db, sql.c_str(), -1, &stmt, NULL);
    if (!stmt) {
        throw std::string("ERROR in image_collection::get_images(): cannot prepare query statement");
    }

    while (sqlite3_step(stmt) == SQLITE_ROW) {
        image_collection::images_row row;

        row.id = sqlite3_column_int(stmt, 0);
        row.name = sqlite_as_string(stmt, 1);
        row.left = sqlite3_column_double(stmt, 2);
        row.top = sqlite3_column_double(stmt, 3);
        row.bottom = sqlite3_column_double(stmt, 4);
        row.right = sqlite3_column_double(stmt, 5);
        row.datetime = sqlite_as_string(stmt, 6);
        row.proj = sqlite_as_string(stmt, 7);
        out.push_back(row);
    }
    sqlite3_finalize(stmt);
    return out;
}

std::string image_collection::distinct_srs() {
    std::string out = "";
    std::string sql = "SELECT DISTINCT proj from images;";
    sqlite3_stmt* stmt;
    sqlite3_prepare_v2(_db, sql.c_str(), -1, &stmt, NULL);
    if (!stmt) {
        throw std::string("ERROR in image_collection::distinct_srs(): cannot prepare query statement");
    }

    if (sqlite3_step(stmt) == SQLITE_ROW) {
        out = sqlite_as_string(stmt, 0);
        if (sqlite3_step(stmt) == SQLITE_ROW) {
            // if more than one row in the results, return empty string
            out = "";
        }
    }
    sqlite3_finalize(stmt);
    return out;
}

bool image_collection::is_aligned() {
    bool aligned = false;
    std::string sql = "SELECT DISTINCT \"left\", \"top\", \"bottom\", \"right\", \"proj\" from images;";
    sqlite3_stmt* stmt;
    sqlite3_prepare_v2(_db, sql.c_str(), -1, &stmt, NULL);
    if (!stmt) {
        throw std::string("ERROR in mage_collection::distinct_srs(): cannot prepare query statement");
    }

    if (sqlite3_step(stmt) == SQLITE_ROW) {
        aligned = true;
        if (sqlite3_step(stmt) == SQLITE_ROW) {
            // if more than one row in the results, return false
            aligned = false;
        }
    }
    sqlite3_finalize(stmt);
    return aligned;
}

std::vector<std::string> image_collection::unroll_archives(std::vector<std::string> descriptors) {
    std::vector<std::string> out;

    for (uint32_t i = 0; i < descriptors.size(); ++i) {
        std::string s = descriptors[i];
        if (s.compare(s.length() - 4, 4, ".zip") == 0 || s.compare(s.length() - 4, 4, ".ZIP") == 0) {
            char** y = VSIReadDirRecursive(("/vsizip/" + s).c_str());
            char** x = y;
            if (x != NULL) {
                while (*x != NULL) {
                    out.push_back("/vsizip/" + filesystem::join(s, *x));
                    ++x;
                }
                CSLDestroy(y);
            }

        } else if (s.compare(s.length() - 3, 3, ".gz") == 0 || s.compare(s.length() - 3, 3, ".GZ") == 0) {
            out.push_back("/vsigzip/" + s);
        } else if (s.compare(s.length() - 4, 4, ".tar") == 0 || s.compare(s.length() - 4, 4, ".TAR") == 0 ||
                   s.compare(s.length() - 7, 7, ".tar.gz") == 0 || s.compare(s.length() - 7, 7, ".TAR.GZ") == 0 ||
                   s.compare(s.length() - 4, 4, ".tgz") == 0 || s.compare(s.length() - 4, 4, ".TGZ") == 0) {
            char** y = VSIReadDirRecursive(("/vsitar/" + s).c_str());
            char** x = y;
            if (x != NULL) {
                while (*x != NULL) {
                    out.push_back("/vsitar/" + filesystem::join(s, *x));
                    ++x;
                }
                CSLDestroy(y);
            }
        } else {
            out.push_back(s);
        }
    }
    return out;
}

std::shared_ptr<image_collection> image_collection::create_from_tables(std::vector<std::string> band_name,
                                                                       std::vector<std::string> image_name, std::vector<std::string> image_proj,
                                                                       std::vector<std::string> image_datetime, std::vector<double> image_left,
                                                                       std::vector<double> image_top, std::vector<double> image_bottom,
                                                                       std::vector<double> image_right, std::vector<std::string> gdalrefs_descriptor,
                                                                       std::vector<uint16_t> gdalrefs_band_num) {
    // Make sure that length of all vectors is identical
    if (band_name.size() != image_name.size() ||
        band_name.size() != image_proj.size() ||
        band_name.size() != image_datetime.size() ||
        band_name.size() != image_left.size() ||
        band_name.size() != image_top.size() ||
        band_name.size() != image_bottom.size() ||
        band_name.size() != image_right.size() ||
        band_name.size() != gdalrefs_descriptor.size() ||
        band_name.size() != gdalrefs_band_num.size()) {
        GCBS_ERROR("Arguments must have identical size.");
        throw std::string("Arguments must have identical size.");
    }
    if (band_name.empty()) {
        GCBS_ERROR("Arguments are empty");
        throw std::string("Arguments are empty");
    }

    // create empty image collection
    std::shared_ptr<image_collection> o = std::make_shared<image_collection>();

    std::map<std::string, uint32_t> image_ids;
    std::map<std::string, uint32_t> band_ids;
    uint32_t cur_image_id;
    uint32_t cur_band_id;

    sqlite3_exec(o->get_db_handle(), "BEGIN TRANSACTION;", NULL, NULL, NULL);  // what if this fails?!
    for (uint32_t i = 0; i < image_name.size(); ++i) {
        sqlite3_exec(o->get_db_handle(), "SAVEPOINT s1;", NULL, NULL, NULL);  // what if this fails?!
        auto itband = band_ids.find(band_name[i]);
        if (itband == band_ids.end()) {
            std::string sql_band_insert = "INSERT INTO bands(name) VALUES ('" + sqlite_escape_singlequotes(band_name[i]) + "')";
            if (sqlite3_exec(o->get_db_handle(), sql_band_insert.c_str(), NULL, NULL, NULL) != SQLITE_OK) {
                GCBS_WARN("Failed to add band '" + band_name[i] + "' for dataset at row " + std::to_string(i) + "; dataset will be skipped");
                sqlite3_exec(o->get_db_handle(), "ROLLBACK TO s1;", NULL, NULL, NULL);  // what if this fails?!
                continue;
            }
            cur_band_id = sqlite3_last_insert_rowid(o->get_db_handle());
            band_ids.insert(std::make_pair(band_name[i], cur_band_id));
        } else {
            cur_band_id = itband->second;
        }

        auto itimage = image_ids.find(image_name[i]);
        if (itimage == image_ids.end()) {
            datetime d = datetime::from_string(image_datetime[i]);
            std::string sql_insert_image = "INSERT INTO images(name, datetime, left, top, bottom, right, proj) VALUES('" + sqlite_escape_singlequotes(image_name[i]) + "','" +
                                           d.to_string() + "'," +
                                           std::to_string(image_left[i]) + "," + std::to_string(image_top[i]) + "," + std::to_string(image_bottom[i]) + "," + std::to_string(image_right[i]) + ",'" + sqlite_escape_singlequotes(image_proj[i]) + "')";
            if (sqlite3_exec(o->get_db_handle(), sql_insert_image.c_str(), NULL, NULL, NULL) != SQLITE_OK) {
                GCBS_WARN("Failed to add image '" + image_name[i] + "' for dataset at row " + std::to_string(i) + "; dataset will be skipped");
                sqlite3_exec(o->get_db_handle(), "ROLLBACK TO s1;", NULL, NULL, NULL);  // what if this fails?!
                continue;
            }
            cur_image_id = sqlite3_last_insert_rowid(o->get_db_handle());
            image_ids.insert(std::make_pair(image_name[i], cur_image_id));  // take care of race conditions if things run parallel at some po
        } else {
            cur_image_id = itimage->second;
        }

        std::string sql_insert_gdalref = "INSERT INTO gdalrefs(descriptor, image_id, band_id, band_num) VALUES('" + sqlite_escape_singlequotes(gdalrefs_descriptor[i]) + "'," + std::to_string(cur_image_id) + "," + std::to_string(cur_band_id) + "," + std::to_string(gdalrefs_band_num[i]) + ");";
        if (sqlite3_exec(o->get_db_handle(), sql_insert_gdalref.c_str(), NULL, NULL, NULL) != SQLITE_OK) {
            GCBS_WARN("Failed to add dataset '" + gdalrefs_descriptor[i] + "'; dataset will be skipped");
            sqlite3_exec(o->get_db_handle(), "ROLLBACK TO s1;", NULL, NULL, NULL);  // what if this fails?!
            continue;
        }
        sqlite3_exec(o->get_db_handle(), "RELEASE s1;", NULL, NULL, NULL);  // what if this fails?!
    }
    sqlite3_exec(o->get_db_handle(), "COMMIT TRANSACTION;", NULL, NULL, NULL);  // what if this fails?!
    return o;
}

uint32_t image_collection::insert_band(uint32_t id, std::string name, std::string type, double offset, double scale, std::string unit, std::string nodata) {
    std::string sql_band_insert = "INSERT INTO bands(id, name, type, offset, scale, unit, nodata) VALUES (" + std::to_string(id) +
                                  ",'" + sqlite_escape_singlequotes(name) + "', '" + type + "'," + std::to_string(offset) + "," +
                                  std::to_string(scale) + ",'" + sqlite_escape_singlequotes(unit) + "','" + nodata + "')";
    if (sqlite3_exec(_db, sql_band_insert.c_str(), NULL, NULL, NULL) != SQLITE_OK) {
        GCBS_ERROR("Failed to insert band into image collection database");
        throw std::string("Failed to insert band into image collection database");
    }
    return sqlite3_last_insert_rowid(_db);
}

uint32_t image_collection::insert_band(std::string name, std::string type, double offset, double scale, std::string unit, std::string nodata) {
    std::string sql_band_insert = "INSERT INTO bands(name, type, offset, scale, unit, nodata) VALUES ('" +
                                  sqlite_escape_singlequotes(name) + "', '" + type + "'," + std::to_string(offset) + "," +
                                  std::to_string(scale) + ",'" + sqlite_escape_singlequotes(unit) + "','" + nodata + "')";
    if (sqlite3_exec(_db, sql_band_insert.c_str(), NULL, NULL, NULL) != SQLITE_OK) {
        GCBS_ERROR("Failed to insert band into image collection database");
        throw std::string("Failed to insert band into image collection database");
    }
    return sqlite3_last_insert_rowid(_db);
}

uint32_t image_collection::insert_image(uint32_t id, std::string name, double left, double top, double bottom, double right, std::string datetime, std::string proj) {
    datetime = datetime::from_string(datetime).to_string();
    std::string sql_insert_image = "INSERT INTO images(id, name, datetime, left, top, bottom, right, proj) VALUES(" + std::to_string(id) + ",'" + sqlite_escape_singlequotes(name) + "','" +
                                   datetime + "'," + std::to_string(left) + "," + std::to_string(top) + "," + std::to_string(bottom) + "," + std::to_string(right) + ",'" + sqlite_escape_singlequotes(proj) + "')";
    if (sqlite3_exec(_db, sql_insert_image.c_str(), NULL, NULL, NULL) != SQLITE_OK) {
        GCBS_ERROR("Failed to insert image into image collection database");
        throw std::string("Failed to insert image into image collection database");
    }
    return (sqlite3_last_insert_rowid(_db));
}

uint32_t image_collection::insert_image(std::string name, double left, double top, double bottom, double right, std::string datetime, std::string proj) {
    datetime = datetime::from_string(datetime).to_string();
    std::string sql_insert_image = "INSERT INTO images(name, datetime, left, top, bottom, right, proj) VALUES('" + sqlite_escape_singlequotes(name) + "','" +
                                   datetime + "'," + std::to_string(left) + "," + std::to_string(top) + "," + std::to_string(bottom) + "," + std::to_string(right) + ",'" + sqlite_escape_singlequotes(proj) + "')";
    if (sqlite3_exec(_db, sql_insert_image.c_str(), NULL, NULL, NULL) != SQLITE_OK) {
        GCBS_ERROR("Failed to insert image into image collection database");
        throw std::string("Failed to insert image into image collection database");
    }
    return (sqlite3_last_insert_rowid(_db));
}

void image_collection::insert_dataset(uint32_t image_id, uint32_t band_id, std::string descriptor, uint32_t band_num) {
    std::string sql_insert_gdalref = "INSERT INTO gdalrefs(descriptor, image_id, band_id, band_num) VALUES('" + sqlite_escape_singlequotes(descriptor) + "'," + std::to_string(image_id) + "," + std::to_string(band_id) + "," + std::to_string(band_num) + ");";
    if (sqlite3_exec(_db, sql_insert_gdalref.c_str(), NULL, NULL, NULL) != SQLITE_OK) {
        GCBS_ERROR("Failed to insert dataset into image collection database");
        throw std::string("Failed to insert dataset into image collection database");
    }
}

void image_collection::insert_collection_md(std::string key, std::string value) {
    std::string sql_insert = "INSERT INTO collection_md(key, value) VALUES('" + sqlite_escape_singlequotes(key) + "','" + sqlite_escape_singlequotes(value) + "');";
    if (sqlite3_exec(_db, sql_insert.c_str(), NULL, NULL, NULL) != SQLITE_OK) {
        GCBS_ERROR("Failed to insert collection metadata into image collection database");
        throw std::string("Failed to insert collection metadata into image collection database");
    }
}

void image_collection::insert_band_md(uint32_t band_id, std::string key, std::string value) {
    std::string sql_insert = "INSERT INTO band_md(band_id, key, value) VALUES(" + std::to_string(band_id) + ",'" + sqlite_escape_singlequotes(key) + "','" + sqlite_escape_singlequotes(value) + "');";
    if (sqlite3_exec(_db, sql_insert.c_str(), NULL, NULL, NULL) != SQLITE_OK) {
        GCBS_ERROR("Failed to insert band metadata into image collection database");
        throw std::string("Failed to insert band metadata into image collection database");
    }
}

void image_collection::insert_image_md(uint32_t image_id, std::string key, std::string value) {
    std::string sql_insert = "INSERT INTO image_md(image_id, key, value) VALUES(" + std::to_string(image_id) + ",'" + sqlite_escape_singlequotes(key) + "','" + sqlite_escape_singlequotes(value) + "');";
    if (sqlite3_exec(_db, sql_insert.c_str(), NULL, NULL, NULL) != SQLITE_OK) {
        GCBS_ERROR("Failed to insert image metadata into image collection database");
        throw std::string("Failed to insert image metadata into image collection database");
    }
}

void image_collection::transaction_start() {
    sqlite3_exec(_db, "BEGIN TRANSACTION;", NULL, NULL, NULL);  // what if this fails?!
}

void image_collection::transaction_end() {
    sqlite3_exec(_db, "COMMIT TRANSACTION;", NULL, NULL, NULL);  // what if this fails?!
}

bool image_collection::is_empty() {
    if (count_bands() == 0 || count_images() == 0 || count_gdalrefs() == 0) return true;
    return false;
}

std::string image_collection::sqlite_as_string(sqlite3_stmt* stmt, uint16_t col) {
    const unsigned char* a = sqlite3_column_text(stmt, col);
    if (!a) {
        return std::string("");
    } else {
        return std::string(reinterpret_cast<const char*>(a));
    }
}

std::string image_collection::sqlite_escape_singlequotes(std::string s) {
    if (s.empty()) return s;
    std::size_t pos = 0;
    while((pos = s.find("'", pos)) != std::string::npos) {
        s.replace(pos, 1, "''");
        pos += 2;
    }
    return s;
}

}  // namespace gdalcubes
