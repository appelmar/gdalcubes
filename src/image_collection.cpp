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

#include "image_collection.h"

#include <regex>
#include "config.h"
#include "external/date.h"
#include "filesystem.h"
#include "utils.h"

image_collection::image_collection(collection_format format) : _format(format), _filename(""), _db(nullptr) {
    if (sqlite3_open_v2("", &_db, SQLITE_OPEN_READWRITE | SQLITE_OPEN_CREATE | SQLITE_OPEN_FULLMUTEX, NULL) != SQLITE_OK) {
        std::string msg = "ERROR in image_collection::create(): cannot create temporary image collection file.";
        throw msg;
    }

    // Enable foreign key constraints
    sqlite3_db_config(_db, SQLITE_DBCONFIG_ENABLE_FKEY, 1, NULL);

    // Create tables

    // key value metadata table
    std::string sql_schema_md = "CREATE TABLE md(key TEXT PRIMARY KEY, value TEXT);";
    if (sqlite3_exec(_db, sql_schema_md.c_str(), NULL, NULL, NULL) != SQLITE_OK) {
        throw std::string("ERROR in image_collection::create(): cannot create image collection schema (i).");
    }
    std::string sql_insert_format = "INSERT INTO md(key, value) VALUES('collection_format','" + _format.json().dump() + "');";
    if (sqlite3_exec(_db, sql_insert_format.c_str(), NULL, NULL, NULL) != SQLITE_OK) {
        throw std::string("ERROR in image_collection::create(): cannot insert collection format to database.");
    }

    // Create Bands
    std::string sql_schema_bands = "CREATE TABLE bands (id INTEGER PRIMARY KEY, name TEXT, type VARCHAR(16), offset NUMERIC DEFAULT 0.0, scale NUMERIC DEFAULT 1.0, unit VARCHAR(16) DEFAULT '', nodata VARCHAR(16) DEFAULT '');";
    if (sqlite3_exec(_db, sql_schema_bands.c_str(), NULL, NULL, NULL) != SQLITE_OK) {
        throw std::string("ERROR in collection_format::apply(): cannot create image collection schema (ii).");
    }
    if (!_format.json().count("bands")) {
        throw std::string("ERROR in collection_format::apply(): image collection format does not contain any bands.");
    }

    if (_format.json()["bands"].size() == 0) {
        throw std::string("ERROR in collection_format::apply(): image collection format does not contain any bands.");
    }

    uint16_t band_id = 0;
    for (auto it = _format.json()["bands"].begin(); it != _format.json()["bands"].end(); ++it) {
        std::string sql_insert_band;
        if (it.value().count("nodata")) {
            sql_insert_band = "INSERT INTO bands(id, name, nodata) VALUES(" + std::to_string(band_id) + ",'" + it.key() + "','" + std::to_string(it.value()["nodata"].get<double>()) + "');";
        } else {
            sql_insert_band = "INSERT INTO bands(id, name) VALUES(" + std::to_string(band_id) + ",'" + it.key() + "');";
        }
        ++band_id;
        if (sqlite3_exec(_db, sql_insert_band.c_str(), NULL, NULL, NULL) != SQLITE_OK) {
            throw std::string("ERROR in collection_format::apply(): cannot insert bands to collection database.");
        }
    }

    // Create image table
    std::string sql_schema_images = "CREATE TABLE images (id INTEGER PRIMARY KEY, name TEXT, left NUMERIC, top NUMERIC, bottom NUMERIC, right NUMERIC, datetime TEXT, proj TEXT, UNIQUE(name));CREATE INDEX idx_image_names ON images(name);";
    if (sqlite3_exec(_db, sql_schema_images.c_str(), NULL, NULL, NULL) != SQLITE_OK) {
        throw std::string("ERROR in collection_format::apply(): cannot create image collection schema (iii).");
    }

    std::string sql_schema_gdalrefs =
        "CREATE TABLE gdalrefs (image_id INTEGER, band_id INTEGER, descriptor TEXT, band_num INTEGER, FOREIGN KEY (image_id) REFERENCES images(id) ON DELETE CASCADE, PRIMARY KEY (image_id, band_id), FOREIGN KEY (band_id) REFERENCES bands(id) ON DELETE CASCADE);"
        "CREATE INDEX idx_gdalrefs_bandid ON gdalrefs(band_id);"
        "CREATE INDEX idx_gdalrefs_imageid ON gdalrefs(image_id);";
    if (sqlite3_exec(_db, sql_schema_gdalrefs.c_str(), NULL, NULL, NULL) != SQLITE_OK) {
        throw std::string("ERROR in collection_format::apply(): cannot create image collection schema (iv).");
    }
}

image_collection::image_collection(std::string filename) : _format(), _filename(filename), _db(nullptr) {
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
    std::string sql_select_format = "SELECT value FROM md WHERE key='collection_format';";
    sqlite3_stmt* stmt;
    sqlite3_prepare_v2(_db, sql_select_format.c_str(), -1, &stmt, NULL);
    if (!stmt) {
        std::string msg = "ERROR in image_collection::image_collection(): cannot extract collection format from existing image collection file.";
        throw msg;
    }
    if (sqlite3_step(stmt) == SQLITE_DONE) {
        std::string msg = "ERROR in image_collection::image_collection(): cannot extract collection format from existing image collection file.";
        throw msg;
    } else {
        _format.load_string(sqlite_as_string(stmt, 0));
    }
    sqlite3_finalize(stmt);
}

std::shared_ptr<image_collection> image_collection::create(collection_format format, std::vector<std::string> descriptors, bool strict) {
    std::shared_ptr<image_collection> o = std::make_shared<image_collection>(format);
    o->add(descriptors, strict);
    return o;
}

struct image_band {
    GDALDataType type;
    std::string unit;
    double scale;
    double offset;
    std::string nodata;
};

void image_collection::add(std::vector<std::string> descriptors, bool strict) {
    std::vector<std::regex> regex_band_pattern;

    /* TODO: The following will fail if other applications create image collections and assign ids to bands differently. A better solution would be to load band ids, names, and nums from the database bands table directly
     */

    /**
     * TODO: enable .zip, .gz, .tar, .tar.gz...
     */
    std::vector<std::string> band_name;
    std::vector<uint16_t> band_num;
    std::vector<uint16_t> band_ids;
    std::vector<bool> band_complete;

    uint16_t band_id = 0;
    for (nlohmann::json::iterator it = _format.json()["bands"].begin(); it != _format.json()["bands"].end(); ++it) {
        band_name.push_back(it.key());
        regex_band_pattern.push_back(std::regex(it.value()["pattern"].get<std::string>()));
        if (it.value().count("band")) {
            band_num.push_back(it.value()["band"].get<int>());
        } else {
            band_num.push_back(1);
        }
        band_ids.push_back(band_id);
        band_complete.push_back(false);
        ++band_id;
    }

    std::string global_pattern = "";
    if (_format.json().count("pattern")) global_pattern = _format.json()["pattern"].get<std::string>();
    std::regex regex_global_pattern(global_pattern);

    if (!_format.json()["images"].count("pattern"))
        throw std::string("ERROR in image_collection::add(): image collection format does not contain a composition rule for images.");
    std::regex regex_images(_format.json()["images"]["pattern"].get<std::string>());

    // @TODO: Make datetime optional, e.g., for DEMs
    std::string datetime_format = "%Y-%m-%d";
    if (!_format.json().count("datetime") || !_format.json()["datetime"].count("pattern")) {
        throw std::string("ERROR in image_collection::add(): image collection format does not contain a rule to derive date/time.");
    }
    std::regex regex_datetime(_format.json()["datetime"]["pattern"].get<std::string>());
    if (_format.json()["datetime"].count("format")) {
        datetime_format = _format.json()["datetime"]["format"].get<std::string>();
    }

    uint32_t counter = -1;
    std::shared_ptr<progress> p = config::instance()->get_default_progress_bar()->get();
    p->set(0);  // explicitly set to zero to show progress bar immediately
    for (auto it = descriptors.begin(); it != descriptors.end(); ++it) {
        ++counter;

        // if .zip, .tar.gz, .gz, or .tar

        p->set((double)counter / (double)descriptors.size());
        if (!global_pattern.empty()) {  // prevent unnecessary GDALOpen calls
            if (!std::regex_match(*it, regex_global_pattern)) {
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
        char* proj4;
        if (dataset->GetGeoTransform(affine_in) != CE_None) {
            // No affine transformation, maybe GCPs?
            GDALClose((GDALDatasetH)dataset);
            if (strict) throw std::string("ERROR in image_collection::add(): GDAL cannot derive affine transformation for '" + *it + "'. GCPs or unreferenced images are currently not supported.");
            GCBS_WARN("Failed to derive affine transformation from " + *it);
            continue;
        } else {
            bbox.left = affine_in[0];
            bbox.right = affine_in[0] + affine_in[1] * dataset->GetRasterXSize() + affine_in[2] * dataset->GetRasterYSize();
            bbox.top = affine_in[3];
            bbox.bottom = affine_in[3] + affine_in[4] * dataset->GetRasterXSize() + affine_in[5] * dataset->GetRasterYSize();
            OGRSpatialReference srs_in;
            srs_in.SetFromUserInput(dataset->GetProjectionRef());
            srs_in.exportToProj4(&proj4);
            bbox.transform(proj4, "EPSG:4326");
        }

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
        GDALClose((GDALDatasetH)dataset);

        // TODO: check consistency for all files of an image?!
        // -> add parameter checks=true / false

        std::cmatch res_image;
        if (!std::regex_match(it->c_str(), res_image, regex_images)) {
            if (strict) throw std::string("ERROR in image_collection::add(): image composition rule failed for " + std::string(*it));
            GCBS_WARN("Skipping " + *it + " due to failed image composition rule");
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
            std::cmatch res_datetime;
            if (!std::regex_match(it->c_str(), res_datetime, regex_datetime)) {  // not sure to continue or throw an exception here...
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
                                           std::to_string(bbox.left) + "," + std::to_string(bbox.top) + "," + std::to_string(bbox.bottom) + "," + std::to_string(bbox.right) + ",'" + proj4 + "')";
            if (sqlite3_exec(_db, sql_insert_image.c_str(), NULL, NULL, NULL) != SQLITE_OK) {
                if (strict) throw std::string("ERROR in image_collection::add(): cannot add image to images table.");
                GCBS_WARN("Skipping " + *it + " due to failed image table insert");
                continue;
            }
            image_id = sqlite3_last_insert_rowid(_db);  // take care of race conditions if things run parallel at some point

        } else {
            image_id = sqlite3_column_int(stmt, 0);
            // TODO: if checks, compare l,r,b,t, datetime,proj4 from images table with current GDAL dataset
        }
        sqlite3_finalize(stmt);

        // Insert into gdalrefs table

        for (uint16_t i = 0; i < band_name.size(); ++i) {
            if (std::regex_match(*it, regex_band_pattern[i])) {
                // TODO: if checks, check whether bandnum exists in GDALdataset
                // TODO: if checks, compare band type, offset, scale, unit, etc. with current GDAL dataset

                if (!band_complete[i]) {
                    std::string sql_band_update = "UPDATE bands SET type='" + utils::string_from_gdal_type(bands[band_num[i] - 1].type) + "'," +
                                                  "scale=" + std::to_string(bands[band_num[i] - 1].scale) + "," +
                                                  "offset=" + std::to_string(bands[band_num[i] - 1].offset) + "," +
                                                  "unit='" + bands[band_num[i] - 1].unit + "' WHERE name='" + band_name[i] + "';";
                    // if has no data

                    if (sqlite3_exec(_db, sql_band_update.c_str(), NULL, NULL, NULL) != SQLITE_OK) {
                        if (strict) throw std::string("ERROR in image_collection::add(): cannot update band table.");
                        GCBS_WARN("Skipping " + *it + " due to failed band table update");
                        continue;
                    }
                    band_complete[i] = true;
                }

                std::string sql_insert_gdalref = "INSERT INTO gdalrefs(descriptor, image_id, band_id, band_num) VALUES('" + *it + "'," + std::to_string(image_id) + "," + std::to_string(band_ids[i]) + "," + std::to_string(band_num[i]) + ");";
                if (sqlite3_exec(_db, sql_insert_gdalref.c_str(), NULL, NULL, NULL) != SQLITE_OK) {
                    if (strict) throw std::string("ERROR in image_collection::add(): cannot add dataset to gdalrefs table.");
                    GCBS_WARN("Skipping " + *it + "  due to failed gdalrefs insert");
                    break;  // break only works because there is nothing after the loop.
                }
            }
        }
        CPLFree(proj4);
    }
    p->set(1);
    p->finalize();
}

void image_collection::add(std::string descriptor, bool strict) {
    std::vector<std::string> x{descriptor};
    return add(x, strict);
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
                                                                                 std::vector<std::string> bands, std::string order_by) {
    bounds_2d<double> range_trans = (srs == "EPSG:4326") ? range.s : range.s.transform(srs, "EPSG:4326");
    std::string sql =
        "SELECT images.name, gdalrefs.descriptor, images.datetime, bands.name, gdalrefs.band_num "
        "FROM images INNER JOIN gdalrefs ON images.id = gdalrefs.image_id INNER JOIN bands ON gdalrefs.band_id = bands.id WHERE "
        "images.datetime >= '" +
        range.t0.to_string(datetime_unit::SECOND) + "' AND images.datetime <= '" + range.t1.to_string(datetime_unit::SECOND) +
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
    if (!order_by.empty()) {  // explicitly test order by column if given to avoid SQL injection
        if (order_by == "images.name" ||
            order_by == "gdalrefs.descriptor" ||
            order_by == "images.datetime" ||
            order_by == "bands.name" ||
            order_by == "gdalrefs.band_num")
            sql += " ORDER BY " + order_by;

        else {
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
        r.image_name = sqlite_as_string(stmt, 0);
        r.descriptor = sqlite_as_string(stmt, 1);
        r.datetime = sqlite_as_string(stmt, 2);
        r.band_name = sqlite_as_string(stmt, 3);
        r.band_num = sqlite3_column_int(stmt, 4);

        out.push_back(r);
    }
    sqlite3_finalize(stmt);
    return out;
}

std::vector<image_collection::bands_row> image_collection::get_bands() {
    std::vector<image_collection::bands_row> out;

    std::string sql = "SELECT id, name, type, offset,scale, unit, nodata FROM bands ORDER BY name;";  // changing the order my have consequences to data read implementations
    sqlite3_stmt* stmt;
    sqlite3_prepare_v2(_db, sql.c_str(), -1, &stmt, NULL);
    if (!stmt) {
        throw std::string("ERROR in image_collection::get_bands(): cannot prepare query statement");
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
        out.push_back(row);
    }
    sqlite3_finalize(stmt);
    return out;
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

std::string image_collection::sqlite_as_string(sqlite3_stmt* stmt, uint16_t col) {
    const unsigned char* a = sqlite3_column_text(stmt, col);
    if (!a) {
        return std::string("");
    } else {
        return std::string(reinterpret_cast<const char*>(a));
    }
}