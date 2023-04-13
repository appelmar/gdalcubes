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

#ifndef IMAGE_COLLECTION_H
#define IMAGE_COLLECTION_H

#include <ogr_spatialref.h>

#include "collection_format.h"
#include "coord_types.h"
#include "datetime.h"

// SQLite forward declarations
class sqlite3;
class sqlite3_stmt;

namespace gdalcubes {

/**
 * @note copy construction and assignment are deleted because the sqlite must not be shared (handle will be closed in destructor). Instrad, use
 * std::shared_ptr<image_collection> to share the whole image collection resource if needed.
  */
class image_collection {
   public:
    struct bands_row {
        uint32_t id;
        std::string name;
        GDALDataType type;
        double offset;
        double scale;
        std::string unit;
        std::string nodata;
        uint32_t image_count;
    };

    struct images_row {
        uint32_t id;
        std::string name;
        double left;
        double top;
        double bottom;
        double right;
        std::string datetime;
        std::string proj;
    };

    struct gdalrefs_row {
        uint32_t image_id;
        uint32_t band_id;
        std::string descriptor;
        uint16_t band_num;
    };

   public:
    /**
     * Default constructor, creates an empty image collection
     * in a temporary SQLite database
     */
    image_collection();

    /**
     * Constructs an empty image collection with given format
     * @param format
     */
    image_collection(collection_format format);

    /**
     * Opens an existing image collection from a file
     * @param filename
     */
    image_collection(std::string filename);

    ~image_collection();
    image_collection(const image_collection&) = delete;
    void operator=(const image_collection&) = delete;

    // move constructor
    image_collection(image_collection&& A) : _format(A._format), _filename(A._filename), _db(A._db) {}

    static std::shared_ptr<image_collection> create(collection_format format, std::vector<std::string> descriptors, bool strict = true);
    static std::shared_ptr<image_collection> create(std::vector<std::string> descriptors, std::vector<std::string> date_time,
                                                    std::vector<std::string> band_names = {}, bool use_subdatasets = false, bool one_band_per_file = false);
    static std::shared_ptr<image_collection> create();

    std::string to_string();

    void add_with_collection_format(std::vector<std::string> descriptors, bool strict = true);
    void add_with_collection_format(std::string descriptor, bool strict = true);

    void add_with_datetime(std::vector<std::string> descriptors, std::vector<std::string> date_time, std::vector<std::string> band_names = {}, bool use_subdatasets = false);
    void add_with_datetime_bands(std::vector<std::string> descriptors, std::vector<std::string> date_time,
                                                   std::vector<std::string> band_names, bool use_subdatasets = false);

    void write(const std::string filename);

    /**
     * Stores an image collection as a new temporary database, same as image_collection::write("")
     */
    inline void temp_copy() { return write(""); }

    /**
     * Check if the database schema is complete, i.e., the database is ready for e.g. adding images.
     * @return true, if the database is OK, false otherwise (e.g. if a table is missing or there are no bands)
     * @note NOT YET IMPLEMENTED
     */
    bool is_valid();

    /**
    * Check if the database has no images and/or bands
    * @return true, if the database is empty, i.e. does not point to any bands or datasets
    */
    bool is_empty();

    /**
     * Check if the collection has gdalrefs for all image X band combinations, i.e. whether there is no image with
     * missing bands.
     * @return true, if the collection is complete, false otherwise.
     * @note NOT YET IMPLEMENTED
     */
    bool is_complete();

    inline bool is_temporary() { return _filename.empty(); }

    /**
     * Removes images and corresponding gdal dataset references captured before start or after end from the database.
     * @param start Posix time representing the start datetime of the range
     * @param end Posix time representing the end datetime of the range
     * @note This operations works in-place, it will overwrite the opened database file. Consider calling `temp_copy()` to create
     * a temporary copy of the database before.
     * @note NOT YET IMPLEMENTED
     */
    void filter_datetime_range(date::sys_seconds start, date::sys_seconds end);

    /**
     * @see image_collection::filter_datetime_range(std::tm start, std::tm end)
     * @param start start datetime as ISO8601 string
     * @param end end datetime as ISO8601 string
     */
    void filter_datetime_range(std::string start, std::string end) {
        date::sys_seconds tstart = datetime::tryparse("%Y-%m-%dT%H:%M:%S", start);
        date::sys_seconds tend = datetime::tryparse("%Y-%m-%dT%H:%M:%S", end);
        filter_datetime_range(tstart, tend);
    }

    /**
     * Removes unneeded bands and the corresponding gdal dataset references.
     * @param bands std::vector of band names, all bands not in the list will be removed from the collection.
     * @note This operations works in-place, it will overwrite the opened database file. Consider calling `temp_copy()` to create
     * a temporary copy of the database before.
     * @note NOT YET IMPLEMENTED
     */
    void filter_bands(std::vector<std::string> bands);

    /**
    * Removes images that are located outside of a given rectangular range.
    * @param range the coordinates of the rectangular spatial range
    * @param proj Projection string is a format readable by GDAL (typically proj4 or proj5).
    * @note This operations works in-place, it will overwrite the opened database file. Consider calling `temp_copy()` to create
    * a temporary copy of the database before.
    * @note NOT YET IMPLEMENTED
    */
    void filter_spatial_range(bounds_2d<double> range, std::string proj);

    uint16_t count_bands();
    uint32_t count_images();
    uint32_t count_gdalrefs();

    struct find_range_st_row {
        find_range_st_row() : image_id(0), image_name(""), descriptor(""), datetime(""), band_name(""), band_num(1) {}
        uint32_t image_id;
        std::string image_name;
        std::string descriptor;
        std::string datetime;
        std::string band_name;
        uint16_t band_num;
        std::string srs;
    };
    std::vector<find_range_st_row> find_range_st(bounds_st range, std::string srs,
                                                 std::vector<std::string> bands, std::vector<std::string> order_by = {});
    inline std::vector<find_range_st_row> find_range_st(bounds_st range, std::string srs, std::vector<std::string> order_by = {}) {
        return find_range_st(range, srs, std::vector<std::string>(), order_by);
    };

    /**
     * Return available bands of an image collection. Bands without
     * correspoding datasets are omitted.
     * @return
     */
    std::vector<image_collection::bands_row> get_available_bands();

    /**
     * Return all bands, including bands of the collection format
     * without corresponding datasets in the collection
     * @return
     */
    std::vector<image_collection::bands_row> get_all_bands();

    std::vector<image_collection::gdalrefs_row> get_gdalrefs();

    std::vector<image_collection::images_row> get_images();

    /**
     * Helper function to create image collections from full tables
     *
     * @note Rows with identical image names are expected to have identical datetime, left, right, bottom, top, and proj values
     * @note Currently, it is not possible to specify additional band metadata (offset, scale, unit, etc.)
     * @note All vector arguments represent columns of a table and must have identical sizes.
     */
    static std::shared_ptr<image_collection> create_from_tables(std::vector<std::string> band_name,
                                                                std::vector<std::string> image_name,
                                                                std::vector<std::string> image_proj,
                                                                std::vector<std::string> image_datetime,
                                                                std::vector<double> image_left,
                                                                std::vector<double> image_top,
                                                                std::vector<double> image_bottom,
                                                                std::vector<double> image_right,
                                                                std::vector<std::string> gdalrefs_descriptor,
                                                                std::vector<uint16_t> gdalrefs_band_num);

    // Lowe level functions to create image collections programmatically

    uint32_t insert_band(uint32_t id, std::string name, std::string type = "", double offset = 0.0, double scale = 1.0, std::string unit = "", std::string nodata = "");
    uint32_t insert_band(std::string name, std::string type = "", double offset = 0.0, double scale = 1.0, std::string unit = "", std::string nodata = "");

    uint32_t insert_image(uint32_t id, std::string name, double left, double top, double bottom, double right, std::string datetime, std::string proj);
    uint32_t insert_image(std::string name, double left, double top, double bottom, double right, std::string datetime, std::string proj);

    void insert_dataset(uint32_t image_id, uint32_t band_id, std::string descriptor, uint32_t band_num = 1);

    void insert_collection_md(std::string key, std::string value);
    void insert_band_md(uint32_t band_id, std::string key, std::string value);
    void insert_image_md(uint32_t image_id, std::string key, std::string value);

    void transaction_start();
    void transaction_end();

    /**
     * Derive the size of a pixel for one or all bands in bytes
     * @param band band identifier, if emtpy the sum of all bands is used
     * @return the size of a pixel with one or all bands or zero if the specified band does not exist.
     */
    uint16_t pixel_size_bytes(std::string band = "");

    /**
     * Derive the spatial and temporal extent of all images
     * @return the spatial and temporal extent
     */
    bounds_st extent();

    inline std::string get_filename() { return _filename; }

    /**
     * @brief Check whether all images in a collection have the same SRS and spatial extent
     * @return true, if the image collection is aligned
     */
    bool is_aligned();

    /**
     * Check whether all images in a collection have the same SRS and if yes, return the SRS
     * @return SRS (usually proj4 string) if all images have the same SRS, empty string otherwise
     */
    std::string distinct_srs();

    /**
     * @brief Replace .zip .gz .tar files from given list by their content as virtual file system GDAL descriptors
     * @see https://www.gdal.org/gdal_virtual_file_systems.html
     * @param descriptors input list of filenames
     * @note This function is not recursive, i.e., it will not unroll .zip files within .zip files etc.
     * @return list of filenames with unrolled archive and or compressed files using GDAL VSI
     */
    static std::vector<std::string> unroll_archives(std::vector<std::string> descriptors);

    /**
     * Return a pointer to the SQlite handle.
     * Ownership is not transferred, i.e. do NOT call sqlite3_close() on the returned object.
     * @return
     */
    sqlite3* get_db_handle();




   protected:
    collection_format _format;
    std::string _filename;
    sqlite3* _db;

    static std::string sqlite_as_string(sqlite3_stmt* stmt, uint16_t col);


    static std::string sqlite_escape_singlequotes(std::string s);

    /**
     * Add a single image to the collection, where one GDAL dataset has spatial dimensions and variables / spectral bands
     * @param descriptor GDAL dataset descriptor
     */
    void add_spaceband_image(const std::string& descriptor);

    /**
     * Add a single image to the collection, where one GDAL dataset has spatial dimensions and time encoded as bands but only one variable
     * @param descriptor GDAL dataset descriptor
     */
    void add_spacetime_image(const std::string& descriptor);
};

}  // namespace gdalcubes

#endif  //IMAGE_COLLECTION_H
