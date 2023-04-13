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

#ifndef CUBE_H
#define CUBE_H

#include <mutex>
#include <set>

#include "config.h"
#include "view.h"

namespace gdalcubes {

typedef std::array<uint32_t, 4> cube_coordinate_btyx;
typedef std::array<uint32_t, 3> cube_coordinate_tyx;
typedef std::array<uint32_t, 4> cube_size_btyx;
typedef std::array<uint32_t, 3> cube_size_tyx;
typedef std::array<uint32_t, 4> chunk_coordinate_btyx;
typedef std::array<uint32_t, 3> chunk_coordinate_tyx;
typedef std::array<uint32_t, 4> chunk_size_btyx;
typedef std::array<uint32_t, 3> chunk_size_tyx;

typedef uint32_t chunkid_t;

/**
 * Data structure for defining packed data exports (dividing and adding data values by a scale and offset value
 * respectively, before conversion to smaller integer types during the export.
 */
struct packed_export {
    // Enum definitions for output types, PACK_NONE will not
    // apply any packing. PACK_FLOAT will convert 8 byte doubles
    // to 4 byte floats but IGNORE provided scale and offset values
    enum class packing_type {
        PACK_NONE,
        PACK_UINT8,
        PACK_UINT16,
        PACK_UINT32,
        PACK_INT16,
        PACK_INT32,
        PACK_FLOAT32
    };
    packing_type type;

    /**
     * A vector of per-band scale factors. Must have size() == 1 or size() == nbands
     * If the vector contains only one element, the same scale factor is assumed
     * for all bands.
     */
    std::vector<double> scale;

    /**
     * A vector of per-band offset values. Must have size() == 1 or size() == nbands
     * If the vector contains only one element, the same offset value is assumed
     * for all bands.
     */
    std::vector<double> offset;

    /**
     * A vector of new nodata values. Must have size() == 1 or size() == nbands
     * If the vector contains only one element, the same nodata value is assumed
     * for all bands.
     */
    std::vector<double> nodata;

    static packed_export make_none() {
        packed_export out;
        out.type = packing_type::PACK_NONE;
        return out;
    }

    static packed_export make_float32() {
        packed_export out;
        out.type = packing_type::PACK_FLOAT32;
        out.scale = {1};
        out.offset = {0};
        out.nodata = {std::numeric_limits<int16_t>::lowest()};
        return out;
    }

    static packed_export make_uint16(double scale, double offset, double nodata) {
        packed_export out;
        out.type = packing_type::PACK_UINT16;
        out.scale = {scale};
        out.offset = {offset};
        out.nodata = {nodata};
        return out;
    }

    static packed_export make_uint16(std::vector<double> scale, std::vector<double> offset, std::vector<double> nodata) {
        packed_export out;
        out.type = packing_type::PACK_UINT16;
        out.scale = scale;
        out.offset = offset;
        out.nodata = nodata;
        return out;
    }

    static packed_export make_int16(double scale, double offset) {
        packed_export out;
        out.type = packing_type::PACK_UINT16;
        out.scale = {scale};
        out.offset = {offset};
        return out;
    }

    static packed_export make_int16(std::vector<double> scale, std::vector<double> offset) {
        packed_export out;
        out.type = packing_type::PACK_INT16;
        out.scale = scale;
        out.offset = offset;
        return out;
    }

    static packed_export make_uint8(double scale, double offset) {
        packed_export out;
        out.type = packing_type::PACK_UINT16;
        out.scale = {scale};
        out.offset = {offset};
        return out;
    }

    static packed_export make_uint8(std::vector<double> scale, std::vector<double> offset) {
        packed_export out;
        out.type = packing_type::PACK_UINT8;
        out.scale = scale;
        out.offset = offset;
        return out;
    }
};

class cube;

class chunk_data;

/**
 * @brief Virtual base class for processing a data cube chunk-wise, i.e. applying the same function over all chunks in a cube.
 */
class chunk_processor {
   public:
    virtual ~chunk_processor() {}

    /**
     * @brief Return the maximum number of threads used to process chunks in parallel.
     * @return Maximum number of parallel threads
     */
    virtual uint32_t max_threads() = 0;

    /**
     * Apply a function f over all chunks of a given data cube c
     * @param c data cube
     * @param f function to be applied over all chunks of c, the function gets the chunk_id, the actual chunk data, and a mutex object as input. The latter is only needed
     * for parallel chunk processing, e.g., for synchronous file writing.
     */
    virtual void
    apply(std::shared_ptr<cube> c, std::function<void(chunkid_t, std::shared_ptr<chunk_data>, std::mutex &)> f) = 0;
};

/**
 * @brief Implementation of the chunk_processor class for single-thread sequential chunk processing
 */
class chunk_processor_singlethread : public chunk_processor {
   public:
    /**
    * @copydoc chunk_processor::max_threads
    */
    uint32_t max_threads() override {
        return 1;
    }

    /**
     * @copydoc chunk_processor::apply
     */
    void apply(std::shared_ptr<cube> c,
               std::function<void(chunkid_t, std::shared_ptr<chunk_data>, std::mutex &)> f) override;
};

/**
 * @brief Implementation of the chunk_processor class for multithreaded parallel chunk processing
 */
class chunk_processor_multithread : public chunk_processor {
   public:
    /**
    * @copydoc chunk_processor::max_threads
    */
    uint32_t max_threads() override {
        return _nthreads;
    }

    /**
     * @brief Construct a multithreaded chunk processor
     * @param nthreads number of threads
     */
    chunk_processor_multithread(uint16_t nthreads) : _nthreads(nthreads) {}

    /**
    * @copydoc chunk_processor::apply
    */
    void apply(std::shared_ptr<cube> c,
               std::function<void(chunkid_t, std::shared_ptr<chunk_data>, std::mutex &)> f) override;

    /**
     * Query the number of threads to be used in parallel chunk processing
     * @return the number of threads
     */
    inline uint16_t get_threads() { return _nthreads; }

   private:
    uint16_t _nthreads;
};

/**
 * @brief A simple structure for band information
 */
struct band {
    band(std::string name) : name(name), no_data_value(std::to_string(NAN)), offset(0), scale(1), unit(""), type("float64") {}

    std::string name;
    std::string no_data_value;
    double offset;
    double scale;
    std::string unit;
    std::string type;
};

/**
 * @brief Collection class for image bands, accessible by name or id
 */
class band_collection {
   public:
    /**
     * Adds a band to a collection
     * @param b band information for the new band
     */
    void add(band b) {
        if (!has(b.name)) {
            _bands.push_back(b);
            _band_idx[b.name] = _bands.size() - 1;
        }
    }

    /**
     * Checks whether the collection contains a band with given name
     * @param name unique band name
     * @return true if the collection contains a band with given name, false otherwise
     */
    inline bool has(std::string name) {
        return _band_idx.find(name) != _band_idx.end();
    }

    /**
     * Query a band from the collection by name
     * @param name name of the band
     * @return band information
     */
    inline band get(std::string name) {
        return _bands[_band_idx[name]];
    }

    /**
     * Query the id of a band from its id
     * @param name name of the band
     * @return if of the corresponding band
     */
    inline uint32_t get_index(std::string name) {
        return _band_idx[name];
    }

    /**
     * Query a band from the collection by id
     * @param id id of the band (zero-based)
     * @return band information
     */
    inline band get(uint32_t i) {
        return _bands[i];
    }

    /**
     * Count the number of bands in the collection
     * @return the number of bands in the collection
     */
    inline uint32_t count() {
        return _bands.size();
    }

   private:
    std::map<std::string, uint32_t> _band_idx;
    std::vector<band> _bands;
};

/**
 * @brief A class for storing actual data of one chunk
 *
 * This class is typically used with smart pointers as
 * std::shared_ptr<chunk_data>
 */
class chunk_data {
   public:
    /**
     * @brief Default constructor that creates an empty chunk
     */
    chunk_data() : _buf(nullptr), _size({{0, 0, 0, 0}}) {}

    ~chunk_data() {
        if (_buf && _size[0] * _size[1] * _size[2] * _size[3] > 0) std::free(_buf);
    }

    /**
     * @brief Count the number of bands in the chunk
     * @return number of bands
     */
    inline uint16_t count_bands() { return _size[0]; }

    /**
     * @brief Count the number of pixels in the chunk
     * @return number of pixels \f$(nx \times ny \times nt)\f$
     */
    inline uint32_t count_values() { return _size[1] * _size[2] * _size[3]; }

    /**
     * @brief Returns the total memory size in bytes that is consumed by the chunk
     * @return size of the chunk in bytes
     */
    uint64_t total_size_bytes() {
        return empty() ? 0 : sizeof(double) * _size[0] * _size[1] * _size[2] * _size[3];
    }

    /**
     * @brief Check whether there is data in the buffer
     * @return true, if there is no data in the buffer (either size == 0, or buf == nullptr)
     */
    bool empty() {
        if (_size[0] * _size[1] * _size[2] * _size[3] == 0) return true;
        if (!_buf) return true;
        return false;
    }

    /**
     * @brief Check if the buffer contains only NAN values
     * @return true, if all values are NAN, or the chunk is empty
     */
    bool all_nan() {
        if (empty()) return true;
        for (uint32_t i=0; i< _size[0] * _size[1] * _size[2] * _size[3]; ++i) {
            if (!std::isnan(((double*)_buf)[i])) {
                return false;
            }
        }
        return true;
    }

    /**
     * @brief Access the raw buffer where the data is stored in memory
     *
     * This method is dangerous and provides direct access to the data buffer. Use with caution and never free any memory /
     * remove / add vector elements if you don't know exactly what you do.
     * @return void pointer pointing to the data buffer
     */
    inline void *buf() { return _buf; }

    /**
     * @brief (Re)set the raw buffer where the data is stored in memory
     *
     * This method is dangerous and provides direct access to the data buffer. Use with caution and never free any memory /
     * remove / add vector elements if you don't know exactly what you do.
     *
     * @param b new buffer object, this class takes the ownership, i.e., eventually std::frees memory automatically in the destructor.
     */
    inline void buf(void *b) {
        if (_buf && _size[0] * _size[1] * _size[2] * _size[3] > 0) std::free(_buf);
        _buf = b;
    }

    /**
     * @brief Query the size of the contained data
     *
     * The result is an array of size 4 representing the size of dimensions in the order (bands, time, y, x).
     * @return std::array<uint32, 4> with number of cells with regard to bands, time, y, and x
     */
    inline coords_nd<uint32_t, 4> size() { return _size; }

    /**
     * @brief Set the size of the chunk data
     *
     * This method is dangerous, use with caution and never change the size without changing the buffer accordingly.
     *
     * @param s new size
     */
    inline void size(coords_nd<uint32_t, 4> s) { _size = s; }

    /**
     * @brief Write a single chunk to a netCDF file
     *
     * @details
     * The produced netCDF file will NOT contain the spatial reference
     * and dimension metadata, i.e. only the size of dimensions and the actual data
     * @param path output path
     * @param netCDF compression level (0 = no compression, 1 = fastest, 9 = smallest)
     * @force force creation of a netCDF file even if the chunk is empty or NaN only
     */
    void write_ncdf(std::string path, uint8_t compression_level = 0, bool force = false);

    /**
     * @brief Read a single chunk from a netCDF file
     * @param path input path
     */
    void read_ncdf(std::string path);

   private:
    void *_buf;
    chunk_size_btyx _size;
};

/**
 * @brief Base class for all data cube types.
 *
 * This class is virtual but provides a few default methods
 * such as chunk-wise NetCDF export, conversion of chunk-local, cube-local, and spatial coordinates / datetime.
 *
 *
 * @note Specific subclass instances should be created by using static create() methods instead uf using the constructors of
 * directly. This makes sure that we always work with smart pointers and that connections between cubes are maintained
 * consistently.
 *
 * @note Cubes should primarily be used as smart pointers, i.e. std::shared_ptr<xxx_cube>.
 */
class cube : public std::enable_shared_from_this<cube> {
   public:
    /**
     * @brief Create an empty data cube
     */
    cube() : _st_ref(nullptr), _chunk_size({16, 256, 256}), _bands() {}

    /**
     * @brief Create an empty data cube with given spacetime reference
     * @param st_ref space time reference (extent, size, SRS) of the cube
     */
    cube(std::shared_ptr<cube_stref> st_ref) : _st_ref(st_ref), _chunk_size(), _bands(), _pre(), _succ() {
        _chunk_size = {16, 256, 256};

        // TODO: add bands
    }

    virtual ~cube() {}

    /**
     * @brief Find the chunk that contains a given point
     * @param p point in spacetime coordinates
     * @return a unique chunk identifier (uint32_t)
     */
    chunkid_t find_chunk_that_contains(coords_st p) const {
        uint32_t cumprod = 1;
        chunkid_t id = 0;

        // 1. Convert map coordinates to view-based coordinates
        coords_nd<uint32_t, 3> s = _st_ref->cube_coords(p);

        // 2. Find out which chunk contains the integer view coordinates
        id += cumprod * (s[2] / _chunk_size[2]);
        cumprod *= (uint32_t)std::ceil(((double)_st_ref->nx()) / ((double)_chunk_size[2]));
        id += cumprod * (s[1] / _chunk_size[1]);
        cumprod *= (uint32_t)std::ceil(((double)(_st_ref->ny()) / ((double)_chunk_size[1])));
        id += cumprod * (s[0] / _chunk_size[0]);
        //cumprod *= (uint32_t)std::ceil(((double)(_c->view().nt()) / ((double)_size_t);

        return id;
    }

    /**
     * @brief Calculate the limits of a given chunk
     * @param c chunk coordinates
     * @return the limits in t,y, and x, as integer data cube coordinates
     */
    bounds_nd<uint32_t, 3> chunk_limits(chunk_coordinate_tyx c) const {
        cube_coordinate_tyx out_vcoords_low;
        cube_coordinate_tyx out_vcoords_high;

        out_vcoords_low[0] = c[0] * _chunk_size[0];
        out_vcoords_low[1] = c[1] * _chunk_size[1];
        out_vcoords_low[2] = c[2] * _chunk_size[2];
        out_vcoords_high[0] = out_vcoords_low[0] + _chunk_size[0] - 1;
        out_vcoords_high[1] = out_vcoords_low[1] + _chunk_size[1] - 1;
        out_vcoords_high[2] = out_vcoords_low[2] + _chunk_size[2] - 1;

        // Shrink to view
        if (out_vcoords_high[0] >= _st_ref->nt())
            out_vcoords_high[0] = _st_ref->nt() - 1;
        if (out_vcoords_low[0] >= _st_ref->nt())
            out_vcoords_low[0] = _st_ref->nt() - 1;

        if (out_vcoords_high[1] >= _st_ref->ny())
            out_vcoords_high[1] = _st_ref->ny() - 1;
        if (out_vcoords_low[1] >= _st_ref->ny())
            out_vcoords_low[1] = _st_ref->ny() - 1;

        if (out_vcoords_high[2] >= _st_ref->nx())
            out_vcoords_high[2] = _st_ref->nx() - 1;
        if (out_vcoords_low[2] >= _st_ref->nx())
            out_vcoords_low[2] = _st_ref->nx() - 1;

        bounds_nd<uint32_t, 3> out;
        out.low = out_vcoords_low;
        out.high = out_vcoords_high;
        return out;
    };

    /**
     * @brief Calculate the limits of a given chunk
     * @param id chunk id
     * @return the limits in t,y, and x, as integer data cube coordinates
     */
    bounds_nd<uint32_t, 3> chunk_limits(chunkid_t id) const {
        coords_nd<uint32_t, 3> out_vcoords_low;
        coords_nd<uint32_t, 3> out_vcoords_high;

        // Transform to global coordinates based on chunk id
        uint32_t cumprod = 1;
        int32_t n;

        n = (uint32_t)std::ceil(((double)(_st_ref->nx())) / ((double)_chunk_size[2]));
        out_vcoords_low[2] = (id / cumprod) % n;
        out_vcoords_low[2] *= _chunk_size[2];
        out_vcoords_high[2] = out_vcoords_low[2] + _chunk_size[2] - 1;
        cumprod *= n;

        n = (uint32_t)std::ceil(((double)(_st_ref->ny())) / ((double)_chunk_size[1]));
        out_vcoords_low[1] = (id / cumprod) % n;
        out_vcoords_low[1] *= _chunk_size[1];
        out_vcoords_high[1] = out_vcoords_low[1] + _chunk_size[1] - 1;
        cumprod *= n;

        n = (uint32_t)std::ceil(((double)(_st_ref->nt())) / ((double)_chunk_size[0]));
        out_vcoords_low[0] = (id / cumprod) % n;
        out_vcoords_low[0] *= _chunk_size[0];
        out_vcoords_high[0] = out_vcoords_low[0] + _chunk_size[0] - 1;
        cumprod *= n;

        // Shrink to view
        if (out_vcoords_high[0] >= _st_ref->nt())
            out_vcoords_high[0] = _st_ref->nt() - 1;
        if (out_vcoords_low[0] >= _st_ref->nt())
            out_vcoords_low[0] = _st_ref->nt() - 1;

        if (out_vcoords_high[1] >= _st_ref->ny())
            out_vcoords_high[1] = _st_ref->ny() - 1;
        if (out_vcoords_low[1] >= _st_ref->ny())
            out_vcoords_low[1] = _st_ref->ny() - 1;

        if (out_vcoords_high[2] >= _st_ref->nx())
            out_vcoords_high[2] = _st_ref->nx() - 1;
        if (out_vcoords_low[2] >= _st_ref->nx())
            out_vcoords_low[2] = _st_ref->nx() - 1;

        bounds_nd<uint32_t, 3> out;
        out.low = out_vcoords_low;
        out.high = out_vcoords_high;
        return out;
    };

    /**
     * @brief Default string description method for data cubes
     * @return std::string with a human-readable description of the data cube
     */
    virtual std::string to_string() {
        // TODO: implement default to_string method here

        return std::string("cube::to_string() has no default method yet.");
    }

    /**
     * @brief Count the total number of chunks of the cube
     * @return Total number of chunks
     */
    inline uint32_t count_chunks() const {
        return count_chunks_x() * count_chunks_y() * count_chunks_t();
    }

    /**
    * @brief Count the number of chunks of the cube in x direction
    * @return Total number of chunk in x direction
    */
    inline uint32_t count_chunks_x() const {
        return std::ceil((double)_st_ref->nx() / (double)_chunk_size[2]);
    }

    /**
    * @brief Count the number of chunks of the cube in y direction
    * @return Total number of chunk in y direction
    */
    inline uint32_t count_chunks_y() const {
        return std::ceil((double)_st_ref->ny() / (double)_chunk_size[1]);
    }

    /**
    * @brief Count the number of chunks of the cube in t direction
    * @return Total number of chunk in t direction
    */
    inline uint32_t count_chunks_t() const {
        return std::ceil((double)_st_ref->nt() / (double)_chunk_size[0]);
    }

    /**
     * @brief Convert a given one-dimensional chunk id to chunk coordinates
     * @param id chunk id
     * @return chunk coordinates in t,y, and x directions
     */
    chunk_coordinate_tyx chunk_coords_from_id(chunkid_t id) {
        chunk_coordinate_tyx out;

        uint32_t cumprod = 1;
        int32_t n;

        n = count_chunks_x();
        out[2] = (id / cumprod) % n;
        cumprod *= n;
        n = count_chunks_y();
        out[1] = (id / cumprod) % n;
        cumprod *= n;
        n = count_chunks_t();
        out[0] = (id / cumprod) % n;
        cumprod *= n;

        return out;
    }

    /**
     * @brief Convert chunk coordinates to a one-dimensional chunk id
     * @param c chunk coordinates in t,y, and x directions
     * @return chunk id
     */
    chunkid_t chunk_id_from_coords(chunk_coordinate_tyx c) {
        return c[0] * count_chunks_y() * count_chunks_x() + c[1] * count_chunks_x() + c[2];
    }

    /**
     * @brief Derive the true size of a specific chunk
     *
     * The size may be different to the _chunk_size setting at the boundary regions.
     * @return
     */
    coords_nd<uint32_t, 3> chunk_size(chunkid_t id) const {
        bounds_nd<uint32_t, 3> vb = chunk_limits(id);
        coords_nd<uint32_t, 3> out;
        out[0] = vb.high[0] - vb.low[0] + 1;
        out[1] = vb.high[1] - vb.low[1] + 1;
        out[2] = vb.high[2] - vb.low[2] + 1;
        return out;
    };

    /**
     * @brief Derive the spatiotemporal extent / bounds of a given chunk
     * @param id Chunk identifier
     * @return Spatiotemporal extent
     */
    bounds_st bounds_from_chunk(chunkid_t id) const {
        bounds_st out_st;

        bounds_nd<uint32_t, 3> vb = chunk_limits(id);
        coords_st low = _st_ref->map_coords(vb.low);
        vb.high[0] += 1;
        vb.high[1] += 1;
        vb.high[2] += 1;
        coords_st high = _st_ref->map_coords(vb.high);

        out_st.s.left = low.s.x;
        out_st.s.right = high.s.x;
        out_st.s.bottom = high.s.y; // Important: high has lower y coordinates (because integer cube coordinates go top -> bottom)
        out_st.s.top = low.s.y; // Important: low has higher y coordinates (because integer cube coordinates go top -> bottom)
        out_st.t0 = low.t;
        out_st.t1 = high.t;

        return out_st;
    }

    /**
     * @brief Get the chunk size of the cube
     *
     * @note This function returns the global chunk size setting of a cube. Specific chunks at the boundary of a cube
     * may be smaller and store less cells.
     * @see cube::chunk_size(chunkid_t id) to calculate the size of a specific chunk
     * @return chunk size in the order (datetime, y, x)
     */
    inline cube_size_tyx chunk_size() { return _chunk_size; }

    /**
     * @brief Get the spatiotemporal reference (extent, size, projection) of a cube
     * @return a cube_st_reference object
     */
    inline std::shared_ptr<cube_stref> st_reference() { return _st_ref; }

    /**
    * @brief Set the spatiotemporal reference (extent, size, projection) of a cube
     *
     * @note This is dangerous, as it will not update the reference of derived cubes automatically, use only before constructing derived cubes
     *
     * @param st_ref new cube_st_reference object (or from a derived class)
    */
    inline void st_reference(std::shared_ptr<cube_stref_regular> st_ref) {
        _st_ref = st_ref;
    }

    /**
     * @brief Get size of a cube
     * @return cube size / number of cells in the order (bands, datetime, y, x)
     */
    inline cube_size_btyx size() { return {_bands.count(), size_t(), size_y(), size_x()}; }

    /**
     * @brief Get the number of bands of the cube
     * @return  Integer number of cells
     */
    inline uint32_t size_bands() { return _bands.count(); }

    /**
     * @brief Get the number of cells in the temporal dimension
     * @return  Integer number of cells
     */
    inline uint32_t size_t() { return (_st_ref != nullptr) ? _st_ref->nt() : 0; }

    /**
    * @brief Get the number of cells in the y dimension
    * @return  Integer number of cells
    */
    inline uint32_t size_y() { return (_st_ref != nullptr) ? _st_ref->ny() : 0; }

    /**
    * @brief Get the number of cells in the x dimension
    * @return  Integer number of cells
    */
    inline uint32_t size_x() { return (_st_ref != nullptr) ? _st_ref->nx() : 0; }

    /**
     * @brief Read chunk data
     * Virtual function that is called to read actual data of a given chunk.
     * Please make sure that this function is called only when the data is really needed (lazy evaluation).
     *
     * @param id the id of the requested chunk
     * @return a smart pointer to chunk data
     */
    virtual std::shared_ptr<chunk_data> read_chunk(chunkid_t id) = 0;

    /**
     * @brief Write a data cube as a set of GeoTIFF files under a given directory
     *
     * @param dir directory where to store the files
     * @param p chunk processor instance, defaults to the global configuration
     */
    void write_chunks_gtiff(std::string dir,
                            std::shared_ptr<chunk_processor> p = config::instance()->get_default_chunk_processor());

    /**
     * @brief Writes a data cube as a collection of GeoTIFF files
     *
     * Each time slice of the cube will be written to one GeoTIFF file.

     * @param dir output path
     * @param prefix name prefix for created files
     * @param overviews generate GDAL overview images
     * @param cog write cloud-optmized GeoTIFFs (forces overviews = true)
     * @param overview_resampling resampling algorithm used to generate overviews (see https://gdal.org/programs/gdaladdo.html)
     * @param additional creation_options key value pairs passed to GDAL as GTiff creation options (see https://gdal.org/drivers/raster/gtiff.html)
     * @param packing reduce size of output tile with packing (apply scale + offset and use smaller integer data types)
     * @param drop_empty_slices if true, empty time slices will be skipped; not yet implemented
     * @param p chunk processor instance, defaults to the global configuration
     *
     * @note Overview levels will be chosen by halving the number of pixels until the larger dimension
     * has less than 256 pixels.
     *
     * @note argument `drop_empty_slices` is not yet implemented.
     *
     * @note Depending on `overviews` and `cog` GeoTIFFs created and postprocessed stepwise:
     * 1. time slices of cubes of the cube are exported as tiled normal GeoTIFFs
     * 2. Overviews are generated (internal)
     * 3. A new copy of the TIF with overviews is created, moving the IFDs of overviews to the beginning of the file (using gdal_translate -co COPY_SRC_OVERVIEWS=YES)
     * It seems not possible to generate COGs fomr in memory data without temporal copy, However,
     * improvements are very welcome (maybe the COG driver of GDAL > 3.1 helps?).
     */
    void write_tif_collection(std::string dir, std::string prefix = "",
                              bool overviews = false, bool cog = false,
                              std::map<std::string, std::string> creation_options = std::map<std::string, std::string>(),
                              std::string overview_resampling = "NEAREST",
                              packed_export packing = packed_export::make_none(),
                              bool drop_empty_slices = false,
                              std::shared_ptr<chunk_processor> p = config::instance()->get_default_chunk_processor());

    /**
     * Write a data cube as a single netCDF file
     * @param path path of the target file
     * @param compression_level deflate level, 0 = no compression, 1 = fast, 9 = small
     * @param with_VRT additional write VRT files for time slices that are easier to display, e.g., in QGIS
     * @param write_bounds boolean, if true, variables time_bnds, y_bnds, x_bnds per dimension values will be added
     * @param packing reduce size of output tile with packing (apply scale + offset and use smaller integer data types)
     * @param drop_empty_slices if true, empty time slices will be skipped
     * @param p chunk processor instance, defaults to the global configuration
     *
     * @note argument `drop_empty_slices` is not yet implemented.
     */
    void write_netcdf_file(std::string path, uint8_t compression_level = 0,
                           bool with_VRT = false, bool write_bounds = true, packed_export packing = packed_export::make_none(),
                           bool drop_empty_slices = false,
                           std::shared_ptr<chunk_processor> p = config::instance()->get_default_chunk_processor());

    void write_chunks_netcdf(std::string dir, std::string name = "", uint8_t compression_level = 0,
                             std::shared_ptr<chunk_processor> p = config::instance()->get_default_chunk_processor());

    std::shared_ptr<chunk_data> to_double_array(std::shared_ptr<chunk_processor> p = config::instance()->get_default_chunk_processor());

    void write_single_chunk_netcdf(chunkid_t id, std::string path, uint8_t compression_level = 0);

    /**
     * @brief Writes a data cube as a collection of PNG files
     *
     * Each time slice of the cube will be written to one GeoTIFF file.
     * @param dir output path
     * @param prefix name prefix for created files
     * @param bands vector with band names; must contain either no, one, or three elements; If empty, the first band is used to produce a grayscale image
     * @param zlim_min minimum color value per channel (if vector has size 1 but three bands are provided, the same value is used for all bands)
     * @param zlim_max maximum color value per channel (if vector has size 1 but three bands are provided, the same value is used for all bands)
     * @param gamma gamma correction factor to be applied
     * @param na_color = {}
     * @param na_transparent if true, the resulting image will cintain an alpha channel, where NA values will be assigned alpha = 0%
     * @param creation_options key value pairs passed to GDAL as PNG creation options (see https://gdal.org/drivers/raster/png.html)
     * @param drop_empty_slices if true, empty time slices will be skipped; not yet implemented
     * @param p chunk processor instance, defaults to the global configuration
     *
     * @note argument `drop_empty_slices` is not yet implemented.
     *
     * @note This function can be used to create single-band grayscale or three-bands color RGB plots.
     * If bands is empty, the first band of the data cube will be used to produce a single-band grayscale image. If
     * three bands are provided, they are interpreted in the order R,G,B. Depsite the number of provided bands,
     * the color scale defined by zlim_min and zlim_max may be provided as single element vectors in which case the values are used for all color channels, or as
     * three element vectors with separate scales for each RGB channel. If na_color is not empty, it can be used to set the color values of NA pixels either as a single gray value or as three separate RGB values.
     * In the case of producing a single band grayscale image but three separate vlaues as na_color, the output image will contain three channels with gray values (R=G=B) for non-NA pixels and the na_color for NA values.
     * Notice that if na_transparent is true, na_color is ignored.
     * Further creation options can be used, e.g., to control the compression level of PNG files.
     */
    void write_png_collection(std::string dir, std::string prefix = "",
                              std::vector<std::string> bands = {}, std::vector<double> zlim_min = {0}, std::vector<double> zlim_max = {255}, double gamma = 1, std::vector<uint8_t> na_color = {}, bool na_transparent = true,
                              std::map<std::string, std::string> creation_options = std::map<std::string, std::string>(),
                              bool drop_empty_slices = false,
                              std::shared_ptr<chunk_processor> p = config::instance()->get_default_chunk_processor());

    /**
     * Get the cube's bands
     * @return all bands of the cube object as band_collection
     */
    inline band_collection bands() {
        return _bands;
    }

    /**
     * @brief Add a data cube as a parent cube to keep track of cube connections
     * @param c the data cube this cube depends on
     */
    inline void add_parent_cube(std::shared_ptr<cube> c) {
        _pre.push_back(std::weak_ptr<cube>(c));
    }

    /**
     * @brief Add a child data cube to keep track of cubes connections
     * @param c derived data cube
     */
    inline void add_child_cube(std::shared_ptr<cube> c) {
        _succ.push_back(std::weak_ptr<cube>(c));
    }

    /**
     * Abstract function to create a JSON representation of a cube object
     * @return a JSON object which can be used to recreate it with cube_factory
     * @see cube_factory
     */
    virtual json11::Json make_constructible_json() = 0;

   protected:
    /**
     * Spacetime reference of a cube, including extent, size, and projection
     */
    std::shared_ptr<cube_stref> _st_ref;  // TODO: replace with unique_ptr

    /**
     * @brief Size of chunks / number of cells per chunk in the order (datetime, y, x)
     */
    cube_size_tyx _chunk_size;

    /**
     * @brief The collection's bands
     */
    band_collection _bands;

    /**
     * @brief List of cube instances this cube takes as input (predecessor)
     */
    std::vector<std::weak_ptr<cube>> _pre;

    /**
     * @brief List of cube instances that take this object as input (successors)
     */
    std::vector<std::weak_ptr<cube>> _succ;


    /**
     * @brief Optimization flag defining whether pixel values depend on pixel values located somewhere else
     */
    bool _optim_is_local; // TODO
    bool _optim_is_chunk_local;  // TODO
    bool _optim_reccomend_cache_input; // TODO

};

}  // namespace gdalcubes

#endif  //CUBE_H
