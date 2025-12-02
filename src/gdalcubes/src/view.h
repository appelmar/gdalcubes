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

#ifndef VIEW_H
#define VIEW_H


#include <cstdint> // 2023-01-12: GCC 13 compatibility
#include "coord_types.h"
//#include "datetime.h"
#include "external/json11/json11.hpp"



namespace gdalcubes {

/**
 * A utility structure to work with different aggregation
 * algorithms
 */
struct aggregation {
    enum class aggregation_type {
        AGG_NONE,
        AGG_MIN,
        AGG_MAX,
        AGG_MEAN,
        AGG_MEDIAN,
        AGG_FIRST,
        AGG_LAST,
        AGG_IMAGE_COUNT,
        AGG_VALUE_COUNT
    };

    static aggregation_type from_string(std::string s);

    static std::string to_string(aggregation_type a);
};

/**
 * @brief Utility structure to work with different resampling
 * algorithms and their different types in GDAL
 */
struct resampling {
    /**
     * @brief An enumeration listing all available resampling types
     */
    enum class resampling_type {
        RSMPL_NEAR,
        RSMPL_BILINEAR,
        RSMPL_CUBIC,
        RSMPL_CUBICSPLINE,
        RSMPL_LANCZOS,
        RSMPL_AVERAGE,
        RSMPL_MODE,
        RSMPL_MAX,
        RSMPL_MIN,
        RSMPL_MED,
        RSMPL_Q1,
        RSMPL_Q3
    };

    /**
     * @brief Get the resampling type from its name
     * @param s string representation of a resampling algorithm
     * @return the corresponding resampling type entry
     */
    static resampling_type from_string(std::string s);

    /**
     * @brief Get the name of a resampling method
     * @param r resampling type
     * @return name string of the resampling method
     */
    static std::string to_string(resampling_type r);

    /**
     * @brief Convert a resampling type to the corresponding GDAL type to be used in RasterIO
     * @note RasterIO does not support all available resampling types
     * @param r resampling type
     * @return GDALRIOResampleAlg
     */
    static GDALRIOResampleAlg to_gdal_rasterio(resampling_type r);
};

// Forwartds declarations of all spacetime reference classes
class cube_stref_regular;
class cube_stref_labeled_time;

struct duration;
class datetime;
enum class datetime_unit;

class cube_stref {
   public:
    virtual uint32_t nx() = 0;
    virtual uint32_t ny() = 0;
    virtual uint32_t nt() = 0;
    virtual double left() = 0;
    virtual double right() = 0;
    virtual double bottom() = 0;
    virtual double top() = 0;
    virtual datetime t0() = 0;
    virtual datetime t1() = 0;
    virtual std::string srs() = 0;
    virtual OGRSpatialReference srs_ogr() const = 0;

    // Spatiotemporal support of grid cells
    virtual duration dt() = 0;
    virtual datetime_unit dt_unit() = 0;
    virtual int32_t dt_interval() = 0;
    virtual double dx() = 0;
    virtual double dy() = 0;

    virtual coords_st map_coords(coords_nd<uint32_t, 3> p) = 0;
    virtual coords_nd<uint32_t, 3> cube_coords(coords_st p) = 0;
    virtual datetime datetime_at_index(uint32_t index) = 0;
    virtual uint32_t index_at_datetime(datetime t) = 0;

    virtual std::shared_ptr<cube_stref> copy() = 0;

    virtual bool has_regular_space() = 0;
    virtual bool has_regular_time() = 0;

    // TODO: add x_at_index, y_at_index, index_at_x, index_at_y

    // TODO: add dimension label functions x_labels
    //    virtual std::vector<double> x_labels();
    //    virtual std::vector<double> y_labels();
    //    virtual std::vector<datetime> t_labels();

    // TODO: add subset function?
    //virtual std::unique_ptr<cube_stref> subset_t() ... // TODO

    //virtual std::unique_ptr<cube_stref> copy() = 0;

    static std::string type_string(std::shared_ptr<cube_stref> obj);
};

/**
 * @brief Spatial and temporal reference for data cubes
 */
class cube_stref_regular : public cube_stref {
   public:
    virtual ~cube_stref_regular() {}

    virtual void set_x_axis(double min, double max, double delta);
    virtual void set_x_axis(double min, double max, uint32_t n);
    virtual void set_x_axis(double min, uint32_t n, double delta);
    virtual void set_x_axis(uint32_t n, double max, double delta);




    virtual void set_y_axis(double min, double max, double delta);
    virtual void set_y_axis(double min, double max, uint32_t n) ;
    virtual void set_y_axis(double min, uint32_t n, double delta);
    virtual void set_y_axis(uint32_t n, double max, double delta);


    virtual void set_t_axis(datetime min, datetime max, duration delta);


    virtual void set_t_axis(datetime min, datetime max, uint32_t n);
    virtual uint32_t nx() override { return _nx; }
    virtual uint32_t ny() override { return _ny; }
    virtual double dx() override { return (_win.right - _win.left) / _nx; }
    virtual bool has_regular_space() override {return true;}
    virtual bool has_regular_time() override {return true; }
    virtual double dy() override { return (_win.top - _win.bottom) / _ny; }
    virtual double left() override { return _win.left; }
    virtual double right() override { return _win.right; }
    virtual double bottom() override { return _win.bottom; }
    virtual double top() override { return _win.top; }
    virtual std::string srs() override { return _srs; }
    virtual void srs(std::string srs) { _srs = srs; }

    /**
     * Return the spatial reference system / projection
     * @return OGRSpatialReference object
     */
    virtual OGRSpatialReference srs_ogr() const override;
    

    /**
     * Getter / setter for the lower boundary of the cube's temporal extent (start datetime)
     * @return reference to the object's t0 object
     */
    virtual datetime t0() override { return _t0; }

    virtual datetime t1() override { return _t1; }

    virtual uint32_t nt() override;


    virtual bounds_2d<double> win() { return _win; }
    virtual void win(bounds_2d<double> win) { _win = win; }

    virtual duration dt() override { return _dt; }
    virtual datetime_unit dt_unit() override { return _dt.dt_unit; }
    virtual void dt_unit(datetime_unit unit) { _dt.dt_unit = unit; }
    virtual int32_t dt_interval() override { return _dt.dt_interval; }
    virtual void dt_interval(int32_t interval) { _dt.dt_interval = interval; }


    /**
    * Convert integer cube-based coordinates to spacetime coordinates
    * @note cube-based coordinates are in the order (t,y,x), (0,0,0) corresponds to the earliest date (t0) for the
    * upper left pixel.
    * @note Output coordinates will have the projection / SRS as in cube_st_reference::proj()
    * @see cube_st_reference::view_coords()
    * @param p cube-based coordinates
    * @return spacetime coordinates
    */
    virtual coords_st map_coords(coords_nd<uint32_t, 3> p) override;

    /**
     * Convert spacetime coordinates to integer cube-based coordinates
     * @note cube-based coordinates are in the order (t,y,x), (0,0,0) corresponds to the earliest date (t0) for the
     * lower left pixel.
     * @note the function assumes input coordinates have the  projection / SRS as in cube_st_reference::proj()
     * @see cube_st_reference::map_coords()
     * @param p spacetime coordinates
     * @return cube-based coordinates
     */
    virtual coords_nd<uint32_t, 3> cube_coords(coords_st p) override ;

    virtual datetime datetime_at_index(uint32_t index) override {
        return _t0 + _dt * index;
    }

    virtual uint32_t index_at_datetime(datetime t) override {
        return (uint32_t)((t - _t0) / _dt);
    }

    virtual std::shared_ptr<cube_stref> copy() override;

    friend bool operator==(const cube_stref_regular &l, const cube_stref_regular &r);

    inline friend bool operator!=(const cube_stref_regular &l, const cube_stref_regular &r) { return !(l == r); }

   protected:
    /**
     * @brief Spatial reference system / projection
     *
     * The string must be readable for OGRSpatialReference::SetFromUserInput, i.e.,
     * it can be "EPSG:xxx", WKT, or PROJ.4
     *
     */
    std::string _srs;

    /**
     * @brief Spatial window
     */
    bounds_2d<double> _win;

    datetime _t0;
    datetime _t1;
    uint32_t _nx;
    uint32_t _ny;

    duration _dt;
};

/**
 * A data cube view includes the spacetime reference of a cube (extent, resolution, projection) and
 * optional resampling and aggregation algorithms that are applied when original images are
 * read from an image_collection_cube. Aggregation refers to how multiple values for the same
 * cube cell from different images are combined whereas resampling refers to the algorithm used
 * to warp / reproject images to the cube geometry.
 */
class cube_view : public cube_stref_regular {
   public:
    cube_view() : _resampling(resampling::resampling_type::RSMPL_NEAR), _aggregation(aggregation::aggregation_type::AGG_FIRST) {}
    /**
         * Deserializes a cube_view object from a JSON file.
         * @param filename Path to the json file on disk
         * @return A cube_view object
         */
    static cube_view read_json(std::string filename);

    /**
        * Deserializes a cube_view object from a JSON file.
        * @param str JSON string
        * @return A cube_view object
        */
    static cube_view read_json_string(std::string str);

    /**
        * Serializes a cube_view object as a JSON file.
        * @param filename output file
        */
    void write_json(std::string filename);

    /**
          * Serializes a cube_view object as a JSON string.
          * @return JSON string
          */
    std::string write_json_string();

    /**
         * Getter / setter for aggregation method
         * @return reference to the object's aggregation field
         */
    inline aggregation::aggregation_type &aggregation_method() { return _aggregation; }

    /**
        * Getter / setter for resampling method
        * @return reference to the object's resampling field
        */
    inline resampling::resampling_type &resampling_method() { return _resampling; }

   private:
    static cube_view read(json11::Json j);

    resampling::resampling_type _resampling;
    aggregation::aggregation_type _aggregation;
};

class cube_stref_labeled_time : public cube_stref_regular {
   public:
    virtual bool has_regular_space() override {
        return true;
    }

    virtual bool has_regular_time() override {
        return false;
    }

//    virtual void t0(datetime t0) override {
//        // DO NOTHING, t0 is derived from labels
//    }

    virtual datetime t0() override {
        return _t_values[0];
    }

//    virtual void t1(datetime t1) override {
//        // DO NOTHING, t0 is derived from labels
//    }

    virtual datetime t1() override {
        return _t_values[_t_values.size() - 1];
    }

//    virtual void nt(uint32_t n) override {
//        // DO NOTHING, nt is derived from labels
//    }

    virtual uint32_t nt() override {
        return _t_values.size();
    }

//    virtual void dt(duration dt) override {
//        _t0.unit(dt.dt_unit);
//        _t1.unit(dt.dt_unit);
//        _dt = dt;
//    }

    virtual duration dt() override {
        return _dt;
    }

    void set_time_labels(std::vector<datetime> t);
    void set_time_labels(std::vector<std::string> t);

    std::vector<datetime> get_time_labels() {
        return _t_values;
    }

    std::vector<std::string> get_time_labels_as_string();
    virtual coords_st map_coords(coords_nd<uint32_t, 3> p) override ;
    virtual coords_nd<uint32_t, 3> cube_coords(coords_st p) override;

    virtual datetime datetime_at_index(uint32_t index) override {
        return _t_values[index];
    }

    virtual uint32_t index_at_datetime(datetime t) override;

    friend bool operator==(const cube_stref_labeled_time &l, const cube_stref_labeled_time &r);

    inline friend bool operator!=(const cube_stref_labeled_time &l, const cube_stref_labeled_time &r) { return !(l == r); }

    virtual std::shared_ptr<cube_stref> copy() override;

   protected:
    std::vector<datetime> _t_values;
    std::map<datetime, uint32_t> _t_index;
};

}  // namespace gdalcubes

#endif  //VIEW_H
