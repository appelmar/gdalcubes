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

#ifndef VIEW_H
#define VIEW_H

#include <algorithm>
#include <cstdint> // 2023-01-12: GCC 13 compatibility
#include "coord_types.h"
#include "datetime.h"
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

    static aggregation_type from_string(std::string s) {
        std::transform(s.begin(), s.end(), s.begin(), ::tolower);
        if (s == "none") {
            return aggregation_type::AGG_NONE;
        } else if (s == "min") {
            return aggregation_type::AGG_MIN;
        } else if (s == "max") {
            return aggregation_type::AGG_MAX;
        } else if (s == "mean") {
            return aggregation_type::AGG_MEAN;
        } else if (s == "median") {
            return aggregation_type::AGG_MEDIAN;
        } else if (s == "first") {
            return aggregation_type::AGG_FIRST;
        } else if (s == "last") {
            return aggregation_type::AGG_LAST;
        } else if (s == "count_images") {
            return aggregation_type::AGG_IMAGE_COUNT;
        } else if (s == "count_values") {
            return aggregation_type::AGG_VALUE_COUNT;
        }
        return aggregation_type::AGG_NONE;
    }

    static std::string to_string(aggregation_type a) {
        switch (a) {
            case aggregation_type::AGG_NONE:
                return "none";
            case aggregation_type::AGG_MIN:
                return "min";
            case aggregation_type::AGG_MAX:
                return "max";
            case aggregation_type::AGG_MEAN:
                return "mean";
            case aggregation_type::AGG_MEDIAN:
                return "median";
            case aggregation_type::AGG_FIRST:
                return "first";
            case aggregation_type::AGG_LAST:
                return "last";
            case aggregation_type::AGG_IMAGE_COUNT:
                return "count_images";
            case aggregation_type::AGG_VALUE_COUNT:
                return "count_values";
            default:
                return "none";
        }
    }
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
    static resampling_type from_string(std::string s) {
        std::transform(s.begin(), s.end(), s.begin(), ::tolower);
        if (s == "near" || s == "nearest") {
            return resampling_type::RSMPL_NEAR;
        } else if (s == "bilinear") {
            return resampling_type::RSMPL_BILINEAR;
        } else if (s == "cubic") {
            return resampling_type::RSMPL_CUBIC;
        } else if (s == "cubicspline") {
            return resampling_type::RSMPL_CUBICSPLINE;
        } else if (s == "lanczos") {
            return resampling_type::RSMPL_LANCZOS;
        } else if (s == "average" || s == "mean") {
            return resampling_type::RSMPL_AVERAGE;
        } else if (s == "mode") {
            return resampling_type::RSMPL_MODE;
        } else if (s == "max") {
            return resampling_type::RSMPL_MAX;
        } else if (s == "min") {
            return resampling_type::RSMPL_MIN;
        } else if (s == "med" || s == "median") {
            return resampling_type::RSMPL_MED;
        } else if (s == "q1") {
            return resampling_type::RSMPL_Q1;
        } else if (s == "q3") {
            return resampling_type::RSMPL_Q3;
        }
        return resampling_type::RSMPL_NEAR;
    }

    /**
         * @brief Get the name of a resampling method
         * @param r resampling type
         * @return name string of the resampling method
         */
    static std::string to_string(resampling_type r) {
        switch (r) {
            case resampling_type::RSMPL_NEAR:
                return "near";
            case resampling_type::RSMPL_BILINEAR:
                return "bilinear";
            case resampling_type::RSMPL_CUBIC:
                return "cubic";
            case resampling_type::RSMPL_CUBICSPLINE:
                return "cubicspline";
            case resampling_type::RSMPL_LANCZOS:
                return "lanczos";
            case resampling_type::RSMPL_AVERAGE:
                return "average";
            case resampling_type::RSMPL_MODE:
                return "mode";
            case resampling_type::RSMPL_MAX:
                return "max";
            case resampling_type::RSMPL_MIN:
                return "min";
            case resampling_type::RSMPL_MED:
                return "med";
            case resampling_type::RSMPL_Q1:
                return "q1";
            case resampling_type::RSMPL_Q3:
                return "q3";
            default:
                return "near";
        }
    }

    /**
         * @brief Convert a resampling type to the corresponding GDAL type to be used in RasterIO
         * @note RasterIO does not support all available resampling types
         * @param r resampling type
         * @return GDALRIOResampleAlg
         */
    static GDALRIOResampleAlg to_gdal_rasterio(resampling_type r) {
        switch (r) {
            case resampling_type::RSMPL_NEAR:
                return GRIORA_NearestNeighbour;
            case resampling_type::RSMPL_BILINEAR:
                return GRIORA_Bilinear;
            case resampling_type::RSMPL_CUBIC:
                return GRIORA_Cubic;
            case resampling_type::RSMPL_CUBICSPLINE:
                return GRIORA_CubicSpline;
            case resampling_type::RSMPL_LANCZOS:
                return GRIORA_Lanczos;
            case resampling_type::RSMPL_AVERAGE:
                return GRIORA_Average;
            case resampling_type::RSMPL_MODE:
                return GRIORA_Mode;
            case resampling_type::RSMPL_MAX:
            case resampling_type::RSMPL_MIN:
            case resampling_type::RSMPL_MED:
            case resampling_type::RSMPL_Q1:
            case resampling_type::RSMPL_Q3:
            default:
                return GRIORA_NearestNeighbour;  // Not yet defined in gdal.h
        }
    }
};

// Forwartds declarations of all spacetime reference classes
class cube_stref_regular;
class cube_stref_labeled_time;

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

    static std::string type_string(std::shared_ptr<cube_stref> obj) {
        if (std::dynamic_pointer_cast<cube_stref_labeled_time>(obj) != nullptr) {
            return "cube_stref_labeled_time";
        }
        if (std::dynamic_pointer_cast<cube_stref_regular>(obj) != nullptr) {
            return "cube_stref_regular";
        }
        return "";
    }
};

/**
 * @brief Spatial and temporal reference for data cubes
 */
class cube_stref_regular : public cube_stref {
   public:
    virtual ~cube_stref_regular() {}

    virtual void set_x_axis(double min, double max, double delta) {
        _win.left = min;
        _win.right = max;
        _nx = (uint32_t)std::ceil((_win.right - _win.left) / delta);
        double exp_x = _nx * delta - (_win.right - _win.left);
        _win.right += exp_x / 2;
        _win.left -= exp_x / 2;
        if (std::fabs(exp_x) > std::numeric_limits<double>::epsilon()) {
            GCBS_INFO("Size of the cube in x direction does not align with dx, extent will be enlarged by " +
                      std::to_string(exp_x / 2) + " at both sides.");
        }
    }
    virtual void set_x_axis(double min, double max, uint32_t n) {
        _win.left = min;
        _win.right = max;
        _nx = n;
    }
    virtual void set_x_axis(double min, uint32_t n, double delta) {
        _win.left = min;
        _nx = n;
        _win.right = min + n * delta;
    }
    virtual void set_x_axis(uint32_t n, double max, double delta) {
        _win.right = max;
        _nx = n;
        _win.left = max - n * delta;
    }




    virtual void set_y_axis(double min, double max, double delta) {
        _win.bottom= min;
        _win.top = max;
        _ny = (uint32_t)std::ceil((_win.top - _win.bottom) / delta);
        double exp_y = _ny * delta - (_win.top - _win.bottom);
        _win.top += exp_y / 2;
        _win.bottom -= exp_y / 2;
        if (std::fabs(exp_y) > std::numeric_limits<double>::epsilon()) {
            GCBS_INFO("Size of the cube in y direction does not align with dy, extent will be enlarged by " +
                      std::to_string(exp_y / 2) + " at both sides.");
        }
    }
    virtual void set_y_axis(double min, double max, uint32_t n) {
        _win.bottom = min;
        _win.top = max;
        _ny = n;
    }
    virtual void set_y_axis(double min, uint32_t n, double delta) {
        _win.bottom = min;
        _ny = n;
        _win.top = min + n * delta;

    }
    virtual void set_y_axis(uint32_t n, double max, double delta) {
        _win.top = max;
        _ny = n;
        _win.bottom = max - n * delta;
    }

    virtual void set_t_axis(datetime min, datetime max, duration delta) {
        if (min.unit() != max.unit()) { // If different units, use coarser
            min.unit(std::max(min.unit(), max.unit()));
            max.unit(std::max(min.unit(), max.unit()));
            GCBS_WARN("Different datetime units given for start / end dates, using " + datetime::unit_to_string(min.unit()) );
        }
        datetime_unit tu = min.unit();
        datetime_unit u = delta.dt_unit;

        if (tu > u) {
            // This part makes sure that if t1 and t0 have a coarser datetime unit than dt,
            // t0 and t1 are extendec accordingly.
            // For example, t0 = "2000", t1 = "2000", dt = "P1D" would lead to t0 = "2000-01-01", t1 = "2000-12-31"
            date::sys_seconds min_new =
                date::sys_days{date::year(min.year()) /
                    ((tu <= datetime_unit::MONTH)? date::month(min.month()) : date::month(1)) /
                    ((tu <= datetime_unit::DAY)? date::day(min.dayofmonth()) : date::day(1))} +
                    ((tu <= datetime_unit::HOUR)? std::chrono::hours{min.hours()}  :std::chrono::hours{0}) +
                    ((tu <= datetime_unit::MINUTE)? std::chrono::minutes{min.minutes()}  :std::chrono::minutes{0}) +
                    ((tu <= datetime_unit::SECOND)? std::chrono::seconds{min.seconds()}  :std::chrono::seconds{0});

            date::sys_seconds max_new;

            if (tu > datetime_unit::DAY) {
                max_new = date::sys_days{date::year(max.year()) /
                               ((tu <= datetime_unit::MONTH)? date::month(max.month()) : date::month(12)) /
                               date::last} + std::chrono::hours{23} + std::chrono::minutes{59} + std::chrono::seconds{59};
            }
            else {
                max_new =
                    date::sys_days{date::year(max.year()) / date::month(max.month())  / date::day(max.dayofmonth())} +
                    ((tu <= datetime_unit::HOUR)? std::chrono::hours{max.hours()}  :std::chrono::hours{23}) +
                    ((tu <= datetime_unit::MINUTE)? std::chrono::minutes{max.minutes()}  :std::chrono::minutes{59}) +
                    ((tu <= datetime_unit::SECOND)? std::chrono::seconds{max.seconds()}  :std::chrono::seconds{59});
            }
            _t0 = datetime(min_new, u);
            _t1 = datetime(max_new, u);
        }
        else {
            _t0 = min;
            _t1 = max;
            _t0.unit(u);
            _t1.unit(u);
        }


        duration dtotal = _t1 - _t0;  // + 1 if include end date2
        dtotal.dt_interval += 1;
        if (dtotal % delta != 0) {
            duration end_duration;  // end duration has one (day / month / unit) less than dt)
            end_duration.dt_interval = delta.dt_interval - 1;
            end_duration.dt_unit = u;
            _t1 = (_t0 + delta * (dtotal / delta)) + end_duration;
        }
        _dt = delta;

        if (u == datetime_unit::YEAR) {
            auto p0 = date::sys_days{date::year(_t0.year()) / date::month(1) / date::day(1)} +
                      std::chrono::hours{0} + std::chrono::minutes{0} + std::chrono::seconds{0};
            _t0 = datetime(p0, u);

            auto p1 = date::sys_days{date::year(_t1.year()) / date::month(12) / date::day(31)} +
                      std::chrono::hours{23} + std::chrono::minutes{59} + std::chrono::seconds{59};
            _t1 = datetime(p1, u);
        }
        else if (u == datetime_unit::MONTH) {
            auto p0 = date::sys_days{date::year(_t0.year()) / date::month(_t0.month()) / date::day(1)} +
                      std::chrono::hours{0} + std::chrono::minutes{0} + std::chrono::seconds{0};
            _t0 = datetime(p0, u);

            auto p1 = date::sys_days{date::year(_t1.year()) / date::month(_t1.month()) / date::last} +
                      std::chrono::hours{23} + std::chrono::minutes{59} + std::chrono::seconds{59};
            _t1 = datetime(p1, u);
        }
        else if (u == datetime_unit::DAY) {
            auto p0 = date::sys_days{date::year(_t0.year()) / date::month(_t0.month()) / date::day(_t0.dayofmonth())} +
                      std::chrono::hours{0} + std::chrono::minutes{0} + std::chrono::seconds{0};
            _t0 = datetime(p0, u);

            auto p1 = date::sys_days{date::year(_t1.year()) / date::month(_t1.month()) / date::day(_t1.dayofmonth())} +
                      std::chrono::hours{23} + std::chrono::minutes{59} + std::chrono::seconds{59};
            _t1 = datetime(p1, u);
        }
        else if (u == datetime_unit::HOUR) {
            auto p0 = date::sys_days{date::year(_t0.year()) / date::month(_t0.month()) / date::day(_t0.dayofmonth())} +
                      std::chrono::hours{_t0.hours()} + std::chrono::minutes{0} + std::chrono::seconds{0};
            _t0 = datetime(p0, u);

            auto p1 = date::sys_days{date::year(_t1.year()) / date::month(_t1.month()) / date::day(_t1.dayofmonth())} +
                      std::chrono::hours{_t1.hours()} + std::chrono::minutes{59} + std::chrono::seconds{59};
            _t1 = datetime(p1, u);
        }
        else if (u == datetime_unit::MINUTE) {
            auto p0 = date::sys_days{date::year(_t0.year()) / date::month(_t0.month()) / date::day(_t0.dayofmonth())} +
                      std::chrono::hours{_t0.hours()} + std::chrono::minutes{_t0.minutes()} + std::chrono::seconds{0};
            _t0 = datetime(p0, u);

            auto p1 = date::sys_days{date::year(_t1.year()) / date::month(_t1.month()) / date::day(_t1.dayofmonth())} +
                      std::chrono::hours{_t1.hours()} + std::chrono::minutes{_t1.minutes()} + std::chrono::seconds{59};
            _t1 = datetime(p1, u);
        }
        else if (u == datetime_unit::SECOND) { // not needed because SECOND is lowest
            auto p0 = date::sys_days{date::year(_t0.year()) / date::month(_t0.month()) / date::day(_t0.dayofmonth())} +
                      std::chrono::hours{_t0.hours()} + std::chrono::minutes{_t0.minutes()} + std::chrono::seconds{_t0.seconds()};
            _t0 = datetime(p0, u);

            auto p1 = date::sys_days{date::year(_t1.year()) / date::month(_t1.month()) / date::day(_t1.dayofmonth())} +
                      std::chrono::hours{_t1.hours()} + std::chrono::minutes{_t1.minutes()} + std::chrono::seconds{_t1.seconds()};
            _t1 = datetime(p1, u);
        }


        // check whether min/max have been modified (string conversion is somewhat ugly here)
        std::string t0str = _t0.to_string();
        std::string t1str = _t1.to_string();
        if (t0str != min.to_string() ||
            t1str != max.to_string()) {
            GCBS_INFO("Temporal extent of the cube does not align with dt and has been extended to " + _t0.to_string() + "/" + _t1.to_string());
        }

    }


    virtual void set_t_axis(datetime min, datetime max, uint32_t n) {
        if (min.unit() != max.unit()) { // If different units, use coarser
            min.unit(std::max(min.unit(), max.unit()));
            max.unit(std::max(min.unit(), max.unit()));
            GCBS_WARN("Different datetime units given for start / end dates, using " + datetime::unit_to_string(min.unit()) );
        }
        _t0 = min;
        _t1 = max;
        duration d = (_t1 - _t0) + 1;
        duration dnew = dt();
        if (dnew.dt_interval == 0) {  // if dt has not been set
            dnew.dt_unit = d.dt_unit;
            // alternatively, a "reasonable" value should be derived here
        }
        dnew.dt_interval = (int32_t)std::ceil((double)d.dt_interval / (double)n);
        _dt = dnew;
        if (d.dt_interval % n != 0) {
            _t1 = _t0 + _dt * (n - 1);
        }
    }
//    virtual void set_t_axis(double min, uint32_t n, duration delta) {
//        // NOT YET IMPLEMENTED
//    }
//    virtual void set_t_axis(uint32_t n, double max, duration delta) {
//        // NOT YET IMPLEMENTED
//    }

    virtual uint32_t nx() override { return _nx; }
    //virtual void nx(uint32_t nx) { _nx = nx; }

    virtual uint32_t ny() override { return _ny; }
    //virtual void ny(uint32_t ny) { _ny = ny; }

    virtual double dx() override { return (_win.right - _win.left) / _nx; }

    virtual bool has_regular_space() override {
        return true;
    }

    virtual bool has_regular_time() override {
        return true;
    }

    /**
    * Set the size of cells in x dimension
    * @note if the width of the spatial window is not a multiple of the new dx, the window will be widened at both ends
    * @param x size of cells in x dimension
    */
//    virtual void dx(double dx) {
//        _nx = (uint32_t)std::ceil((_win.right - _win.left) / dx);
//        double exp_x = _nx * dx - (_win.right - _win.left);
//        _win.right += exp_x / 2;
//        _win.left -= exp_x / 2;
//        if (std::fabs(exp_x) > std::numeric_limits<double>::epsilon()) {
//            GCBS_INFO("Size of the cube in x direction does not align with dx, extent will be enlarged by " +
//                      std::to_string(exp_x / 2) + " at both sides.");
//        }
//    }

    virtual double dy() override { return (_win.top - _win.bottom) / _ny; }
//    virtual void dy(double dy) {
//        _ny = (uint32_t)std::ceil((_win.top - _win.bottom) / dy);
//        double exp_y = _ny * dy - (_win.top - _win.bottom);
//        _win.top += exp_y / 2;
//        _win.bottom -= exp_y / 2;
//        if (std::fabs(exp_y) > std::numeric_limits<double>::epsilon()) {
//            GCBS_INFO("Size of the cube in y direction does not align with dy, extent will be enlarged by " +
//                      std::to_string(exp_y / 2) + " at both sides.");
//        }
//    }

    virtual double left() override { return _win.left; }
    //virtual void left(double left) { _win.left = left; }

    virtual double right() override { return _win.right; }
    //virtual void right(double right) { _win.right = right; }

    virtual double bottom() override { return _win.bottom; }
    //virtual void bottom(double bottom) { _win.bottom = bottom; }

    virtual double top() override { return _win.top; }
   // virtual void top(double top) { _win.top = top; }

    virtual std::string srs() override { return _srs; }
    virtual void srs(std::string srs) { _srs = srs; }

    /**
         * Return the spatial reference system / projection
         * @return OGRSpatialReference object
         */
    virtual OGRSpatialReference srs_ogr() const override {
        OGRSpatialReference s;
        s.SetFromUserInput(_srs.c_str());
        return s;
    }

    /**
         * Getter / setter for the lower boundary of the cube's temporal extent (start datetime)
         * @return reference to the object's t0 object
         */
    virtual datetime t0() override { return _t0; }
    //virtual void t0(datetime t0) { _t0 = t0; }

    virtual datetime t1() override { return _t1; }
    //virtual void t1(datetime t1) { _t1 = t1; }

    virtual uint32_t nt() override {
        if (_t1 == _t0) return 1;
        duration d = (_t1 - _t0) + 1;
        return (d % _dt == 0) ? d / _dt : (1 + (d / _dt));
    }

//    virtual void nt(uint32_t n) {
//        duration d = (_t1 - _t0) + 1;
//        duration dnew = dt();
//        if (dnew.dt_interval == 0) {  // if dt has not been set
//            dnew.dt_unit = d.dt_unit;
//            // alternatively, a "reasonable" should be derived here
//        }
//        dnew.dt_interval = (int32_t)std::ceil((double)d.dt_interval / (double)n);
//        _dt = dnew;
//        if (d.dt_interval % n != 0) {
//            _t1 = _t0 + _dt * (n - 1);
//            GCBS_INFO(
//                "Temporal size of the cube does not align with nt, end date/time of the cube will be extended to " +
//                _t1.to_string() + ".");
//        }
//        //
//        //        if (nt() == n - 1) {  // in some cases (e.g. d == 9M, n==4), we must extend the temporal extent of the view
//        //            _t1 = _t1 + dt();
//        //            GCBS_WARN("Extent in t direction is indivisible by nt, end date/time will be set to " + _t1.to_string());
//        //        }
//        assert(nt() == n);
//    }

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
    virtual coords_st map_coords(coords_nd<uint32_t, 3> p) override {
        coords_st s;
        s.s.x = _win.left + p[2] * dx();
        s.s.y = _win.top - p[1] * dy();
        s.t = _t0 + _dt * p[0];
        return s;
    }

    /**
         * Convert spacetime coordinates to integer cube-based coordinates
         * @note cube-based coordinates are in the order (t,y,x), (0,0,0) corresponds to the earliest date (t0) for the
         * lower left pixel.
         * @note the function assumes input coordinates have the  projection / SRS as in cube_st_reference::proj()
         * @see cube_st_reference::map_coords()
         * @param p spacetime coordinates
         * @return cube-based coordinates
         */
    virtual coords_nd<uint32_t, 3> cube_coords(coords_st p) override {
        coords_nd<uint32_t, 3> s;
        s[2] = (uint32_t)((p.s.x - _win.left) / dx());
        s[1] = (uint32_t)((_win.top - p.s.y) / dy());
        s[0] = (uint32_t)((p.t - _t0) / _dt);
        return s;
    }

    virtual datetime datetime_at_index(uint32_t index) override {
        return _t0 + _dt * index;
    }

    virtual uint32_t index_at_datetime(datetime t) override {
        return (uint32_t)((t - _t0) / _dt);
    }

    virtual std::shared_ptr<cube_stref> copy() override {
        std::shared_ptr<cube_stref_regular> x = std::make_shared<cube_stref_regular>();
        x->_win.left = _win.left;
        x->_win.right = _win.right;
        x->_win.top = _win.top;
        x->_win.bottom = _win.bottom;
        x->_nx = _nx;
        x->_ny = _ny;
        x->_dt = _dt;
        x->_t0 = _t0;
        x->_t1 = _t1;
        x->_srs = _srs;
        return x;
    }

    friend bool operator==(const cube_stref_regular &l, const cube_stref_regular &r) {
        if (!(l._win.left == r._win.left &&
              l._win.right == r._win.right &&
              l._win.top == r._win.top &&
              l._win.bottom == r._win.bottom &&
              l._nx == r._nx &&
              l._ny == r._ny &&
              l._t0 == r._t0 &&
              l._dt == r._dt))
            return false;

        // compare SRS
        OGRSpatialReference a = l.srs_ogr();
        OGRSpatialReference b = r.srs_ogr();

        if (!a.IsSame(&b))
            return false;

        return true;
    }

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

    void set_time_labels(std::vector<datetime> t) {
        // TODO: if t.empty()
        _t_values = t;

        // TODO: what if not sorted?
        for (uint32_t i = 0; i < _t_values.size(); ++i) {
            _t_index.insert(std::make_pair(t[i], i));
        }
    }

    void set_time_labels(std::vector<std::string> t) {
        // TODO: what if not sorted?
        for (uint32_t i = 0; i < _t_values.size(); ++i) {
            _t_index.insert(std::make_pair(datetime::from_string(t[i]), i));
        }
    }

    std::vector<datetime> get_time_labels() {
        return _t_values;
    }

    std::vector<std::string> get_time_labels_as_string() {
        std::vector<std::string> out;
        for (uint32_t i = 0; i < _t_values.size(); ++i) {
            out.push_back(_t_values[i].to_string());
        }
        return out;
    }

    virtual coords_st map_coords(coords_nd<uint32_t, 3> p) override {
        coords_st s;
        s.s.x = _win.left + p[2] * dx();
        s.s.y = _win.top - p[1] * dy();
        s.t = datetime_at_index(p[0]);
        return s;
    }

    virtual coords_nd<uint32_t, 3> cube_coords(coords_st p) override {
        coords_nd<uint32_t, 3> s;
        s[2] = (uint32_t)((p.s.x - _win.left) / dx());
        s[1] = (uint32_t)((_win.top - p.s.y) / dy());
        s[0] = index_at_datetime(p.t);
        return s;
    }

    virtual datetime datetime_at_index(uint32_t index) override {
        return _t_values[index];
    }

    virtual uint32_t index_at_datetime(datetime t) override {
        auto res = _t_index.find(t);
        if (res == _t_index.end()) {
            GCBS_ERROR("Data cubes does not contain time slice for requested datetime");
            throw std::string("Data cubes does not contain time slice for requested datetime");
        }
        return res->second;
    }

    friend bool operator==(const cube_stref_labeled_time &l, const cube_stref_labeled_time &r) {
        if (!(l._win.left == r._win.left &&
              l._win.right == r._win.right &&
              l._win.top == r._win.top &&
              l._win.bottom == r._win.bottom &&
              l._nx == r._nx &&
              l._ny == r._ny &&
              l._dt == r._dt &&
              l._t_values == r._t_values))
            return false;

        // compare SRS
        OGRSpatialReference a = l.srs_ogr();
        OGRSpatialReference b = r.srs_ogr();

        if (!a.IsSame(&b))
            return false;

        return true;
    }

    inline friend bool operator!=(const cube_stref_labeled_time &l, const cube_stref_labeled_time &r) { return !(l == r); }

    virtual std::shared_ptr<cube_stref> copy() override {
        std::shared_ptr<cube_stref_labeled_time> x = std::make_shared<cube_stref_labeled_time>();
        x->_win.left = _win.left;
        x->_win.right = _win.right;
        x->_win.top = _win.top;
        x->_win.bottom = _win.bottom;
        x->_nx = _nx;
        x->_ny = _ny;
        x->_dt = _dt;
        x->_t0 = _t0;
        x->_t1 = _t1;
        x->_srs = _srs;

        x->_t_values = _t_values;
        x->_t_index = _t_index;
        return x;
    }

   protected:
    std::vector<datetime> _t_values;
    std::map<datetime, uint32_t> _t_index;
};

}  // namespace gdalcubes

#endif  //VIEW_H
