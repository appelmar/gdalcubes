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



#include <fstream>
#include <algorithm>
#include "filesystem.h"
#include "view.h"

namespace gdalcubes {



aggregation::aggregation_type aggregation::from_string(std::string s) {
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


std::string aggregation::to_string(aggregation::aggregation_type a) {
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


resampling::resampling_type resampling::from_string(std::string s) {
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

std::string resampling::to_string(resampling::resampling_type r) {
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

GDALRIOResampleAlg resampling::to_gdal_rasterio(resampling::resampling_type r) {
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



std::string cube_stref::type_string(std::shared_ptr<cube_stref> obj) {
    if (std::dynamic_pointer_cast<cube_stref_labeled_time>(obj) != nullptr) {
        return "cube_stref_labeled_time";
    }
    if (std::dynamic_pointer_cast<cube_stref_regular>(obj) != nullptr) {
        return "cube_stref_regular";
    }
    return "";
}


void cube_stref_regular::set_x_axis(double min, double max, double delta) {
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


void cube_stref_regular::set_x_axis(double min, double max, uint32_t n) {
    _win.left = min;
    _win.right = max;
    _nx = n;
}
void cube_stref_regular::set_x_axis(double min, uint32_t n, double delta) {
    _win.left = min;
    _nx = n;
    _win.right = min + n * delta;
}
void cube_stref_regular::set_x_axis(uint32_t n, double max, double delta) {
    _win.right = max;
    _nx = n;
    _win.left = max - n * delta;
}


void cube_stref_regular::set_y_axis(double min, double max, double delta) {
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
void cube_stref_regular::set_y_axis(double min, double max, uint32_t n) {
    _win.bottom = min;
    _win.top = max;
    _ny = n;
}
void cube_stref_regular::set_y_axis(double min, uint32_t n, double delta) {
    _win.bottom = min;
    _ny = n;
    _win.top = min + n * delta;

}
void cube_stref_regular::set_y_axis(uint32_t n, double max, double delta) {
    _win.top = max;
    _ny = n;
    _win.bottom = max - n * delta;
}




void cube_stref_regular::set_t_axis(datetime min, datetime max, duration delta) {
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


void cube_stref_regular::set_t_axis(datetime min, datetime max, uint32_t n) {
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

OGRSpatialReference cube_stref_regular::srs_ogr() const  {
    OGRSpatialReference s;
    s.SetFromUserInput(_srs.c_str());
    return s;
}

uint32_t cube_stref_regular::nt() {
    if (_t1 == _t0) return 1;
    duration d = (_t1 - _t0) + 1;
    return (d % _dt == 0) ? d / _dt : (1 + (d / _dt));
}


coords_st cube_stref_regular::map_coords(coords_nd<uint32_t, 3> p) {
    coords_st s;
    s.s.x = _win.left + p[2] * dx();
    s.s.y = _win.top - p[1] * dy();
    s.t = _t0 + _dt * p[0];
    return s;
}


coords_nd<uint32_t, 3> cube_stref_regular::cube_coords(coords_st p) {
    coords_nd<uint32_t, 3> s;
    s[2] = (uint32_t)((p.s.x - _win.left) / dx());
    s[1] = (uint32_t)((_win.top - p.s.y) / dy());
    s[0] = (uint32_t)((p.t - _t0) / _dt);
    return s;
}


std::shared_ptr<cube_stref> cube_stref_regular::copy() {
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



bool operator==(const cube_stref_regular &l, const cube_stref_regular &r) {
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



void cube_stref_labeled_time::set_time_labels(std::vector<datetime> t) {
    // TODO: if t.empty()
    _t_values = t;

    // TODO: what if not sorted?
    for (uint32_t i = 0; i < _t_values.size(); ++i) {
        _t_index.insert(std::make_pair(t[i], i));
    }
}

 void cube_stref_labeled_time::set_time_labels(std::vector<std::string> t) {
    // TODO: what if not sorted?
    for (uint32_t i = 0; i < _t_values.size(); ++i) {
        _t_index.insert(std::make_pair(datetime::from_string(t[i]), i));
    }
}


std::vector<std::string> cube_stref_labeled_time::get_time_labels_as_string() {
    std::vector<std::string> out;
    for (uint32_t i = 0; i < _t_values.size(); ++i) {
        out.push_back(_t_values[i].to_string());
    }
    return out;
}

coords_st  cube_stref_labeled_time::map_coords(coords_nd<uint32_t, 3> p) {
    coords_st s;
    s.s.x = _win.left + p[2] * dx();
    s.s.y = _win.top - p[1] * dy();
    s.t = datetime_at_index(p[0]);
    return s;
}

coords_nd<uint32_t, 3> cube_stref_labeled_time::cube_coords(coords_st p) {
    coords_nd<uint32_t, 3> s;
    s[2] = (uint32_t)((p.s.x - _win.left) / dx());
    s[1] = (uint32_t)((_win.top - p.s.y) / dy());
    s[0] = index_at_datetime(p.t);
    return s;
}

uint32_t cube_stref_labeled_time::index_at_datetime(datetime t) {
    auto res = _t_index.find(t);
    if (res == _t_index.end()) {
        GCBS_ERROR("Data cubes does not contain time slice for requested datetime");
        throw std::string("Data cubes does not contain time slice for requested datetime");
    }
    return res->second;
}



bool operator==(const cube_stref_labeled_time &l, const cube_stref_labeled_time &r) {
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

std::shared_ptr<cube_stref> cube_stref_labeled_time::copy()  {
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



cube_view cube_view::read(json11::Json j) {
    cube_view v;

    v._dt = duration::from_string(j["time"]["dt"].string_value());

    std::string st0 = j["time"]["t0"].string_value();
    std::string st1 = j["time"]["t1"].string_value();

    v._t0 = datetime::from_string(st0);
    v._t1 = datetime::from_string(st1);

    // No matter how the start and end datetime are given, use the unit of the datetime interval!
    v._t0.unit(v._dt.dt_unit);
    v._t1.unit(v._dt.dt_unit);

    if (!j["space"].is_null()) {
        v._win.left = j["space"]["left"].number_value();
        v._win.right = j["space"]["right"].number_value();
        v._win.top = j["space"]["top"].number_value();
        v._win.bottom = j["space"]["bottom"].number_value();

        if (!j["space"]["nx"].is_null() && !j["space"]["ny"].is_null()) {
            v._nx = j["space"]["nx"].int_value();
            v._ny = j["space"]["ny"].int_value();
        }
        // if only one of nx and ny is given, use the aspect ratio of the spatial extent to derive the other
        else if (!j["space"]["nx"].is_null() && j["space"]["ny"].is_null()) {
            v._nx = j["space"]["nx"].int_value();
            v._ny = (uint32_t)((double)v._nx * (v._win.top - v._win.bottom) / (v._win.right - v._win.left));
        } else if (j["space"]["nx"].is_null() && !j["space"]["ny"].is_null()) {
            v._ny = j["space"]["ny"].int_value();
            v._nx = (uint32_t)((double)v._ny * (v._win.right - v._win.left) / (v._win.top - v._win.bottom));
        } else {
            throw std::string("ERROR in cube_view::read_json_string(): at least one of nx or ny must be given.");
        }

        v._srs = j["space"]["srs"].string_value();

    } else if (!j["tile"].is_null()) {
        const double EARTH_RADIUS_METERS = 6378137;
        const uint16_t TILE_SIZE_PX = 256;
        const double EARTH_CIRCUMFERENCE_METERS = 2 * M_PI * EARTH_RADIUS_METERS;
        const double WEBMERCATOR_BOUNDS_LEFT = -EARTH_CIRCUMFERENCE_METERS / 2.0;  // -20037508.342789244
        //const double WEBMERCATOR_BOUNDS_LOWER = -EARTH_CIRCUMFERENCE_METERS / 2.0;  // -20037508.342789244
        //const double WEBMERCATOR_BOUNDS_RIGHT = EARTH_CIRCUMFERENCE_METERS / 2.0;   // 20037508.342789244
        const double WEBMERCATOR_BOUNDS_UPPER = EARTH_CIRCUMFERENCE_METERS / 2.0;  // 20037508.342789244

        uint32_t x = j["tile"]["x"].int_value();
        uint32_t y = j["tile"]["y"].int_value();
        uint32_t z = j["tile"]["z"].int_value();

        bounds_2d<double> win;
        win.left = WEBMERCATOR_BOUNDS_LEFT + x * EARTH_CIRCUMFERENCE_METERS / std::pow(2, z);
        win.right = WEBMERCATOR_BOUNDS_LEFT + (x + 1) * EARTH_CIRCUMFERENCE_METERS / std::pow(2, z);
        win.top = WEBMERCATOR_BOUNDS_UPPER - y * EARTH_CIRCUMFERENCE_METERS / std::pow(2, z);
        win.bottom = WEBMERCATOR_BOUNDS_UPPER - (y + 1) * EARTH_CIRCUMFERENCE_METERS / std::pow(2, z);

        v._win = win;
        v._nx = TILE_SIZE_PX;
        v._ny = TILE_SIZE_PX;

        v._srs = "EPSG:3857";
    } else {
        throw std::string("ERROR in cube_view::read(): expected either 'space' or 'tile' in JSON cube view");
    }

    if (j["resampling"].is_null()) {
        v._resampling = resampling::resampling_type::RSMPL_NEAR;
    } else {
        v._resampling = resampling::from_string(j["resampling"].string_value());
    }

    if (j["aggregation"].is_null()) {
        v._aggregation = aggregation::aggregation_type::AGG_NONE;
    } else {
        v._aggregation = aggregation::from_string(j["aggregation"].string_value());
    }

    return v;
}

cube_view cube_view::read_json(std::string filename) {
    if (!filesystem::exists(filename))
        throw std::string("ERROR in cube_view::read_json(): image_collection_cube view file does not exist.");
    std::ifstream i(filename);

    std::stringstream buf;
    buf << i.rdbuf();

    std::string err;  // TODO: do something with error
    json11::Json j = json11::Json::parse(buf.str(), err);
    return read(j);
}

cube_view cube_view::read_json_string(std::string str) {
    std::istringstream i(str);
    std::string err;  // TODO: do something with error
    json11::Json j = json11::Json::parse(str, err);
    return read(j);
}

void cube_view::write_json(std::string filename) {
    json11::Json j = json11::Json::object{
        {"space", json11::Json::object{{"nx", (int)_nx}, {"ny", (int)_ny}, {"left", _win.left}, {"right", _win.right}, {"top", _win.top}, {"bottom", _win.bottom}, {"srs", _srs}}},
        {"time", json11::Json::object{{"dt", dt().to_string()}, {"t0", _t0.to_string()}, {"t1", _t1.to_string()}}},
        {"aggregation", aggregation::to_string(_aggregation)},
        {"resampling", resampling::to_string(_resampling)}};

    std::ofstream o(filename, std::ofstream::out);
    if (!o.good()) {
        throw std::string("ERROR in cube_view::write_json(): cannot write to file.");
    }

    o << std::setw(4) << j.dump() << std::endl;
    o.close();
}

std::string cube_view::write_json_string() {
    json11::Json j = json11::Json::object{
        {"space", json11::Json::object{{"nx", (int)_nx}, {"ny", (int)_ny}, {"left", _win.left}, {"right", _win.right}, {"top", _win.top}, {"bottom", _win.bottom}, {"srs", _srs}}},
        {"time", json11::Json::object{{"dt", dt().to_string()}, {"t0", _t0.to_string()}, {"t1", _t1.to_string()}}},
        {"aggregation", aggregation::to_string(_aggregation)},
        {"resampling", resampling::to_string(_resampling)}};
    std::ostringstream o;
    return j.dump();
}

}  // namespace gdalcubes
