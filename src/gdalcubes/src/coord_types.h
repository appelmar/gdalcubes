/*
    MIT License

    Copyright (c) 2020 Marius Appel <marius.appel@hs-bochum.de>

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

#ifndef COORD_TYPES_H
#define COORD_TYPES_H

//#include <ogr_spatialref.h>

#include <array>
#include <limits>
#include <cstdint>
#include "datetime.h"

namespace gdalcubes {

template <typename Ta>
struct bounds_2d {
    Ta left, bottom, top, right;
    static bool intersects(bounds_2d<Ta> a, bounds_2d<Ta> b) {
        return (
            a.right >= b.left &&
            a.left <= b.right &&
            a.top >= b.bottom &&
            a.bottom <= b.top);
    }
    static bool within(bounds_2d<Ta> a, bounds_2d<Ta> b) {
        return (
            a.left >= b.left &&
            a.right <= b.right &&
            a.top <= b.top &&
            a.bottom >= b.bottom);
    }
    static bool outside(bounds_2d<Ta> a, bounds_2d<Ta> b) {
        return !intersects(a, b);
    }
    static bounds_2d<Ta> union2d(bounds_2d<Ta> a, bounds_2d<Ta> b) {
        bounds_2d<Ta> out;
        out.left = fmin(a.left, b.left);
        out.right = fmax(a.right, b.right);
        out.bottom = fmin(a.bottom, b.bottom);
        out.top = fmax(a.top, b.top);
        return out;
    }

    static bounds_2d<Ta> intersection(bounds_2d<Ta> a, bounds_2d<Ta> b) {
        bounds_2d<Ta> out;
        out.left = fmax(a.left, b.left);
        out.right = fmin(a.right, b.right);
        out.bottom = fmin(a.bottom, b.bottom);
        out.top = fmax(a.top, b.top);
        return out;
    }

    bounds_2d<Ta> transform(std::string srs_from, std::string srs_to) {
        if (srs_from == srs_to) {
            return *this;
        }
        OGRSpatialReference srs_in;
        OGRSpatialReference srs_out;
        srs_in.SetFromUserInput(srs_from.c_str());
        srs_out.SetFromUserInput(srs_to.c_str());

        if (srs_in.IsSame(&srs_out)) {
            return *this;
        }

        OGRCoordinateTransformation* coord_transform = OGRCreateCoordinateTransformation(&srs_in, &srs_out);

        Ta x[4] = {left, left, right, right};
        Ta y[4] = {top, bottom, top, bottom};

        if (coord_transform == NULL || !coord_transform->Transform(4, x, y)) {
            throw std::string("ERROR: coordinate transformation failed (from " + srs_from + " to " + srs_to + ").");
        }

        Ta xmin = std::numeric_limits<Ta>::is_integer ? std::numeric_limits<Ta>::max() : std::numeric_limits<Ta>::max();
        Ta ymin = std::numeric_limits<Ta>::is_integer ? std::numeric_limits<Ta>::max() : std::numeric_limits<Ta>::max();
        Ta xmax = std::numeric_limits<Ta>::is_integer ? std::numeric_limits<Ta>::min() : -std::numeric_limits<Ta>::max();
        Ta ymax = std::numeric_limits<Ta>::is_integer ? std::numeric_limits<Ta>::min() : -std::numeric_limits<Ta>::max();
        for (uint8_t k = 0; k < 4; ++k) {
            if (x[k] < xmin) xmin = x[k];
            if (y[k] < ymin) ymin = y[k];
            if (x[k] > xmax) xmax = x[k];
            if (y[k] > ymax) ymax = y[k];
        }

        left = xmin;
        right = xmax;
        top = ymax;
        bottom = ymin;

        OCTDestroyCoordinateTransformation(coord_transform);

        return *this;
    }
};

template <typename T>
struct coords_2d {
    T x, y;
};

//template <typename T> using coords_nd = std::vector<T>;
template <typename T, uint16_t N>
using coords_nd = std::array<T, N>;

template <typename T, uint16_t N>
struct bounds_nd {
    coords_nd<T, N> low;
    coords_nd<T, N> high;
};

struct coords_st {
    coords_2d<double> s;
    datetime t;
};

struct bounds_st {
    bounds_2d<double> s;
    datetime t0;
    datetime t1;
};

}  // namespace gdalcubes

#endif  // COORD_TYPES_H
