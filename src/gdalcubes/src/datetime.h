/*
    MIT License

    Copyright (c) 2021 Marius Appel <marius.appel@uni-muenster.de>

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
 * Simplistic library for datetime objects.
 * The assumption of the library is that durations are only meaningful if given as a single number with a single datetime unit. Whether or not this assumption
 * is valid depends on the application.
 */

#ifndef DATETIME_H
#define DATETIME_H

// gcc 4.9x misses std::time_get and std::time_put, which is used in the date library
#if defined __GNUC__ && __GNUC__ < 5
#define ONLY_C_LOCALE 1
#endif

#include <chrono>
#include <iomanip>
#include <regex>
#include <cstdint> // 2023-01-12: GCC 13 compatibility
#include "error.h"
#include "external/date.h"

namespace gdalcubes {

enum class datetime_unit {
    SECOND = 0,
    MINUTE = 1,
    HOUR = 2,
    DAY = 3,
    WEEK = 4,
    MONTH = 5,
    YEAR = 6,
    NONE = 255
};

struct duration {
    duration() : dt_interval(0), dt_unit(datetime_unit::DAY) {}
    duration(int32_t interval, datetime_unit unit) : dt_interval(interval), dt_unit(unit) {}

    int32_t dt_interval;
    datetime_unit dt_unit;

    std::string to_string();

    static duration from_string(std::string s);

    friend duration operator*(duration l, const int& r) {
        duration out;
        out.dt_unit = l.dt_unit;
        out.dt_interval = l.dt_interval * r;
        return out;
    }

    friend duration operator+(duration l, const int& r) {
        duration out;
        out.dt_unit = l.dt_unit;
        out.dt_interval = l.dt_interval + r;
        return out;
    }

    friend duration operator/(duration l, const int& r) {
        duration out;
        out.dt_unit = l.dt_unit;
        out.dt_interval = l.dt_interval / r;
        return out;
    }

    friend int operator/(duration l, duration& r) {
        if (l.dt_unit != r.dt_unit) {
            throw std::string("ERROR in duration::operator/(): Incompatible datetime duration units");
        }
        if (r.dt_interval == 0) {
            throw std::string("ERROR in duration::operator/(): Division by zero");
        }

        return l.dt_interval / r.dt_interval;
    }

    friend int operator%(duration l, duration& r) {
        if (l.dt_unit != r.dt_unit) {
            throw std::string("ERROR in duration::operator%(): Incompatible datetime duration units");
        }
        if (r.dt_interval == 0) {
            throw std::string("ERROR in duration::operator%(): Division by zero");
        }

        return l.dt_interval % r.dt_interval;
    }

    friend duration operator-(duration l, const int& r) {
        duration out;
        out.dt_unit = l.dt_unit;
        out.dt_interval = l.dt_interval - r;
        return out;
    }

    friend bool operator==(const duration& l, const duration& r) {
        return (l.dt_unit == r.dt_unit && l.dt_interval == r.dt_interval);
    }

    inline friend bool operator!=(const duration& l, const duration& r) { return !(l == r); }

    duration convert(datetime_unit u);
};

/**
 * @brief Simplistic datetime class
 *
 * This class represents date and time with regard to a unit / granularity.
 * Internally it stores a posixtime object and a numeric month and year.
 * The unit / granularity is derived automatically from given strings.
 *
 */
class datetime {
   public:
    datetime() : _p(), _unit(datetime_unit::DAY) {}
    datetime(date::sys_seconds p) : _p(p), _unit(datetime_unit::DAY) {}
    datetime(date::sys_seconds p, datetime_unit u) : _p(p), _unit(u) {}

    /**
     * Convert to a numeric datetime format as in 20180401123059.
     * This function does **not** convert the datetime to a timestamp or similar
     * @return
     */
    double to_double();

    long seconds() {
        auto daypoint = date::floor<date::days>(_p);
        auto tod = date::make_time(_p - daypoint);
        return tod.seconds().count();
    }

    long minutes() {
        auto daypoint = date::floor<date::days>(_p);
        auto tod = date::make_time(_p - daypoint);
        return tod.minutes().count();
    }

    long hours() {
        auto daypoint = date::floor<date::days>(_p);
        auto tod = date::make_time(_p - daypoint);
        return tod.hours().count();
    }

    int year() {
        auto daypoint = date::floor<date::days>(_p);
        auto ymd = date::year_month_day(daypoint);
        return int(ymd.year());
    }

    unsigned int month() {
        auto daypoint = date::floor<date::days>(_p);
        auto ymd = date::year_month_day(daypoint);
        return unsigned(ymd.month());
    }

    unsigned int dayofmonth() {
        auto daypoint = date::floor<date::days>(_p);
        auto ymd = date::year_month_day(daypoint);
        return unsigned(ymd.day());
    }

    int dayofyear() {
        auto daypoint = date::floor<date::days>(_p);
        auto ymd = date::year_month_day(daypoint);
        auto ymd1 = date::year_month_day(ymd.year(), date::month(1), date::day(1));
        return (date::sys_days(ymd) - date::sys_days(ymd1)).count() + 1;
    }

    int dayofweek() {
        auto daypoint = date::floor<date::days>(_p);
        auto ymd = date::year_month_day(daypoint);
        return (date::weekday(date::sys_days(ymd)) - date::Sunday).count();  // days since sunday
    }

    double epoch_time() {
        return (double)(_p.time_since_epoch().count());
    }

    /**
     * @brief Convert to a simple datetime string
     * @return
     */
    std::string to_string() {
        std::stringstream os;
        os << date::format(datetime_format_for_unit(_unit).c_str(), _p);
        return os.str();
    }

    /**
     * Produce a datetime string up to the given unit.
     * The string will possibly include e.g. seconds although the unit is actually lower e.g. days.
     * This function is implemented to interface tools that require full datetime strings and do not work with
     * e.g. "2001-01".
     * @param u datetime unit
     * @return
     */
    std::string to_string(datetime_unit u) {
        std::stringstream os;
        os << date::format(datetime_format_for_unit(u).c_str(), _p);
        return os.str();
    }

    // Helper function that tries to parse string datetimes according to a given format
    // tries date::parse first and if this does not work std::get_time
    static date::sys_seconds tryparse(std::string format, std::string d);

    // from standard format with variable precision
    static datetime from_string(std::string s);


    static std::string unit_to_string(datetime_unit u) {
        switch(u) {
            case datetime_unit::NONE:
                return "none";
            case datetime_unit::YEAR:
                return "years";
            case datetime_unit::MONTH:
                return "months";
            case datetime_unit::WEEK:
                return "weeks";
            case datetime_unit::DAY:
                return "days";
            case datetime_unit::HOUR:
                return "hours";
            case datetime_unit::MINUTE:
                return "minutes";
            case datetime_unit::SECOND:
                return "seconds";
        }
        return "none";
    }

    /**
     * Optimistic string parser assuming
     * year is given by the first four digits and other (optional)
     * fields are given by two digits in the order YmdHMS. Any kind of separators as well
     * as time zones will simply be skipped.
     * @param s string
     * @return
     */
    static datetime from_YmdHMS_digits(std::string s) ;

    void unit(datetime_unit u) {

        // Reset finer datetime components
        switch (u) {
            case datetime_unit::NONE:
                break;
            case datetime_unit::SECOND:
                break;
            case datetime_unit::MINUTE:
                _p = date::sys_days{date::year(this->year()) / date::month(this->month()) / date::day(this->dayofmonth())} +
                    std::chrono::hours{this->hours()} + std::chrono::minutes{this->minutes()} + std::chrono::seconds{0};
                break;
            case datetime_unit::HOUR:
                _p = date::sys_days{date::year(this->year()) / date::month(this->month()) / date::day(this->dayofmonth())} +
                     std::chrono::hours{this->hours()} + std::chrono::minutes{0} + std::chrono::seconds{0};
                break;
            case datetime_unit::DAY:
                _p = date::sys_days{date::year(this->year()) / date::month(this->month()) / date::day(this->dayofmonth())} +
                     std::chrono::hours{0} + std::chrono::minutes{0} + std::chrono::seconds{0};
                break;
            case datetime_unit::WEEK:
                _p = date::sys_days{date::year(this->year()) / date::month(this->month()) / date::day(this->dayofmonth())} +
                     std::chrono::hours{0} + std::chrono::minutes{0} + std::chrono::seconds{0};
                break;
            case datetime_unit::MONTH:
                _p = date::sys_days{date::year(this->year()) / date::month(this->month()) / date::day(1)} +
                     std::chrono::hours{0} + std::chrono::minutes{0} + std::chrono::seconds{0};
                break;
            case datetime_unit::YEAR:
               _p = date::sys_days{date::year(this->year()) / date::month(1) / date::day(1)} +
                    std::chrono::hours{0} + std::chrono::minutes{0} + std::chrono::seconds{0};
                break;
        }
        _unit = u;

    }
    datetime_unit unit() const { return _unit; }

    friend duration operator-(const datetime& l,
                              const datetime& r) {
        duration out;
        out.dt_unit = std::max(l.unit(), r.unit());  // TODO: warning if not the same unit
        //GCBS_WARN("datetime objects have different units");
        // std::tm ltm = l._p;
        //std::tm rtm = r._p;

        auto dp_l = date::floor<date::days>(l._p);
        auto ymd_l = date::year_month_day(dp_l);
        auto dp_r = date::floor<date::days>(r._p);
        auto ymd_r = date::year_month_day(dp_r);

        switch (out.dt_unit) {
            case datetime_unit::NONE:
                break;
            case datetime_unit::SECOND:
                out.dt_interval = std::chrono::duration_cast<std::chrono::seconds>(l._p - r._p).count();
                break;
            case datetime_unit::MINUTE:
                out.dt_interval = std::chrono::duration_cast<std::chrono::minutes>(l._p - r._p).count();
                break;
            case datetime_unit::HOUR:
                out.dt_interval = std::chrono::duration_cast<std::chrono::hours>(l._p - r._p).count();
                break;
            case datetime_unit::DAY:
                out.dt_interval = (dp_l - dp_r).count();
                break;
            case datetime_unit::WEEK:
                out.dt_interval = (dp_l - dp_r).count() / 7;
                break;
            case datetime_unit::MONTH: {
                // see https://github.com/HowardHinnant/date/wiki/Examples-and-Recipes#deltamonths
                out.dt_interval = (ymd_l.year() / ymd_l.month() - ymd_r.year() / ymd_r.month()).count();
                break;
            }
            case datetime_unit::YEAR:
                out.dt_interval = (ymd_l.year() - ymd_r.year()).count();
                break;
        }
        return out;
    }

    friend bool operator==(const datetime& l, const datetime& r) {
        if (l.unit() != r._unit) return false;  // TODO: warning
        return ((l - r).dt_interval == 0);
    }

    inline friend bool operator!=(const datetime& l, const datetime& r) { return !(l == r); }

    friend bool operator<(const datetime& l, const datetime& r) {
        if (l.unit() != r._unit) return false;

        if (l.unit() == datetime_unit::NONE) return false;
        return ((l - r).dt_interval < 0);
    }
    inline friend bool operator>(const datetime& l, const datetime& r) { return r < l; }
    inline friend bool operator<=(const datetime& l, const datetime& r) { return !(l > r); }
    inline friend bool operator>=(const datetime& l, const datetime& r) { return !(l < r); }

    friend datetime operator+(datetime l, const duration& r) {
        datetime out(l._p);
        auto dp = date::floor<date::days>(l._p);
        auto ymd = date::year_month_day(dp);
        auto tod = l._p - dp;
        out._unit = r.dt_unit;
        switch (r.dt_unit) {
            case datetime_unit::NONE:
                break;
            case datetime_unit::SECOND:
                out._p += std::chrono::seconds{r.dt_interval};
                break;
            case datetime_unit::MINUTE:
                out._p += std::chrono::minutes{r.dt_interval};
                break;
            case datetime_unit::HOUR:
                out._p += std::chrono::hours{r.dt_interval};
                break;
            case datetime_unit::DAY:
                out._p = dp + date::days{r.dt_interval} + tod;
                break;
            case datetime_unit::WEEK:
                out._p = dp + date::days{7 * r.dt_interval} + tod;
                break;
            case datetime_unit::MONTH:
                ymd = ymd + date::months{r.dt_interval};
                if (!ymd.ok()) {
                    ymd = ymd.year() / ymd.month() / date::last;  // use last day of month if needed
                }
                out._p = date::sys_days{ymd} + tod;
                break;
            case datetime_unit::YEAR:
                ymd = ymd + date::years{r.dt_interval};
                if (!ymd.ok()) {
                    ymd = ymd.year() / ymd.month() / date::last;  // use last day of month if needed (e.g. 2016-02-29 + 1 YEAR -> 2017-02-28)
                }
                out._p = date::sys_days{ymd} + tod;
                break;
        }
        return out;
    }

   private:
    date::sys_seconds _p;
    datetime_unit _unit;

    static std::string datetime_format_for_unit(datetime_unit u);
};

}  // namespace gdalcubes

#endif  //DATETIME_H
