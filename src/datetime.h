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
#include "error.h"
#include "external/date.h"

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

    std::string to_string() {
        switch (dt_unit) {
            case datetime_unit::NONE:
                return "";
            case datetime_unit::SECOND:
                return "PT" + std::to_string(dt_interval) + "S";
            case datetime_unit::MINUTE:
                return "PT" + std::to_string(dt_interval) + "M";
            case datetime_unit::HOUR:
                return "PT" + std::to_string(dt_interval) + "H";
            case datetime_unit::DAY:
                return "P" + std::to_string(dt_interval) + "D";
            case datetime_unit::WEEK:
                return "P" + std::to_string(dt_interval) + "W";
            case datetime_unit::MONTH:
                return "P" + std::to_string(dt_interval) + "M";
            case datetime_unit::YEAR:
                return "P" + std::to_string(dt_interval) + "Y";
        }
        return "";
    }

    static duration from_string(std::string s) {
        std::regex rexp("P(T?)([0-9]+)([YMWDHS])");

        std::cmatch res;
        if (!std::regex_match(s.c_str(), res, rexp)) {
            throw std::string("ERROR in duration::from_string(): cannot derive date interval");
        }
        duration d;
        d.dt_interval = std::stoi(res[2]);
        if (!res[1].str().empty()) {
            if (res[3] == "H")
                d.dt_unit = datetime_unit::HOUR;
            else if (res[3] == "M")
                d.dt_unit = datetime_unit::MINUTE;
            else if (res[3] == "S")
                d.dt_unit = datetime_unit::SECOND;
        } else {
            if (res[3] == "Y")
                d.dt_unit = datetime_unit::YEAR;
            else if (res[3] == "M")
                d.dt_unit = datetime_unit::MONTH;
            else if (res[3] == "W")
                d.dt_unit = datetime_unit::WEEK;
            else if (res[3] == "D")
                d.dt_unit = datetime_unit::DAY;
        }
        return d;
    }

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

    duration convert(datetime_unit u) {
        duration out;
        if (u == datetime_unit::NONE || dt_unit == datetime_unit::NONE) {
            GCBS_ERROR("Failed conversion of datetime duration with undefined unit");
            return out;
        }
        out.dt_unit = dt_unit;
        out.dt_interval = dt_interval;
        if (u == dt_unit) {
            return out;
        }
        while (out.dt_unit != u) {
            if (out.dt_unit < u) {
                switch (out.dt_unit) {
                    case datetime_unit::NONE:
                    case datetime_unit::YEAR:
                        GCBS_ERROR("Failed conversion of datetime duration with undefined unit");
                        out.dt_unit = dt_unit;
                        out.dt_interval = dt_interval;
                        return out;
                    case datetime_unit::SECOND:
                        out.dt_unit = datetime_unit::MINUTE;
                        out.dt_interval = (int)std::ceil((double)out.dt_interval / 60.0);
                        break;
                    case datetime_unit::MINUTE:
                        out.dt_unit = datetime_unit::HOUR;
                        out.dt_interval = (int)std::ceil((double)out.dt_interval / 60.0);
                        break;
                    case datetime_unit::HOUR:
                        out.dt_unit = datetime_unit::DAY;
                        out.dt_interval = (int)std::ceil((double)out.dt_interval / 24.0);
                        break;
                    case datetime_unit::DAY:
                        out.dt_unit = datetime_unit::MONTH;
                        out.dt_interval = (int)std::ceil((double)out.dt_interval / 30.0);
                        break;
                    case datetime_unit::WEEK:
                        out.dt_unit = datetime_unit::MONTH;
                        out.dt_interval = (int)std::ceil((double)out.dt_interval * 7 / 30.0);
                        break;
                    case datetime_unit::MONTH:
                        out.dt_unit = datetime_unit::YEAR;
                        out.dt_interval = (int)std::ceil((double)out.dt_interval / 12.0);
                        break;
                }
            } else {
                switch (out.dt_unit) {
                    case datetime_unit::NONE:
                    case datetime_unit::SECOND:
                        GCBS_ERROR("Failed conversion of datetime duration with undefined unit");
                        out.dt_unit = dt_unit;
                        out.dt_interval = dt_interval;
                        return out;
                    case datetime_unit::MINUTE:
                        out.dt_unit = datetime_unit::SECOND;
                        out.dt_interval = (int)std::ceil((double)out.dt_interval * 60.0);
                        break;
                    case datetime_unit::HOUR:
                        out.dt_unit = datetime_unit::MINUTE;
                        out.dt_interval = (int)std::ceil((double)out.dt_interval * 60.0);
                        break;
                    case datetime_unit::DAY:
                        out.dt_unit = datetime_unit::HOUR;
                        out.dt_interval = (int)std::ceil((double)out.dt_interval * 24.0);
                        break;
                    case datetime_unit::WEEK:
                        out.dt_unit = datetime_unit::DAY;
                        out.dt_interval = (int)std::ceil((double)out.dt_interval * 7.0);
                        break;
                    case datetime_unit::MONTH:
                        out.dt_unit = datetime_unit::DAY;
                        out.dt_interval = (int)std::ceil((double)out.dt_interval * 30.0);
                        break;
                    case datetime_unit::YEAR:
                        out.dt_unit = datetime_unit::MONTH;
                        out.dt_interval = (int)std::ceil((double)out.dt_interval * 12.0);
                        break;
                }
            }
        }
        return out;
    }
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
    double to_double() {
        auto daypoint = date::floor<date::days>(_p);
        auto ymd = date::year_month_day(daypoint);
        auto tod = date::make_time(_p - daypoint);

        double out = int(ymd.year());
        if (_unit <= datetime_unit::MONTH) {
            out *= 100;
            out += unsigned(ymd.month());
        }
        if (_unit <= datetime_unit::DAY) {
            out *= 100;
            out += unsigned(ymd.day());
        }
        if (_unit <= datetime_unit::HOUR) {
            out *= 100;
            out += tod.hours().count();
        }
        if (_unit <= datetime_unit::MINUTE) {
            out *= 100;
            out += tod.minutes().count();
        }
        if (_unit <= datetime_unit::SECOND) {
            out *= 100;
            out += tod.seconds().count();
        }
        return out;
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
    static date::sys_seconds tryparse(std::string format, std::string d) {
        bool success = false;
        date::sys_seconds out;  // TODO: set to invalid?!
        if (!success) {
            std::istringstream is(d);
            is >> date::parse(format, out);
            if (bool(is))
                success = true;
        }

#if !defined __GNUC__ || __GNUC__ >= 5 // gcc 4.9x misses std::get_time
        if (!success) {
            std::tm tp;
            tp.tm_sec = 0;
            tp.tm_min = 0;
            tp.tm_hour = 0;
            tp.tm_mday = 1;
            tp.tm_mon = 0;
            tp.tm_year = 0;
            tp.tm_wday = -1;
            tp.tm_yday = -1;

            std::istringstream is(d);
            is >> std::get_time(&tp, format.c_str());  // works only from GCC > 5
            if (!is.fail()) {
                if (tp.tm_yday != -1) {
                    out = date::sys_days{date::year{tp.tm_year + 1900} / 1 / 1} + date::days{tp.tm_yday} +
                          std::chrono::hours{tp.tm_hour} + std::chrono::minutes{tp.tm_min} + std::chrono::seconds{tp.tm_sec};
                } else {
                    out = date::sys_days{date::year{tp.tm_year + 1900} / (tp.tm_mon + 1) / tp.tm_mday} +
                          std::chrono::hours{tp.tm_hour} + std::chrono::minutes{tp.tm_min} + std::chrono::seconds{tp.tm_sec};
                }
                success = true;
            }
        }
#endif

        if (!success) {
            GCBS_ERROR("Cannot parse datetime string '" + d + "' with format '" + format + "'");
            throw std::string("Cannot parse datetime string '" + d + "' with format '" + format + "'");
        }
        return out;
    }

    // from standard format with variable precision (
    static datetime from_string(std::string s) {
        std::istringstream is(s);

        // TODO: Regex does not support ISO weeks / day of year yet
        std::regex regex1("([0-9]{4})(?:-?([0-9]{2})(?:-?([0-9]{2})(?:[T\\s]?([0-9]{2})(?::?([0-9]{2})(?::?([0-9]{2}))?)?)?)?)?");
        datetime out;

        std::cmatch res;
        if (!std::regex_match(s.c_str(), res, regex1)) {
            throw std::string("ERROR in datetime::from_string(): cannot derive datetime from string");
        } else {
            if (res.size() != 7) throw std::string("ERROR in datetime::from_string(): cannot derive datetime from string");
            uint16_t i = 2;
            while (!res[i].str().empty() && i < 7) ++i;
            if (i == 2) {
                //date::sys_days d = date::year(std::stoi(res[1].str()));
                out._p = date::sys_days{date::year(std::stoi(res[1].str())) / date::month(1) / date::day(1)};
                out._unit = datetime_unit::YEAR;
            } else if (i == 3) {
                out._p = date::sys_days{date::year(std::stoi(res[1].str())) / date::month(std::stoi(res[2].str())) / date::day(1)};
                out._unit = datetime_unit::MONTH;
            } else if (i == 4) {
                out._p = date::sys_days{date::year(std::stoi(res[1].str())) / date::month(std::stoi(res[2].str())) / date::day(std::stoi(res[3].str()))};
                out._unit = datetime_unit::DAY;
            } else if (i == 5) {
                out._p = date::sys_days{date::year(std::stoi(res[1].str())) / date::month(std::stoi(res[2].str())) / date::day(std::stoi(res[3].str()))} +
                         std::chrono::hours{std::stoi(res[4].str())};
                out._unit = datetime_unit::HOUR;
            } else if (i == 6) {
                out._p = date::sys_days{date::year(std::stoi(res[1].str())) / date::month(std::stoi(res[2].str())) / date::day(std::stoi(res[3].str()))} +
                         std::chrono::hours{std::stoi(res[4].str())} + std::chrono::minutes{std::stoi(res[5].str())};
                out._unit = datetime_unit::MINUTE;
            } else if (i == 7) {
                out._p = date::sys_days{date::year(std::stoi(res[1].str())) / date::month(std::stoi(res[2].str())) / date::day(std::stoi(res[3].str()))} +
                         std::chrono::hours{std::stoi(res[4].str())} + std::chrono::minutes{std::stoi(res[5].str())} + std::chrono::seconds{std::stoi(res[6].str())};
                out._unit = datetime_unit::SECOND;
            }
        }
        return out;
    }

    datetime_unit& unit() { return _unit; }
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
        //        switch (l.unit()) {
        //            case datetime_unit::NONE:
        //                return false;
        //            case datetime_unit::SECOND:
        //                return ((l - r).dt_interval == 0);
        //            case datetime_unit::MINUTE:
        //                return ((l - r).dt_interval == 0);
        //            case datetime_unit::HOUR:
        //                return ((l - r).dt_interval == 0);
        //            case datetime_unit::DAY:
        //                return ((l - r).dt_interval == 0);
        //            case datetime_unit::WEEK: {
        //                return ((l - r).dt_interval == 0);
        //            }
        //            case datetime_unit::MONTH:
        //                return ((l - r).dt_interval == 0);
        //            case datetime_unit::YEAR:
        //                return ((l - r).dt_interval == 0);
        //        }
        //        return false;
    }

    inline friend bool operator!=(const datetime& l, const datetime& r) { return !(l == r); }

    friend bool operator<(datetime& l, datetime& r) {
        if (l.unit() != r._unit) return false;

        if (l.unit() == datetime_unit::NONE) return false;
        return ((l - r).dt_interval < 0);
    }
    inline friend bool operator>(datetime& l, datetime& r) { return r < l; }
    inline friend bool operator<=(datetime& l, datetime& r) { return !(l > r); }
    inline friend bool operator>=(datetime& l, datetime& r) { return !(l < r); }

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

    static std::string datetime_format_for_unit(datetime_unit u) {
        switch (u) {
            case datetime_unit::NONE:
            case datetime_unit::SECOND:
                return "%Y-%m-%dT%H:%M:%S";
            case datetime_unit::MINUTE:
                return "%Y-%m-%dT%H:%M";
            case datetime_unit::HOUR:
                return "%Y-%m-%dT%H";
            case datetime_unit::DAY:
                return "%Y-%m-%d";
            case datetime_unit::WEEK:
                return "%Y-%m-%dT%H:%M:%S";
            case datetime_unit::MONTH:
                return "%Y-%m";
            case datetime_unit::YEAR:
                return "%Y";
        }
        return "%Y-%m-%dT%H:%M:%S";
    }
};

#endif  //DATETIME_H
