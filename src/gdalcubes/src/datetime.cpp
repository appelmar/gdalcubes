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

#include "datetime.h"
#include <cctype>
#include <sstream>

namespace gdalcubes {

std::string duration::to_string() {
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

duration duration::from_string(std::string s) {
    duration d; // result object
    std::size_t pos = s.find("P");
    if (pos == std::string::npos) {
        throw std::string("ERROR in duration::from_string(): cannot derive date interval");
    }
    ++pos;
    bool is_time = false;
    if (s[pos] == 'T')  {
        ++pos;
        is_time = true;
    }
    std::stringstream numstr;
    while (pos < s.length()) {
        if (std::isdigit(s[pos])) {
            numstr << s[pos];
        }
        else {
            break;
        }
        ++pos;
    }
    if (numstr.str().empty()) {
        throw std::string("ERROR in duration::from_string(): cannot derive date interval");
    }
    d.dt_interval = std::stoi(numstr.str());

    if (pos >= s.length()) {
        throw std::string("ERROR in duration::from_string(): cannot derive date interval");
    }
    if (is_time) {
        if (s[pos] == 'H') {
            d.dt_unit = datetime_unit::HOUR;
        }
        else if (s[pos] == 'M') {
            d.dt_unit = datetime_unit::MINUTE;
        }
        else if (s[pos] == 'S') {
            d.dt_unit = datetime_unit::SECOND;
        }
        else {
            throw std::string("ERROR in duration::from_string(): cannot derive date interval, no valid datetime unit");
        }
    }
    else {
        if (s[pos] == 'Y') {
            d.dt_unit = datetime_unit::YEAR;
        }
        else if (s[pos] == 'M') {
            d.dt_unit = datetime_unit::MONTH;
        }
        else if (s[pos] == 'W') {
            d.dt_unit = datetime_unit::WEEK;
        }
        else if (s[pos] == 'D') {
            d.dt_unit = datetime_unit::DAY;
        }
        else {
            throw std::string("ERROR in duration::from_string(): cannot derive date interval, no valid datetime unit");
        }
    }
    return d;
}

duration duration::convert(gdalcubes::datetime_unit u) {
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

double datetime::to_double() {
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

date::sys_seconds datetime::tryparse(std::string format, std::string d) {
    bool success = false;
    date::sys_seconds out;  // TODO: set to invalid?!
    if (!success) {
        std::istringstream is(d);
        is >> date::parse(format, out);
        if (bool(is))
            success = true;
    }

#if !defined __GNUC__ || __GNUC__ >= 5  // gcc 4.9x misses std::get_time
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




datetime datetime::from_YmdHMS_digits(std::string s) {

    std::vector<uint8_t> digits;
    for(uint32_t i=0; i < s.length(); ++i) {
        if (std::isdigit(s[i]))
            digits.push_back(s[i] - '0');
    }
    datetime out;

    //YYYYMMDD HHMMSS
    if (digits.size() >= 14) {
        uint32_t Y = digits[0] * 1000 + digits[1] * 100 + digits[2] * 10 + digits[3];
        uint32_t m = digits[4] * 10 + digits[5];
        uint32_t d = digits[6] * 10 + digits[7];
        uint32_t H = digits[8] * 10 + digits[9];
        uint32_t M = digits[10] * 10 + digits[11];
        uint32_t S = digits[12] * 10 + digits[13];
        out._p = date::sys_days{date::year(Y) / date::month(m) / date::day(d)} +
                 std::chrono::hours{H} + std::chrono::minutes{M} + std::chrono::seconds{S};
        out._unit = datetime_unit::SECOND;
    }
    else if (digits.size() >= 12) {
        uint32_t Y = digits[0] * 1000 + digits[1] * 100 + digits[2] * 10 + digits[3];
        uint32_t m = digits[4] * 10 + digits[5];
        uint32_t d = digits[6] * 10 + digits[7];
        uint32_t H = digits[8] * 10 + digits[9];
        uint32_t M = digits[10] * 10 + digits[11];
        out._p = date::sys_days{date::year(Y) / date::month(m) / date::day(d)} +
                 std::chrono::hours{H} + std::chrono::minutes{M};
        out._unit = datetime_unit::MINUTE;
    }
    else if (digits.size() >= 10) {
        uint32_t Y = digits[0] * 1000 + digits[1] * 100 + digits[2] * 10 + digits[3];
        uint32_t m = digits[4] * 10 + digits[5];
        uint32_t d = digits[6] * 10 + digits[7];
        uint32_t H = digits[8] * 10 + digits[9];
        out._p = date::sys_days{date::year(Y) / date::month(m) / date::day(d)} +
                 std::chrono::hours{H} ;
        out._unit = datetime_unit::HOUR;
    }
    else if (digits.size() >= 8) {
        uint32_t Y = digits[0] * 1000 + digits[1] * 100 + digits[2] * 10 + digits[3];
        uint32_t m = digits[4] * 10 + digits[5];
        uint32_t d = digits[6] * 10 + digits[7];
        out._p = date::sys_days{date::year(Y) / date::month(m) / date::day(d)};
        out._unit = datetime_unit::DAY;
    }
    else if (digits.size() >= 6) {
        uint32_t Y = digits[0] * 1000 + digits[1] * 100 + digits[2] * 10 + digits[3];
        uint32_t m = digits[4] * 10 + digits[5];
        out._p = date::sys_days{date::year(Y) / date::month(m) / date::day(1)};
        out._unit = datetime_unit::MONTH;
    }
    else if (digits.size() >= 4) {
        uint32_t Y = digits[0] * 1000 + digits[1] * 100 + digits[2] * 10 + digits[3];
        out._p = date::sys_days{date::year(Y) / date::month(1) / date::day(1)};
        out._unit = datetime_unit::YEAR;
    }
    else {
        GCBS_DEBUG("Failed to parse dateetime from string '" + s + "'");
        out._unit = datetime_unit::NONE;
    }
    return out;
}

datetime datetime::from_string(std::string s) {
    return from_YmdHMS_digits(s); // 0.6.4 [experimental]
}

std::string datetime::datetime_format_for_unit(gdalcubes::datetime_unit u) {
    switch (u) {
        case datetime_unit::NONE:
        case datetime_unit::SECOND:
            return "%Y-%m-%dT%H:%M:%S";
        case datetime_unit::MINUTE:
            return "%Y-%m-%dT%H:%M:%S";
        case datetime_unit::HOUR:
            return "%Y-%m-%dT%H:%M:%S";
        case datetime_unit::DAY:
            return "%Y-%m-%d";
        case datetime_unit::WEEK:
            return "%Y-%m-%d";
        case datetime_unit::MONTH:
            return "%Y-%m-%d";
        case datetime_unit::YEAR:
            return "%Y-%m-%d";
    }
    return "%Y-%m-%dT%H:%M:%S";
}

}  // namespace gdalcubes
