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

#include "utils.h"

#include <functional>  // std::hash
#include <iomanip>
#include <mutex>
#include <random>
#include <sstream>
#include <string>

namespace gdalcubes {

std::string utils::generate_unique_filename(uint16_t n, std::string prefix, std::string suffix) {
    static std::random_device rd{};
    static std::mt19937 gen(rd());  //Standard mersenne_twister_engine seeded with rd()
    static const std::string LETTERS = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789";
    static std::uniform_int_distribution<> dis(0, LETTERS.length() - 1);
    static std::mutex mtx;
    mtx.lock();
    std::stringstream ss;
    for (uint16_t i = 0; i < n; ++i) {
        ss << LETTERS[dis(gen)];
    }
    std::string out = prefix + ss.str() + suffix;
    mtx.unlock();
    return out;
}

std::string utils::get_curdatetime() {
    // Current date/time based on current system
    time_t now = time(0);

    // Convert now to tm struct for local timezone
    tm *localtm = localtime(&now);

    std::stringstream out;
    out << (localtm->tm_year + 1900) << "-" << std::setfill('0') << std::setw(2)
        << (localtm->tm_mon + 1) << "-" << std::setfill('0') << std::setw(2)
        << localtm->tm_mday << " " << std::setfill('0') << std::setw(2)
        << localtm->tm_hour << ":" << std::setfill('0') << std::setw(2)
        << localtm->tm_min << ":" << std::setfill('0') << std::setw(2)
        << localtm->tm_sec;
    return out.str();
}

std::string utils::get_curdate() {
    // Current date/time based on current system
    time_t now = time(0);

    // Convert now to tm struct for local timezone
    tm *localtm = localtime(&now);

    std::stringstream out;
    out << (localtm->tm_year + 1900) << "-" << std::setfill('0') << std::setw(2)
        << (localtm->tm_mon + 1) << "-" << std::setfill('0') << std::setw(2)
        << localtm->tm_mday;
    return out.str();
}

GDALDataType utils::gdal_type_from_string(std::string s) {
    if (s == "int16") return GDALDataType::GDT_Int16;
    if (s == "int32") return GDALDataType::GDT_Int32;
    if (s == "uint8") return GDALDataType::GDT_Byte;
    if (s == "uint16") return GDALDataType::GDT_UInt16;
    if (s == "uint32") return GDALDataType::GDT_UInt32;
    if (s == "float64") return GDALDataType::GDT_Float64;
    if (s == "float32") return GDALDataType::GDT_Float32;
    return GDALDataType::GDT_Unknown;
}

std::string utils::string_from_gdal_type(GDALDataType t) {
    switch (t) {
        case GDT_Float64:
            return "float64";
        case GDT_Float32:
            return "float32";
        case GDT_Int16:
            return "int16";
        case GDT_Int32:
            return "int32";
        case GDT_UInt32:
            return "uint32";
        case GDT_UInt16:
            return "uint16";
        case GDT_Byte:
            return "uint8";
        default:
            return "null";
    }
}

std::string utils::dbl_to_string(double x, uint8_t precision) {
    std::ostringstream ss;
    ss << std::fixed;
    ss << std::setprecision(precision);
    ss << x;
    return ss.str();
}

std::string utils::hash(std::string in) {
    std::size_t str_hash = std::hash<std::string>{}(in);
    return std::to_string(str_hash);
}

void utils::env::set(std::map<std::string, std::string> vars) {
    for (auto it = vars.begin(); it != vars.end(); ++it) {
#ifdef _WIN32
        _putenv((it->first + "=" + it->second).c_str());
#else
        setenv(it->first.c_str(), it->second.c_str(), 1);
#endif
        _vars[it->first] = it->second;
    }
}

void utils::env::unset(std::set<std::string> var_names) {
    for (auto it = var_names.begin(); it != var_names.end(); ++it) {
#ifdef _WIN32
        _putenv((*it + "=").c_str());
#else
        unsetenv((*it).c_str());
#endif
        auto i = _vars.find(*it);
        if (i != _vars.end()) {
            _vars.erase(i);
        }
    }
}

void utils::env::unset_all() {
    for (auto it = _vars.begin(); it != _vars.end(); ++it) {
#ifdef _WIN32
        _putenv((it->first + "=").c_str());
#else
        unsetenv(it->first.c_str());
#endif
    }
    _vars.clear();
}


std::string utils::env::get(std::string var_name, std::string default_value) {
    std::string out = default_value;
    char* r = getenv(var_name.c_str());
    if (r) {
        out = r;
    }
    return out;
}

std::string utils::env::to_string() {
    std::string out;
    out = "{";
    for (auto it = _vars.begin(); it != _vars.end(); ++it) {
        out += "{\"" + it->first + "\":\"" + it->second + "\"},";
    }
    out[out.length()-1] = '}';
    return out;
}

}  // namespace gdalcubes
