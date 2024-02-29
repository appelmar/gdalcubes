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
#ifndef UTILS_H
#define UTILS_H

#include <gdal_priv.h>
#include <set>
#include <map>
#include <cstdint> // 2023-01-12: GCC 13 compatibility

namespace gdalcubes {

/**
 * @brief A utility class for commonly used functions
 */
class utils {
   public:
    /**
    * @brief Generate a unique random filename
    * @param n number of characters of the random part
    * @param prefix string to append before the random part
    * @param suffix string to append after the random part
    * @return filename string
    */
    static std::string generate_unique_filename(uint16_t n = 8, std::string prefix = "", std::string suffix = "");

    /**
     * @brief Get the current datetime
     * @see  http://stackoverflow.com/questions/997946/how-to-get-current-time-and-date-in-c
     * @return datetime string
     */
    static std::string get_curdatetime();

    /**
     * @brief Get the current date
     * @return date string
     */
    static std::string get_curdate();

    /**
     * Convert a type name to the corresponding GDAL data type
     * @param s type name
     * @return GDAL data type
     */
    static GDALDataType gdal_type_from_string(std::string s);

    /**
     * @brief Convert a GDAL data type to a string typename used in gdalcubes
     * @param t GDAL data type
     * @return type name string
     */
    static std::string string_from_gdal_type(GDALDataType t);

    static std::string dbl_to_string(double x, uint8_t precision = std::numeric_limits<double>::max_digits10);

    /**
     * A simple (noncryptographic) hash function for strings using std::hash
     * @param in input string
     * @return hashed string
     */
    static std::string hash(std::string in);

    /**
     * A simple class to manage setting / unsetting environment variables,
     * mainly used to set variables for child processes.
     */
    class env {
       public:
        static env& instance() {
            static env _instance;
            return _instance;
        }
        ~env() {}

        void set(std::map<std::string, std::string> vars);
        void unset(std::set<std::string> var_names);
        void unset_all();
        std::string get(std::string var_name, std::string default_value = "");

        // Convert environment variable map to a JSON string
        std::string to_string();

       private:

        std::map<std::string, std::string> _vars;

        env() : _vars() {}
        env( const env& );
        env & operator = (const env &);
    };


};

}  // namespace gdalcubes

#endif  //UTILS_H
