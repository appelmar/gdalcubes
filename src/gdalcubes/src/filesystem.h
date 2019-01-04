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

#ifndef FILESYSTEM_H
#define FILESYSTEM_H

#include <cpl_conv.h>
#include <cpl_string.h>
#include <cpl_vsi.h>
#include <functional>
#include <string>

#ifdef _WIN32
#define DIR_SEPARATOR "\\"
#else
#define DIR_SEPARATOR "/"
#endif

/**
 * @brief A simple wrapper class around GDAL's CPL and VSI interfaces for simple filesystem operations
 */
class filesystem {
   public:
    static bool exists(std::string p) {
        VSIStatBufL s;
        return VSIStatL(p.c_str(), &s) == 0;
    }

    static bool is_directory(std::string p) {
        VSIStatBufL s;
        if (VSIStatL(p.c_str(), &s) != 0)
            return false;  // File / directory does not exist
        return VSI_ISDIR(s.st_mode);
    }
    static bool is_regular_file(std::string p) {
        VSIStatBufL s;
        if (VSIStatL(p.c_str(), &s) != 0)
            return false;  // File / directory does not exist
        return VSI_ISREG(s.st_mode);
    }

    static std::string stem(std::string p) {
        return CPLGetBasename(p.c_str());
    }

    static std::string filename(std::string p) {
        return std::string(CPLGetFilename(p.c_str()));
    }
    static std::string extension(std::string p) {
        return std::string(CPLGetExtension(p.c_str()));
    }

    static std::string directory(std::string p) {
        return std::string(CPLGetPath(p.c_str()));
    }

    static std::string get_working_dir() {
        char* x = CPLGetCurrentDir();
        std::string p;
        if (x) {
            p = join(std::string(x), p);
            CPLFree(x);
        }
        return p;
    }

    static std::string make_absolute(std::string p) {
        if (CPLIsFilenameRelative(p.c_str())) {
            char* x = CPLGetCurrentDir();
            if (x) {
                p = join(std::string(x), p);
                CPLFree(x);
            }
        }
        return p;
    }

    static std::string parent(std::string p) {
        if (!is_directory(p)) {
            return directory(p);
        }
        return std::string(CPLGetPath(CPLCleanTrailingSlash(p.c_str())));
    }

    static std::string join(std::string p1, std::string p2) {
        return p1 + DIR_SEPARATOR + p2;
    }

    static void iterate_directory(std::string p, std::function<void(const std::string&)> f) {
        char** y = VSIReadDir(p.c_str());
        char** x = y;
        if (x != NULL) {
            while (*x != NULL) {
                f(join(p, std::string(*x)));
                ++x;
            }
            CSLDestroy(y);
        }
    }

    static void iterate_directory_recursive(std::string p, std::function<void(const std::string&)> f) {
        char** y = VSIReadDirRecursive(p.c_str());
        char** x = y;
        if (x != NULL) {
            while (*x != NULL) {
                f(join(p, std::string(*x)));
                ++x;
            }
            CSLDestroy(y);
        }
    }

    static void mkdir(std::string p) {
        VSIMkdir(p.c_str(), 0777);
    }

    static void mkdir_recursive(std::string p) {
        //VSIMkdirRecursive(p.c_str(), 0777); // available from GDAL 2.3

        if (p.empty()) return;

        if (is_directory(p)) {
            return;
        }

        std::string par = parent(p);

        if (par == p || par.length() >= p.length()) {
            return;
        }

        if (!exists(par)) {
            mkdir_recursive(par);
        }
        mkdir(p);
    }

    static bool is_relative(std::string p) {
        return CPLIsFilenameRelative(p.c_str()) != 0;
    }

    static bool is_absolute(std::string p) {
        return !is_relative(p);
    }

    static std::string get_tempdir() {
        std::vector<std::string> env_vars = {"TMPDIR", "TMP", "TEMP", "TEMPDIR", "USERPROFILE"};
        for (uint16_t i = 0; i < env_vars.size(); ++i) {
            if (std::getenv(env_vars[i].c_str()) != NULL) {
                if (filesystem::is_directory(std::getenv(env_vars[i].c_str())))
                    return std::string(std::getenv(env_vars[i].c_str()));
            }
        }

#ifdef _WIN32
        if (std::getenv("SYSTEMROOT") != NULL) {
            if (filesystem::is_directory(std::getenv("SYSTEMROOT")))
                return std::string(std::getenv("SYSTEMROOT"));
        }
        return "C:\\Windows";
#else
        return "/tmp";
#endif
    };

    static uint32_t file_size(std::string p) {
        VSIStatBufL s;
        if (VSIStatL(p.c_str(), &s) != 0)
            return 0;  // File / directory does not exist
        return s.st_size;
    }
};

#endif  //FILESYSTEM_H
