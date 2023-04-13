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

#include "filesystem.h"

#include <cpl_conv.h>
#include <cpl_string.h>
#include <cpl_vsi.h>

#include <vector>

namespace gdalcubes {

bool filesystem::exists(std::string p) {
    VSIStatBufL s;
    return VSIStatL(p.c_str(), &s) == 0;
}

bool filesystem::is_directory(std::string p) {
    VSIStatBufL s;
    if (VSIStatL(p.c_str(), &s) != 0)
        return false;  // File / directory does not exist
    return VSI_ISDIR(s.st_mode);
}

bool filesystem::is_regular_file(std::string p) {
    VSIStatBufL s;
    if (VSIStatL(p.c_str(), &s) != 0)
        return false;  // File / directory does not exist
    return VSI_ISREG(s.st_mode);
}

std::string filesystem::stem(std::string p) {
    return CPLGetBasename(p.c_str());
}

std::string filesystem::filename(std::string p) {
    return std::string(CPLGetFilename(p.c_str()));
}

std::string filesystem::extension(std::string p) {
    return std::string(CPLGetExtension(p.c_str()));
}

std::string filesystem::directory(std::string p) {
    return std::string(CPLGetPath(p.c_str()));
}

std::string filesystem::get_working_dir() {
    char* x = CPLGetCurrentDir();
    std::string p;
    if (x) {
        p = join(std::string(x), p);
        CPLFree(x);
    }
    return p;
}

std::string filesystem::make_absolute(std::string p) {
    if (CPLIsFilenameRelative(p.c_str())) {
        char* x = CPLGetCurrentDir();
        if (x) {
            p = join(std::string(x), p);
            CPLFree(x);
        }
    }
    return p;
}

std::string filesystem::parent(std::string p) {
    if (!is_directory(p)) {
        return directory(p);
    }
    return std::string(CPLGetPath(CPLCleanTrailingSlash(p.c_str())));
}

std::string filesystem::join(std::string p1, std::string p2) {
    return p1 + DIR_SEPARATOR + p2;
}

void filesystem::iterate_directory(std::string p, std::function<void(const std::string&)> f) {
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

void filesystem::iterate_directory_recursive(std::string p, std::function<void(const std::string&)> f) {
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

void filesystem::remove(std::string p) {
    if (is_directory(p)) {
        VSIRmdir(p.c_str());
    } else {
        VSIUnlink(p.c_str());
    }
}

void filesystem::mkdir(std::string p) {
    VSIMkdir(p.c_str(), 0777);
}

void filesystem::mkdir_recursive(std::string p) {
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

bool filesystem::is_relative(std::string p) {
    return CPLIsFilenameRelative(p.c_str()) != 0;
}

bool filesystem::is_absolute(std::string p) {
    return !is_relative(p);
}

std::string filesystem::get_tempdir() {
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
}

uint32_t filesystem::file_size(std::string p) {
    VSIStatBufL s;
    if (VSIStatL(p.c_str(), &s) != 0)
        return 0;  // File / directory does not exist
    return s.st_size;
}

void filesystem::move(std::string src, std::string dest) {
    CPLMoveFile(dest.c_str(), src.c_str());
}

void filesystem::copy(std::string src, std::string dest) {
    CPLCopyFile(dest.c_str(), src.c_str());
}

}  // namespace gdalcubes
