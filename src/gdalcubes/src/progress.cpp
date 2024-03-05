/*
    MIT License

    Copyright (c) 2024 Marius Appel <marius.appel@hs-bochum.de>

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

#include <iostream>
#include "progress.h"


namespace gdalcubes {

void progress_simple_stdout::set(double p) {
    _m.lock();
    _set(p);
    _m.unlock();
};
void progress_simple_stdout::increment(double dp) {
    _m.lock();
    _set(_p + dp);
    _m.unlock();
}
void progress_simple_stdout::finalize() {
    _m.lock();
    for (uint16_t i = 0; i < (((int)(100 * _p)) / 10); ++i) {
        std::cout << "=";
    }
    std::cout << ">| DONE." << std::endl;
    _m.unlock();
}

void progress_simple_stdout::_set(double p) {  // not synchronized
    _m.lock();
    _p = p;
    for (uint16_t i = 0; i < (((int)(100 * p)) / 10); ++i) {
        std::cout << "=";
    }
    std::cout << "> (" << std::round(100 * p) << "%)";
    std::cout << "\r";
    std::cout.flush();
    _m.unlock();
}


void progress_simple_stdout_with_time::set(double p)  {
    _m.lock();
    _set(p);
    _m.unlock();
}

void progress_simple_stdout_with_time::increment(double dp)  {
    _m.lock();
    _set(_p + dp);
    _m.unlock();
}
void progress_simple_stdout_with_time::finalize()  {
    _m.lock();
    for (uint16_t i = 0; i < (((int)(100 * _p)) / 10); ++i) {
        std::cout << "=";
    }
    std::cout << ">| DONE (" << _t->time() << "s)." << std::endl;
    _m.unlock();
}

void progress_simple_stdout_with_time::_set(double p) {  // not synchronized
    _p = p;
    for (uint16_t i = 0; i < (((int)(100 * p)) / 10); ++i) {
        std::cout << "=";
    }
    std::cout << "> (" << std::round(100 * p) << "%)";
    std::cout << "\r";
    std::cout.flush();
}








}  // namespace gdalcubes
