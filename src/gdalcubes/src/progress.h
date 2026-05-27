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
#ifndef PROGRESS_H
#define PROGRESS_H

#include <cmath>
#include <mutex>
#include <cstdint> // 2023-01-12: GCC 13 compatibility

#include "progress_virt.h"
#include "timer.h"
#include "utils.h"




namespace gdalcubes {



/**
 * @brief Implementation of progress that streams updates to stdout
 * @see progress
*/
struct progress_simple_stdout : public progress {
    std::shared_ptr<progress> get() override { return std::make_shared<progress_simple_stdout>(); }
    void set(double p) override;
    void increment(double dp) override;
    virtual void finalize() override;

    progress_simple_stdout() : _p(0) {}

   private:
    void _set(double p);

    std::mutex _m;
    double _p;
};

/**
 * @brief Implementation of progress that streams updates to stdout and measures execution time
 * @see progress
*/
struct progress_simple_stdout_with_time : public progress {
    std::shared_ptr<progress> get() override {
        return std::make_shared<progress_simple_stdout_with_time>();
    }
    void set(double p) override;
    void increment(double dp) override;
    virtual void finalize() override;

    progress_simple_stdout_with_time() : _t(nullptr), _p(0) {
        _t = new timer();
    }
    ~progress_simple_stdout_with_time() {
        if (_t != nullptr) delete _t;
    }

   private:
    void _set(double p);

    timer* _t;
    std::mutex _m;
    double _p;
};

}  // namespace gdalcubes

#endif  //PROGRESS_H
