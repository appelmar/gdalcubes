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
#ifndef TIMER_H
#define TIMER_H

#include <chrono>

namespace gdalcubes {

/**
 * A simple class to measure real elapsed time with high resolution.
 */
class timer {
   public:
    /**
     * Create a new timer and start a time measurement.
     */
    timer() { t = std::chrono::high_resolution_clock::now(); }

    /**
     * Start measuring time from now, it is not needed to call this function as the constructor automatically starts
     * a time measurement.
     */
    void start() {
        t = std::chrono::high_resolution_clock::now();
    }

    /**
     * Measure real elapsed time since  either the last call of start() or construction.
     * @return the real elapsed time in seconds
     */
    double time() {
        std::chrono::high_resolution_clock::time_point tnow = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(
            tnow - t);
        return time_span.count();
    }

   private:
    std::chrono::high_resolution_clock::time_point t;
};

}  // namespace gdalcubes

#endif  //TIMER_H
