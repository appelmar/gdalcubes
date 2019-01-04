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
#ifndef PROGRESS_H
#define PROGRESS_H

#include <cmath>
#include <mutex>
#include "timer.h"
#include "utils.h"

/**
 * @brief Virtual base class for progress updates of long running processes
*/
struct progress {
    /**
     * Get a new progress object of the same class
     * @return a new progress object
     */
    virtual std::shared_ptr<progress> get() = 0;

    /**
    * Set the progress to a specific value in [0,1]
    * @param p progress, 1 means 100%
    */
    virtual void set(double p) = 0;

    /**
   * Increment progress by a value of db
   * @param dp progress increment
   */
    virtual void increment(double dp) = 0;

    /**
     * Finalize the progress update such as printing "DONE"
     */
    virtual void finalize(){};
};

/**
 * @brief Implementation of progress which ignores any progress updates
 * @see progress
*/
struct progress_none : public progress {
    std::shared_ptr<progress> get() override { return std::make_shared<progress_none>(); }
    void set(double p) override {}
    void increment(double dp) override {}
};

/**
 * @brief Implementation of progress that streams updates to stdout
 * @see progress
*/
struct progress_simple_stdout : public progress {
    std::shared_ptr<progress> get() override { return std::make_shared<progress_simple_stdout>(); }
    void set(double p) override {
        _m.lock();
        _set(p);
        _m.unlock();
    };
    void increment(double dp) override {
        _m.lock();
        _set(_p + dp);
        _m.unlock();
    }
    virtual void finalize() override {
        _m.lock();
        for (uint16_t i = 0; i < (((int)(100 * _p)) / 10); ++i) {
            std::cout << "=";
        }
        std::cout << ">| DONE." << std::endl;
        _m.unlock();
    }

    progress_simple_stdout() : _p(0) {}

   private:
    void _set(double p) {  // not synchronized
        _m.lock();
        _p = p;
        for (uint16_t i = 0; i < (((int)(100 * p)) / 10); ++i) {
            std::cout << "=";
        }
        std::cout << "> (" << std::round(100 * p) << "%)";
        std::cout << "\r";
        std::cout.flush();
        _m.unlock();
    };

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
    void set(double p) override {
        _m.lock();
        _set(p);
        _m.unlock();
    };

    void increment(double dp) override {
        _m.lock();
        _set(_p + dp);
        _m.unlock();
    }
    virtual void finalize() override {
        _m.lock();
        for (uint16_t i = 0; i < (((int)(100 * _p)) / 10); ++i) {
            std::cout << "=";
        }
        std::cout << ">| DONE (" << _t->time() << "s)." << std::endl;
        _m.unlock();
    }

    progress_simple_stdout_with_time() : _t(nullptr), _p(0) {
        _t = new timer();
    }
    ~progress_simple_stdout_with_time() {
        if (_t != nullptr) delete _t;
    }

   private:
    void _set(double p) {  // not synchronized
        _p = p;
        for (uint16_t i = 0; i < (((int)(100 * p)) / 10); ++i) {
            std::cout << "=";
        }
        std::cout << "> (" << std::round(100 * p) << "%)";
        std::cout << "\r";
        std::cout.flush();
    };

    timer* _t;
    std::mutex _m;
    double _p;
};

#endif  //PROGRESS_H
