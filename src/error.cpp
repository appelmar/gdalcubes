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

#include "error.h"

#include "config.h"

std::mutex logger::_m;
// TODO: move mutex to specific error handler implementations
void logger::error(std::string msg, std::string where, int error_code) {
    _m.lock();
    config::instance()->get_error_handler()(error_level::ERRLVL_ERROR, msg, where, error_code);
    _m.unlock();
}
void logger::warn(std::string msg, std::string where, int error_code) {
    _m.lock();
    config::instance()->get_error_handler()(error_level::ERRLVL_WARNING, msg, where, error_code);
    _m.unlock();
}
void logger::debug(std::string msg, std::string where, int error_code) {
    _m.lock();
    config::instance()->get_error_handler()(error_level::ERRLVL_DEBUG, msg, where, error_code);
    _m.unlock();
}
void logger::fatal(std::string msg, std::string where, int error_code) {
    _m.lock();
    config::instance()->get_error_handler()(error_level::ERRLVL_FATAL, msg, where, error_code);
    _m.unlock();
}
void logger::trace(std::string msg, std::string where, int error_code) {
    _m.lock();
    config::instance()->get_error_handler()(error_level::ERRLVL_TRACE, msg, where, error_code);
    _m.unlock();
}
void logger::info(std::string msg, std::string where, int error_code) {
    _m.lock();
    config::instance()->get_error_handler()(error_level::ERRLVL_INFO, msg, where, error_code);
    _m.unlock();
}
