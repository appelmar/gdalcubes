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

#ifndef CUBE_FACTORY_H
#define CUBE_FACTORY_H

#include <map>
#include <memory>

#include "cube.h"

namespace gdalcubes {

/**
 * @brief Factory to create (nested) cubes from its JSON representation
 */
class cube_factory {
   public:
    static cube_factory* instance() {
        static CG g;
        if (!_instance) {
            _instance = new cube_factory();
        }
        return _instance;
    }

    /**
     * Create a cube from its JSON representation
     * @param j
     * @return
     */
    std::shared_ptr<cube> create_from_json(json11::Json j);

    std::shared_ptr<cube> create_from_json_file(std::string path);

    /**
     * @brief Registers a cube type with a function to create objects of this type from a JSON description.
     *
     * @param type_name unique name for cube type
     * @param generator function to create an object from a json definition
     */
    void register_cube_type(std::string type_name, std::function<std::shared_ptr<cube>(json11::Json&)> generator);

    void register_default();

   private:
    static cube_factory* _instance;
    cube_factory() {
        register_default();
    }
    cube_factory(const cube_factory&);
    cube_factory& operator=(const cube_factory&);
    ~cube_factory() {}
    class CG {
       public:
        ~CG() {
            if (NULL != cube_factory::_instance) {
                delete cube_factory::_instance;
                cube_factory::_instance = NULL;
            }
        }
    };

    std::map<std::string, std::function<std::shared_ptr<cube>(json11::Json&)>> cube_generators;
};

}  // namespace gdalcubes

#endif  //CUBE_FACTORY_H
