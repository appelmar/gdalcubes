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

#ifndef CUBE_FACTORY_H
#define CUBE_FACTORY_H

#include <map>
#include <memory>
#include "cube.h"

/**
 * @brief Factory to create (nested) cubes from its JSON representation
 */
struct cube_factory {
    /**
     * Create a cube from its JSON representation
     * @param j
     * @return
     */
    static std::shared_ptr<cube> create_from_json(nlohmann::json j);
};

#endif  //CUBE_FACTORY_H
