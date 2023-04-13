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

#ifndef JOIN_BANDS_H
#define JOIN_BANDS_H

#include "cube.h"

namespace gdalcubes {

/**
 * @brief A data cube that combines the bands of two identically-shaped data cubes
 */
class join_bands_cube : public cube {
   public:
    /**
     * @brief Create a data cube that combines the bands of two identically-shaped data cubes
     * @note This static creation method should preferably be used instead of the constructors as
     * the constructors will not set connections between cubes properly.
     * @param A input data cube
     * @param B input data cube
     * @return a shared pointer to the created data cube instance
     */
    //    static std::shared_ptr<join_bands_cube> create(std::shared_ptr<cube> A, std::shared_ptr<cube> B, std::string prefix_A = "A", std::string prefix_B = "B") {
    //        std::shared_ptr<join_bands_cube> out = std::make_shared<join_bands_cube>(A, B, prefix_A, prefix_B);
    //        A->add_child_cube(out);
    //        B->add_child_cube(out);
    //        out->add_parent_cube(A);
    //        out->add_parent_cube(B);
    //        return out;
    //    }

    /**
   * @brief Create a data cube that combines the bands of two or more identically-shaped data cubes
   * @note This static creation method should preferably be used instead of the constructors as
   * the constructors will not set connections between cubes properly.
   * @param in_cubes vector with input data cube
   * @param prefixes vector with prefixes to name bands in the output cube
   * @return a shared pointer to the created data cube instance
   */
    static std::shared_ptr<join_bands_cube> create(std::vector<std::shared_ptr<cube>> in_cubes, std::vector<std::string> prefixes = {}) {
        std::shared_ptr<join_bands_cube> out = std::make_shared<join_bands_cube>(in_cubes, prefixes);
        for (auto it = in_cubes.begin(); it != in_cubes.end(); ++it) {
            (*it)->add_child_cube(out);
            out->add_parent_cube(*it);
        }
        return out;
    }

   public:
    join_bands_cube(std::vector<std::shared_ptr<cube>> in_cubes, std::vector<std::string> prefixes = {}) : cube(), _in(in_cubes), _prefix(prefixes) {
        _st_ref = std::make_shared<cube_stref_regular>();

        if (_in.size() < 2) {
            throw std::string("ERROR in join_bands_cube::join_bands_cube(): Expected at least two input data cubes");
        }

        if (!_prefix.empty() && (in_cubes.size() != prefixes.size())) {
            throw std::string("ERROR in join_bands_cube::join_bands_cube(): The number of name prefixes does not match the number of provided input data cubes");
        }

        for (uint16_t i = 1; i < in_cubes.size(); ++i) {
            if (cube_stref::type_string(_in[0]->st_reference()) != cube_stref::type_string(_in[i]->st_reference())) {
                throw std::string("ERROR in join_bands_cube::join_bands_cube(): Incompatible spatial / temporal reference types");
            }

            if (cube_stref::type_string(_in[0]->st_reference()) == "cube_stref_regular") {
                std::shared_ptr<cube_stref_regular> stref_A = std::dynamic_pointer_cast<cube_stref_regular>(_in[0]->st_reference());
                std::shared_ptr<cube_stref_regular> stref_B = std::dynamic_pointer_cast<cube_stref_regular>(_in[i]->st_reference());
                // Check that A and B have identical shape
                if (*(stref_A) != *(stref_B)) {
                    throw std::string("ERROR in join_bands_cube::join_bands_cube(): Cubes have different shape");
                }
            } else if (cube_stref::type_string(_in[0]->st_reference()) == "cube_stref_labeled_time") {
                std::shared_ptr<cube_stref_labeled_time> stref_A = std::dynamic_pointer_cast<cube_stref_labeled_time>(_in[0]->st_reference());
                std::shared_ptr<cube_stref_labeled_time> stref_B = std::dynamic_pointer_cast<cube_stref_labeled_time>(_in[i]->st_reference());
                // Check that A and B have identical shape
                if (*(stref_A) != *(stref_B)) {
                    throw std::string("ERROR in join_bands_cube::join_bands_cube(): Cubes have different shape");
                }
            }

            if (!(_in[0]->chunk_size()[0] == _in[i]->chunk_size()[0] &&
                  _in[0]->chunk_size()[1] == _in[i]->chunk_size()[1] &&
                  _in[0]->chunk_size()[2] == _in[i]->chunk_size()[2])) {
                throw std::string("ERROR in join_bands_cube::join_bands_cube(): Cubes have different chunk sizes");
            }

            if (_prefix.empty()) {
                bool has_name_conflicts = false;
                // check that there are no name conflicts
                for (uint16_t iba = 0; iba < _in[0]->bands().count(); ++iba) {
                    if (_in[i]->bands().has(_in[i]->bands().get(iba).name)) {
                        has_name_conflicts = true;
                        break;
                    }
                }
                if (has_name_conflicts) {
                    GCBS_WARN("Input cubes have bands with identical names, default name prefixes will be added.");
                    for (uint16_t j = 0; j < in_cubes.size(); ++j) {
                        _prefix.push_back("X" + std::to_string(j + 1));
                    }
                }
            } else {
                if (_prefix[0] == _prefix[i]) {
                    GCBS_ERROR("cannot join cubes with identical prefixes");
                    throw std::string("ERROR in join_bands_cube::join_bands_cube(): Cannot join cubes with identical prefixes");
                }
            }
        }

        _st_ref = _in[0]->st_reference();

        _chunk_size[0] = _in[0]->chunk_size()[0];
        _chunk_size[1] = _in[0]->chunk_size()[1];
        _chunk_size[2] = _in[0]->chunk_size()[2];

        for (uint16_t i = 0; i < in_cubes.size(); ++i) {
            for (uint16_t ib = 0; ib < _in[i]->bands().count(); ++ib) {
                band b = _in[i]->bands().get(ib);
                b.name = _prefix.empty() ? (b.name) : (_prefix[i] + "." + b.name);
                _bands.add(b);
            }
        }
    }

   public:
    ~join_bands_cube() {}

    std::shared_ptr<chunk_data> read_chunk(chunkid_t id) override;

    json11::Json make_constructible_json() override {
        json11::Json::object out;
        out["cube_type"] = "join_bands";
        json11::Json::array cubes;
        for (uint16_t i = 0; i < _in.size(); ++i) {
            cubes.push_back(_in[i]->make_constructible_json());
        }
        out["in_cubes"] = cubes;
        out["prefixes"] = _prefix;
        return out;
    }

   private:
    std::vector<std::shared_ptr<cube>> _in;
    std::vector<std::string> _prefix;
};

}  // namespace gdalcubes

#endif  // JOIN_BANDS_H
