/*
    MIT License

    Copyright (c) 2021 Marius Appel <marius.appel@hs-bochum.de>

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

#ifndef RENAME_BANDS_H
#define RENAME_BANDS_H

#include "cube.h"

namespace gdalcubes {

/**
 * @brief A data cube that renames some ar all bands of an input cube
 */
class rename_bands_cube : public cube {
   public:
    /**
     * @brief Create a data cube that renames some ar all bands of an input cube
     * @note This static creation method should preferably be used instead of the constructors as
     * the constructors will not set connections between cubes properly.
     * @param in input data cube
     * @param bands_names map containing pairs (old_name, new_name)
     * @return a shared pointer to the created data cube instance
     */
    static std::shared_ptr<rename_bands_cube> create(std::shared_ptr<cube> in, std::map<std::string, std::string> band_names) {
        std::shared_ptr<rename_bands_cube> out = std::make_shared<rename_bands_cube>(in, band_names);
        in->add_child_cube(out);
        out->add_parent_cube(in);
        return out;
    }

   public:
    rename_bands_cube(std::shared_ptr<cube> in, std::map<std::string, std::string> band_names) : cube(in->st_reference()->copy()), _in_cube(in), _band_names(band_names) {  // it is important to duplicate st reference here, otherwise changes will affect input cube as well
        _chunk_size[0] = _in_cube->chunk_size()[0];
        _chunk_size[1] = _in_cube->chunk_size()[1];
        _chunk_size[2] = _in_cube->chunk_size()[2];

        for (uint16_t i = 0; i < in->bands().count(); ++i) {
            std::string oldname = in->bands().get(i).name;
            if (band_names.find(oldname) != band_names.end()) {
                band b = in->bands().get(i);
                b.name = band_names[oldname];
                _bands.add(b);
                band_names.erase(band_names.find(oldname));

            } else {
                _bands.add(in->bands().get(i));
            }
        }
        // Warnings for invalid bands
        for (auto it = band_names.begin(); it != band_names.end(); ++it) {
            GCBS_WARN("Input cube has no band with name '" + it->first + "'");
        }
    }

   public:
    ~rename_bands_cube() {}

    std::shared_ptr<chunk_data> read_chunk(chunkid_t id) override;

    json11::Json make_constructible_json() override {
        json11::Json::object band_names_json;
        for (auto it = _band_names.begin(); it != _band_names.end(); ++it) {
            band_names_json[it->first] = it->second;
        }
        json11::Json::object out;
        out["cube_type"] = "rename_bands";
        out["band_names"] = band_names_json;
        out["in_cube"] = _in_cube->make_constructible_json();
        return out;
    }

   private:
    std::shared_ptr<cube> _in_cube;
    std::map<std::string, std::string> _band_names;
};

}  // namespace gdalcubes

#endif  //RENAME_BANDS_H
