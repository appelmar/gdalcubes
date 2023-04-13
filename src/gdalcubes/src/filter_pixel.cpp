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

#include "filter_pixel.h"

#include "external/tinyexpr/tinyexpr.h"

#include <cstring>


namespace gdalcubes {

std::shared_ptr<chunk_data> filter_pixel_cube::read_chunk(chunkid_t id) {
    GCBS_TRACE("filter_pixel_cube::read_chunk(" + std::to_string(id) + ")");

    if (id >= count_chunks())
        return  std::make_shared<chunk_data>();  // chunk is outside of the view, we don't need to read anything.

    std::shared_ptr<chunk_data> out = std::make_shared<chunk_data>();
    std::shared_ptr<chunk_data> in = _in_cube->read_chunk(id);
    if (in->empty()) {
        return out;
    }

    // Parse expressions and create symbol table
    std::vector<double> values;
    values.resize(_in_cube->bands().count(), NAN);
    // WARNING: never change size of values from now

    // TODO: add further variables like x, y, t to symbol table

    std::vector<te_variable> vars;
    for (uint16_t i = 0; i < _in_cube->bands().count(); ++i) {
        char* varname = new char[_in_cube->bands().get(i).name.length() + 1];
        std::string temp_name = _in_cube->bands().get(i).name;
        std::transform(temp_name.begin(), temp_name.end(), temp_name.begin(), ::tolower);
        std::strncpy(varname, temp_name.c_str(), temp_name.length() + 1);
        vars.push_back({varname, &values[i]});
    }

    int err;
    te_expr* expr = te_compile(_pred.c_str(), vars.data(), vars.size(), &err);
    if (!expr) {
        std::string msg = "Cannot parse predicate '" + _pred + "': error at token " + std::to_string(err);
        GCBS_ERROR(msg);
        te_free(expr);
        return out;
    }

    out->size({_bands.count(), in->size()[1], in->size()[2], in->size()[3]});
    out->buf(std::calloc(_bands.count() * in->size()[1] * in->size()[2] * in->size()[3], sizeof(double)));

    // We do not need to fill with NAN because we can be sure that it is completeley filled from the input cube
    //double *begin = (double *)out->buf();
    //double *end = ((double *)out->buf()) + size_btyx[0] * size_btyx[1] * size_btyx[2] * size_btyx[3];
    // std::fill(begin, end, NAN);

    for (uint32_t i = 0; i < in->size()[1] * in->size()[2] * in->size()[3]; ++i) {
        for (uint16_t inb = 0; inb < in->size()[0]; ++inb) {
            values[inb] = ((double*)in->buf())[inb * in->size()[1] * in->size()[2] * in->size()[3] + i];
        }
        if (te_eval(expr) != 0) {
            for (uint16_t ib = 0; ib < _bands.count(); ++ib) {
                ((double*)out->buf())[ib * in->size()[1] * in->size()[2] * in->size()[3] + i] = ((double*)in->buf())[ib * in->size()[1] * in->size()[2] * in->size()[3] + i];
            }
        } else {
            for (uint16_t ib = 0; ib < _bands.count(); ++ib) {
                ((double*)out->buf())[ib * in->size()[1] * in->size()[2] * in->size()[3] + i] = NAN;
            }
        }
    }

    te_free(expr);

    for (uint16_t i = 0; i < _in_cube->bands().count(); ++i) {  // only free name of actual variables (not functions)
        delete[] vars[i].name;
    }
    // check if chunk is completely NAN and if yes, return empty chunk
    if (out->all_nan()) {
        out = std::make_shared<chunk_data>();
    }

    return out;
}

bool filter_pixel_cube::parse_predicate() {
    bool res = true;
    std::vector<double> dummy_values;
    std::vector<::te_variable> vars;
    for (uint16_t i = 0; i < _in_cube->bands().count(); ++i) {
        dummy_values.push_back(1.0);
        char* varname = new char[_in_cube->bands().get(i).name.length() + 1];
        std::string temp_name = _in_cube->bands().get(i).name;
        std::transform(temp_name.begin(), temp_name.end(), temp_name.begin(), ::tolower);
        std::strncpy(varname, temp_name.c_str(), temp_name.length() + 1);
        vars.push_back({varname, &dummy_values[i]});
    }

    int err = 0;

    te_expr* expr = te_compile(_pred.c_str(), vars.data(), vars.size(), &err);
    if (expr) {
        te_free(expr);
    } else {
        res = false;
        std::string msg = "Cannot parse predicate'" + _pred + "': error at token " + std::to_string(err);
        GCBS_ERROR(msg);
        // Continue anyway to process all expressions
    }
    for (uint16_t i = 0; i < _in_cube->bands().count(); ++i) {  // only free name of actual variables (not functions)
        delete[] vars[i].name;
    }
    return res;
}

}  // namespace gdalcubes
