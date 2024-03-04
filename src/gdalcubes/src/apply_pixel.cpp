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

#include "apply_pixel.h"

#include "external/tinyexpr/tinyexpr.h"

#include <cstring>

namespace gdalcubes {

std::shared_ptr<chunk_data> apply_pixel_cube::read_chunk(chunkid_t id) {
    GCBS_TRACE("apply_pixel_cube::read_chunk(" + std::to_string(id) + ")");

    if (id >= count_chunks())
        return std::make_shared<chunk_data>();  // chunk is outside of the view, we don't need to read anything.

    std::shared_ptr<chunk_data> out = std::make_shared<chunk_data>();
    std::shared_ptr<chunk_data> in = _in_cube->read_chunk(id);

    out->set_status(in->status());  // propagate chunk status
    if (in->empty()) {
        return out;
    }

    // Parse expressions and create symbol table
    std::vector<double> values;
    values.resize(_in_cube->bands().count() + 9, NAN);
    // CAUTION: never change size of values from now

    // Add further variables like x, y, t to symbol table
    std::vector<te_variable> vars;
    for (uint16_t i = 0; i < _in_cube->bands().count(); ++i) {
        char* varname = new char[_in_cube->bands().get(i).name.length() + 1];
        std::string temp_name = _in_cube->bands().get(i).name;
        std::transform(temp_name.begin(), temp_name.end(), temp_name.begin(), ::tolower);
        std::strncpy(varname, temp_name.c_str(), temp_name.length() + 1);
        vars.push_back({varname, &values[i]});
    }
    vars.push_back({"t0", &values[_in_cube->bands().count() + 0]});
    vars.push_back({"t1", &values[_in_cube->bands().count() + 1]});
    vars.push_back({"left", &values[_in_cube->bands().count() + 2]});
    vars.push_back({"right", &values[_in_cube->bands().count() + 3]});
    vars.push_back({"top", &values[_in_cube->bands().count() + 4]});
    vars.push_back({"bottom", &values[_in_cube->bands().count() + 5]});
    vars.push_back({"ix", &values[_in_cube->bands().count() + 6]});
    vars.push_back({"iy", &values[_in_cube->bands().count() + 7]});
    vars.push_back({"it", &values[_in_cube->bands().count() + 8]});

    std::vector<te_expr*> expr;
    for (uint16_t i = 0; i < _expr.size(); ++i) {
        int err;

        te_expr* x = te_compile(_expr[i].c_str(), vars.data(), vars.size(), &err);
        if (!x) {
            std::string msg = "Cannot parse expression for " + _bands.get(i).name + " '" + _expr[i] + "': error at token " + std::to_string(err);
            GCBS_ERROR(msg);

            // free other expressions
            for (uint16_t j = 0; j < expr.size(); ++j) {
                te_free(expr[j]);
            }
            return out;
        } else {
            expr.push_back(x);
        }
    }

    out->size({_bands.count(), in->size()[1], in->size()[2], in->size()[3]});
    out->buf(std::calloc(_bands.count() * in->size()[1] * in->size()[2] * in->size()[3], sizeof(double)));

    // We do not need to fill with NAN because we can be sure that it is completly filled from the input cube
    //double *begin = (double *)out->buf();
    //double *end = ((double *)out->buf()) + size_btyx[0] * size_btyx[1] * size_btyx[2] * size_btyx[3];
    // std::fill(begin, end, NAN);

    if (_keep_bands) {
        std::memcpy(out->buf(), in->buf(),
                    sizeof(double) * in->size()[0] * in->size()[1] * in->size()[2] * in->size()[3]);
    }

    uint16_t outb = (_keep_bands) ? _in_cube->size_bands() : 0;
    uint16_t expr_idx = 0;
    while (outb < _bands.count()) {
        std::vector<uint16_t> bidx;
        for (auto it = _band_usage[expr_idx].begin(); it != _band_usage[expr_idx].end(); ++it) {
            bidx.push_back(_in_cube->bands().get_index(*it));
        }
        for (uint32_t i = 0; i < in->size()[1] * in->size()[2] * in->size()[3]; ++i) {
            for (uint16_t inb = 0; inb < bidx.size(); ++inb) {
                values[bidx[inb]] = ((double*)in->buf())[bidx[inb] * in->size()[1] * in->size()[2] * in->size()[3] + i];
            }

            // additional variables
            if (_var_usage[expr_idx].count("it") || _var_usage[expr_idx].count("t0") || _var_usage[expr_idx].count("t1")) {
                values[_in_cube->bands().count() + 8] = (double)(_in_cube->chunk_limits(id).low[0] + (i / (in->size()[2] * in->size()[3])));  // _t
            }
            if (_var_usage[expr_idx].count("t0")) {
                values[_in_cube->bands().count() + 0] = (_in_cube->st_reference()->datetime_at_index((int)(values[_in_cube->bands().count() + 8]))).epoch_time();
            }
            if (_var_usage[expr_idx].count("t1")) {
                values[_in_cube->bands().count() + 1] = (_in_cube->st_reference()->datetime_at_index((int)(values[_in_cube->bands().count() + 8] + 1))).epoch_time();
            }

            if (_var_usage[expr_idx].count("ix") || _var_usage[expr_idx].count("left") || _var_usage[expr_idx].count("right")) {
                values[_in_cube->bands().count() + 6] = (double)(_in_cube->chunk_limits(id).low[2] + (i % in->size()[3]));
            }
            if (_var_usage[expr_idx].count("iy") || _var_usage[expr_idx].count("top") || _var_usage[expr_idx].count("bottom")) {
                // values[_in_cube->bands().count() + 7] = (double)(_in_cube->size_y() - 1 - (_in_cube->chunk_limits(id).high[1] - ((i / in->size()[3]) % in->size()[2])));
                values[_in_cube->bands().count() + 7] = (double)(_in_cube->chunk_limits(id).low[1] + ((i / in->size()[3]) % in->size()[2]));
            }

            if (_in_cube->st_reference()->has_regular_space()) {
                if (_var_usage[expr_idx].count("left")) {
                    values[_in_cube->bands().count() + 2] = _in_cube->st_reference()->left() +
                                                            _in_cube->st_reference()->dx() * values[_in_cube->bands().count() + 6];
                }
                if (_var_usage[expr_idx].count("right")) {
                    values[_in_cube->bands().count() + 3] = _in_cube->st_reference()->left() +
                                                            _in_cube->st_reference()->dx() * (values[_in_cube->bands().count() + 6] + 1);
                }
                if (_var_usage[expr_idx].count("top")) {
                    values[_in_cube->bands().count() + 4] = _in_cube->st_reference()->top() -
                                                            _in_cube->st_reference()->dy() * values[_in_cube->bands().count() + 7];
                }
                if (_var_usage[expr_idx].count("bottom")) {
                    values[_in_cube->bands().count() + 5] = _in_cube->st_reference()->top() -
                                                            _in_cube->st_reference()->dy() * (values[_in_cube->bands().count() + 7] + 1);
                }
            }
            ((double*)out->buf())[outb * in->size()[1] * in->size()[2] * in->size()[3] + i] = te_eval(expr[expr_idx]);
        }
        ++outb;
        ++expr_idx;
    }

    // free expressions
    for (uint16_t j = 0; j < expr.size(); ++j) {
        te_free(expr[j]);
    }
    // free varnames
    for (uint16_t i = 0; i < _in_cube->bands().count(); ++i) {
        delete[] vars[i].name;  // delete only names of band variables
    }

    return out;
}

bool apply_pixel_cube::parse_expressions() {
    bool res = true;
    std::vector<double> dummy_values;
    dummy_values.resize(_in_cube->bands().count() + 9, 1.0);
    std::vector<te_variable> vars;
    for (uint16_t i = 0; i < _in_cube->bands().count(); ++i) {
        char* varname = new char[_in_cube->bands().get(i).name.length() + 1];
        std::string temp_name = _in_cube->bands().get(i).name;
        std::transform(temp_name.begin(), temp_name.end(), temp_name.begin(), ::tolower);
        std::strncpy(varname, temp_name.c_str(), temp_name.length() + 1);
        vars.push_back({varname, &dummy_values[i]});
    }

    vars.push_back({"t0", &dummy_values[_in_cube->bands().count() + 0]});
    vars.push_back({"t1", &dummy_values[_in_cube->bands().count() + 1]});
    vars.push_back({"left", &dummy_values[_in_cube->bands().count() + 2]});
    vars.push_back({"right", &dummy_values[_in_cube->bands().count() + 3]});
    vars.push_back({"top", &dummy_values[_in_cube->bands().count() + 4]});
    vars.push_back({"bottom", &dummy_values[_in_cube->bands().count() + 5]});
    vars.push_back({"ix", &dummy_values[_in_cube->bands().count() + 6]});
    vars.push_back({"iy", &dummy_values[_in_cube->bands().count() + 7]});
    vars.push_back({"it", &dummy_values[_in_cube->bands().count() + 8]});

    int err = 0;
    for (uint16_t i = 0; i < _expr.size(); ++i) {
        te_expr* expr = te_compile(_expr[i].c_str(), vars.data(), vars.size(), &err);
        if (expr) {
            te_free(expr);
        } else {
            res = false;
            std::string msg = "Cannot parse expression for " + _bands.get(i).name + " '" + _expr[i] + "': error at token " + std::to_string(err);
            GCBS_ERROR(msg);
            // Continue anyway to process all expressions
        }
    }
    for (uint16_t i = 0; i < _in_cube->bands().count(); ++i) {
        delete[] vars[i].name;  // delete only names of band variables
    }
    return res;
}

}  // namespace gdalcubes
