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

#include "apply_pixel.h"




#ifdef USE_EXPRTK
#include "external/exprtk.hpp"

typedef exprtk::symbol_table<double> symbol_table_t;
typedef exprtk::expression<double> expression_t;
typedef exprtk::parser<double> parser_t;

template <typename T>
struct isnan_func : public exprtk::ifunction<T> {
    isnan_func() : exprtk::ifunction<T>(1) {}
    T operator()(const T& v1) { return T(std::isnan(double(v1))); }
};

template <typename T>
struct isfinite_func : public exprtk::ifunction<T> {
    isfinite_func() : exprtk::ifunction<T>(1) {}
    T operator()(const T& v1) { return T(std::isfinite(double(v1))); }
};

#else
#include "external/tinyexpr/tinyexpr.h"

#endif



std::shared_ptr<chunk_data> apply_pixel_cube::read_chunk(chunkid_t id) {


#ifdef USE_EXPRTK
    GCBS_DEBUG("apply_pixel_cube::read_chunk(" + std::to_string(id) + ")");

    if (id < 0 || id >= count_chunks())
        return std::shared_ptr<chunk_data>();  // chunk is outside of the view, we don't need to read anything.

    std::shared_ptr<chunk_data> out = std::make_shared<chunk_data>();

    // Parse expressions and create symbol table
    std::vector<double> values;
    values.resize(_in_cube->bands().count(), NAN);
    // WARNING: never change size of values from now, because symbol table references to particular elements
    // and changing the size of vectors might cause reallocation

    // TODO: add further variables like x, y, t to symbol table
    symbol_table_t symbol_table;
    for (uint16_t i = 0; i < _in_cube->bands().count(); ++i) {
        symbol_table.add_variable(_in_cube->bands().get(i).name, values[i]);
    }
    symbol_table.add_constants();
    symbol_table.add_constant("NAN", nan(""));
    isnan_func<double> f_isnan;
    isfinite_func<double> f_isfinite;
    symbol_table.add_function("isfinite", f_isfinite);

    std::vector<expression_t> expr;
    for (uint16_t i = 0; i < _bands.count(); ++i) {
        expr.push_back(expression_t());
        expr[i].register_symbol_table(symbol_table);
        parser_t parser;
        if (!parser.compile(_expr[i], expr[i])) {
            for (uint16_t ie = 0; ie < parser.error_count(); ++ie) {
                typedef exprtk::parser_error::type error_t;
                error_t error = parser.get_error(ie);
                std::string msg = "Cannot parse expression for " + _bands.get(i).name + " '" + _expr[i] + "':" + exprtk::parser_error::to_str(error.mode) + "(token '" + error.token.value + "', #" + std::to_string(error.token.position) + ") - " + error.diagnostic;
                GCBS_ERROR(msg);
                return out;
            }
        }
    }

    std::shared_ptr<chunk_data> in = _in_cube->read_chunk(id);
    out->size({_bands.count(), in->size()[1], in->size()[2], in->size()[3]});
    out->buf(std::calloc(_bands.count() * in->size()[1] * in->size()[2] * in->size()[3], sizeof(double)));

    // We do not need to fill with NAN because we can be sure that it is completeley filled from the input cube
    //double *begin = (double *)out->buf();
    //double *end = ((double *)out->buf()) + size_btyx[0] * size_btyx[1] * size_btyx[2] * size_btyx[3];
    // std::fill(begin, end, NAN);

    for (uint16_t outb = 0; outb < _bands.count(); ++outb) {
        std::vector<uint16_t> bidx;
        for (auto it = _band_usage[outb].begin(); it != _band_usage[outb].end(); ++it) {
            bidx.push_back(_in_cube->bands().get_index(*it));
        }
        for (uint32_t i = 0; i < in->size()[1] * in->size()[2] * in->size()[3]; ++i) {
            for (uint16_t inb = 0; inb < bidx.size(); ++inb) {
                values[bidx[inb]] = ((double*)in->buf())[bidx[inb] * in->size()[1] * in->size()[2] * in->size()[3] + i];
            }
            ((double*)out->buf())[outb * in->size()[1] * in->size()[2] * in->size()[3] + i] = expr[outb].value();
        }
    }

    return out;

#else
    GCBS_DEBUG("apply_pixel_cube::read_chunk(" + std::to_string(id) + ")");

    if (id < 0 || id >= count_chunks())
        return std::shared_ptr<chunk_data>();  // chunk is outside of the view, we don't need to read anything.

    std::shared_ptr<chunk_data> out = std::make_shared<chunk_data>();

    // Parse expressions and create symbol table
    std::vector<double> values;
    values.resize(_in_cube->bands().count(), NAN);
    // WARNING: never change size of values from now

    // TODO: add further variables like x, y, t to symbol table

    std::vector<te_variable> vars;
    for (uint16_t i = 0; i < _in_cube->bands().count(); ++i) {
        char * varname = new char [_in_cube->bands().get(i).name.length()+1];
        std::string temp_name = _in_cube->bands().get(i).name;
        std::transform(temp_name.begin(), temp_name.end(), temp_name.begin(), ::tolower);
        std::strcpy (varname, temp_name.c_str());
        vars.push_back({varname, &values[i]});
    }
//
//    symbol_table.add_constants();
//    symbol_table.add_constant("NAN", nan(""));
//    isnan_func<double> f_isnan;
//    isfinite_func<double> f_isfinite;
//    symbol_table.add_function("isfinite", f_isfinite);

    std::vector<te_expr*> expr;
    for (uint16_t i = 0; i < _bands.count(); ++i) {
        int err;

        te_expr* x = te_compile(_expr[i].c_str(), vars.data(), vars.size(), &err);
        if (!x) {
            std::string msg = "Cannot parse expression for " + _bands.get(i).name + " '" + _expr[i] + "': error at token " + std::to_string(err);
            GCBS_ERROR(msg);

            // free other expressions
            for (uint16_t j=0; j<expr.size(); ++j) {
                te_free(expr[j]);
            }
            return out;
        }
        else {
            expr.push_back(x);
        }
    }

    std::shared_ptr<chunk_data> in = _in_cube->read_chunk(id);
    out->size({_bands.count(), in->size()[1], in->size()[2], in->size()[3]});
    out->buf(std::calloc(_bands.count() * in->size()[1] * in->size()[2] * in->size()[3], sizeof(double)));

    // We do not need to fill with NAN because we can be sure that it is completeley filled from the input cube
    //double *begin = (double *)out->buf();
    //double *end = ((double *)out->buf()) + size_btyx[0] * size_btyx[1] * size_btyx[2] * size_btyx[3];
    // std::fill(begin, end, NAN);

    for (uint16_t outb = 0; outb < _bands.count(); ++outb) {
        std::vector<uint16_t> bidx;
        for (auto it = _band_usage[outb].begin(); it != _band_usage[outb].end(); ++it) {
            bidx.push_back(_in_cube->bands().get_index(*it));
        }
        for (uint32_t i = 0; i < in->size()[1] * in->size()[2] * in->size()[3]; ++i) {
            for (uint16_t inb = 0; inb < bidx.size(); ++inb) {
                values[bidx[inb]] = ((double*)in->buf())[bidx[inb] * in->size()[1] * in->size()[2] * in->size()[3] + i];
            }
            ((double*)out->buf())[outb * in->size()[1] * in->size()[2] * in->size()[3] + i] = te_eval(expr[outb]);
        }
    }

    // free expressions
    for (uint16_t j=0; j<expr.size(); ++j) {
        te_free(expr[j]);
    }
    // free varnames
    for (uint16_t i=0; i<vars.size(); ++i) {
        delete[] vars[i].name;
    }

    return out;

#endif
}



bool apply_pixel_cube::parse_expressions() {

#ifdef USE_EXPRTK
    bool res = true;
    std::vector<double> dummy_values;

    // TODO: what if variable / band names are not valid e.g. start with a number?
    symbol_table_t symbol_table;
    for (uint16_t i = 0; i < _in_cube->bands().count(); ++i) {
        dummy_values.push_back(1.0);
        symbol_table.add_variable(_in_cube->bands().get(i).name, dummy_values[i]);
    }
    symbol_table.add_constants();

    std::vector<expression_t> expr;
    for (uint16_t i = 0; i < _bands.count(); ++i) {
        expr.push_back(expression_t());
        expr[i].register_symbol_table(symbol_table);
        parser_t parser;
        if (!parser.compile(_expr[i], expr[i])) {
            res = false;
            for (uint16_t ie = 0; ie < parser.error_count(); ++ie) {
                typedef exprtk::parser_error::type error_t;
                error_t error = parser.get_error(ie);
                std::string msg = "Cannot parse expression for " + _bands.get(i).name + " '" + _expr[i] + "': " + exprtk::parser_error::to_str(error.mode) + "(token '" + error.token.value + "', #" + std::to_string(error.token.position) + ") - " + error.diagnostic;
                GCBS_ERROR(msg);
                // Continue anyway to process all expressions
            }
        }
    }
    return res;
#else
    bool res = true;
    std::vector<double> dummy_values;
    std::vector<te_variable> vars;
    for (uint16_t i = 0; i < _in_cube->bands().count(); ++i) {
        dummy_values.push_back(1.0);
        char * varname = new char [_in_cube->bands().get(i).name.length()+1];
        std::string temp_name = _in_cube->bands().get(i).name;
        std::transform(temp_name.begin(), temp_name.end(), temp_name.begin(), ::tolower);
        std::strcpy (varname, temp_name.c_str());
        vars.push_back({varname, &dummy_values[i]});
    }

    int err = 0;
    for (uint16_t i = 0; i < _bands.count(); ++i) {
        te_expr *expr = te_compile(_expr[i].c_str(), vars.data(), vars.size(), &err);
        if (expr){
            te_free(expr);
        }
        else {
            res = false;
            std::string msg = "Cannot parse expression for " + _bands.get(i).name + " '" + _expr[i] + "': error at token " + std::to_string(err);
            GCBS_ERROR(msg);
            // Continue anyway to process all expressions
        }
    }
    for (uint16_t i=0; i<vars.size(); ++i) {
        delete[] vars[i].name;
    }
    return res;

#endif


}