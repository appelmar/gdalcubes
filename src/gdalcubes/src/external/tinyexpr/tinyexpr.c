/*
 * TINYEXPR - Tiny recursive descent parser and evaluation engine in C
 *
 * Copyright (c) 2015-2018 Lewis Van Winkle
 *
 * http://CodePlea.com
 *
 * This software is provided 'as-is', without any express or implied
 * warranty. In no event will the authors be held liable for any damages
 * arising from the use of this software.
 *
 * Permission is granted to anyone to use this software for any purpose,
 * including commercial applications, and to alter it and redistribute it
 * freely, subject to the following restrictions:
 *
 * 1. The origin of this software must not be misrepresented; you must not
 * claim that you wrote the original software. If you use this software
 * in a product, an acknowledgement in the product documentation would be
 * appreciated but is not required.
 * 2. Altered source versions must be plainly marked as such, and must not be
 * misrepresented as being the original software.
 * 3. This notice may not be removed or altered from any source distribution.
 */


/*
 * THIS FILE INCLUDES THE FOLLOWING MODIFICATIONS (c) 2019 Marius Appel:
 * - added logical infix operators <, <=, >, >=, ==, !=, &&, || and logical not !
 * - added bitwise infix operators &, |, <<, >> and bitwise not ~
 * - added isnan(), isfinite(), and iif() functions
 * - commented out pn() and te_print()
 */



/* COMPILE TIME OPTIONS */

/* Exponentiation associativity:
For a^b^c = (a^b)^c and -a^b = (-a)^b do nothing.
For a^b^c = a^(b^c) and -a^b = -(a^b) uncomment the next line.*/
/* #define TE_POW_FROM_RIGHT */

/* Logarithms
For log = base 10 log do nothing
For log = natural log uncomment the next line. */
/* #define TE_NAT_LOG */

#include "tinyexpr.h"
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdio.h>
#include <limits.h>
#include <time.h> // added by Marius Appel on 2019-09-10

#ifndef NAN
#define NAN (0.0/0.0)
#endif

#ifndef INFINITY
#define INFINITY (1.0/0.0)
#endif


typedef double (*te_fun2)(double, double);


enum {
    TOK_NULL = TE_CLOSURE7+1, TOK_ERROR, TOK_END, TOK_SEP,
    TOK_OPEN, TOK_CLOSE, TOK_NUMBER, TOK_VARIABLE, TOK_INFIX
};


enum {TE_CONSTANT = 1};


typedef struct state {
    const char *start;
    const char *next;
    int type;
    union te_binding binding;
    void *context;

    const te_variable *lookup;
    int lookup_len;
} state;


#define TYPE_MASK(TYPE) ((TYPE)&0x0000001F)

#define IS_PURE(TYPE) (((TYPE) & TE_FLAG_PURE) != 0)
#define IS_FUNCTION(TYPE) (((TYPE) & TE_FUNCTION0) != 0)
#define IS_CLOSURE(TYPE) (((TYPE) & TE_CLOSURE0) != 0)
#define ARITY(TYPE) ( ((TYPE) & (TE_FUNCTION0 | TE_CLOSURE0)) ? ((TYPE) & 0x00000007) : 0 )
#define NEW_EXPR(type, ...) new_expr((type), (const te_expr*[]){__VA_ARGS__})

static te_expr *new_expr(const int type, const te_expr *parameters[]) {
    const int arity = ARITY(type);
    const int psize = sizeof(void*) * arity;
    //const int size = (sizeof(te_expr) - sizeof(void*)) + psize + (IS_CLOSURE(type) ? sizeof(void*) : 0);

    // The following lines have been changed in order to make sure that all expressions 
    // have the same size and that pointers to parameters have correct size (to avoid UBSAN warnings)
    // (modified on 2025-03-18 by Marius Appel)
    // int size = (sizeof(te_expr) - sizeof(void*)) + psize + (IS_CLOSURE(type) ? sizeof(void*) : 0);
    // if (size < sizeof(te_expr)) {
    //   size = sizeof(te_expr);
    // }
    
    int size = sizeof(te_expr);

    te_expr *ret = malloc(size);
    
    //memset(ret, 0, size); 
    if (arity && parameters) {
        ret->parameters = malloc(psize);
        memcpy(ret->parameters, parameters, psize);
    }
    else {
      ret->parameters = malloc(sizeof(void*));
    }
    ret->type = type;
    ret->binding.bound = 0;
    return ret;
}


void te_free_parameters(te_expr *n) {
    if (!n) return;
    switch (TYPE_MASK(n->type)) {
        case TE_FUNCTION7: case TE_CLOSURE7: te_free(n->parameters[6]);     /* Falls through. */
        case TE_FUNCTION6: case TE_CLOSURE6: te_free(n->parameters[5]);     /* Falls through. */
        case TE_FUNCTION5: case TE_CLOSURE5: te_free(n->parameters[4]);     /* Falls through. */
        case TE_FUNCTION4: case TE_CLOSURE4: te_free(n->parameters[3]);     /* Falls through. */
        case TE_FUNCTION3: case TE_CLOSURE3: te_free(n->parameters[2]);     /* Falls through. */
        case TE_FUNCTION2: case TE_CLOSURE2: te_free(n->parameters[1]);     /* Falls through. */
        case TE_FUNCTION1: case TE_CLOSURE1: te_free(n->parameters[0]);
    }
}


void te_free(te_expr *n) {
    if (!n) return;
    te_free_parameters(n);
    if (n->parameters) free(n->parameters);
    free(n);
}


// new functions (added by Marius Appel on Jan 23, 2019)
static double lt(double a, double b) { return a < b; }
static double lte(double a, double b) { return a <= b; }
static double gte(double a, double b) { return a >= b; }
static double gt(double a, double b) { return a > b; }
static double eq(double a, double b) { return a == b; }
static double neq(double a, double b) { return a != b; }
static double land(double a, double b) { return (int)(a) && (int)(b); }
static double lor(double a, double b) { return (int)(a) || (int)(b); }
static double lnot(double a) { return !(int)(a); }
static double shr(double a, double b) { return (int)(a) >> (int)(b); }
static double shl(double a, double b) { return (int)(a) << (int)(b); }
static double band(double a, double b) { return (int)(a) & (int)(b); }
static double bor(double a, double b) { return (int)(a) | (int)(b); }
static double bnot(double a) { return ~(int)(a); }
static double is_nan(double a) { return isnan(a); }
static double is_finite(double a) { return isfinite(a); }
static double iif(double a, double b, double c) { return (int)(a)? b : c;}
// end of new functions

// added by Marius Appel on 2019-09-10
static double year(double t) {
    time_t tt = (long)t;
    struct tm * ptm;
    ptm = gmtime(&tt);
    return (double)(ptm->tm_year + 1900);
}

static double month(double t) {
    time_t tt = (long)t;
    struct tm * ptm;
    ptm = gmtime(&tt);
    return (double)ptm->tm_mon;
}

static double dom(double t) {
    time_t tt = (long)t;
    struct tm * ptm;
    ptm = gmtime(&tt);
    return (double)ptm->tm_mday;
}
static double doy(double t) {
    time_t tt = (long)t;
    struct tm * ptm;
    ptm = gmtime(&tt);
    return (double)ptm->tm_yday;
}
static double dow(double t) {
    time_t tt = (long)t;
    struct tm * ptm;
    ptm = gmtime(&tt);
    return (double)ptm->tm_wday;
}
static double hours(double t) {
    time_t tt = (long)t;
    struct tm * ptm;
    ptm = gmtime(&tt);
    return (double)ptm->tm_hour;
}
static double minutes(double t) {
    time_t tt = (long)t;
    struct tm * ptm;
    ptm = gmtime(&tt);
    return (double)ptm->tm_min;
}
static double seconds(double t) {
    time_t tt = (long)t;
    struct tm * ptm;
    ptm = gmtime(&tt);
    return (double)ptm->tm_sec;
}
// end of new functions

static double pi(void) {return 3.14159265358979323846;}
static double e(void) {return 2.71828182845904523536;}
static double fac(double a) {/* simplest version of fac */
    if (a < 0.0)
        return NAN;
    if (a > UINT_MAX)
        return INFINITY;
    unsigned int ua = (unsigned int)(a);
    unsigned long int result = 1, i;
    for (i = 1; i <= ua; i++) {
        if (i > ULONG_MAX / result)
            return INFINITY;
        result *= i;
    }
    return (double)result;
}
static double ncr(double n, double r) {
    if (n < 0.0 || r < 0.0 || n < r) return NAN;
    if (n > UINT_MAX || r > UINT_MAX) return INFINITY;
    unsigned long int un = (unsigned int)(n), ur = (unsigned int)(r), i;
    unsigned long int result = 1;
    if (ur > un / 2) ur = un - ur;
    for (i = 1; i <= ur; i++) {
        if (result > ULONG_MAX / (un - ur + i))
            return INFINITY;
        result *= un - ur + i;
        result /= i;
    }
    return result;
}
static double npr(double n, double r) {return ncr(n, r) * fac(r);}

static const te_function functions[] = {
    /* must be in alphabetical order */
    {"abs", (funcptr)fabs,     TE_FUNCTION1 | TE_FLAG_PURE, 0},
    {"acos", (funcptr)acos,    TE_FUNCTION1 | TE_FLAG_PURE, 0},
    {"asin", (funcptr)asin,    TE_FUNCTION1 | TE_FLAG_PURE, 0},
    {"atan", (funcptr)atan,    TE_FUNCTION1 | TE_FLAG_PURE, 0},
    {"atan2", (funcptr)atan2,  TE_FUNCTION2 | TE_FLAG_PURE, 0},
    {"ceil", (funcptr)ceil,    TE_FUNCTION1 | TE_FLAG_PURE, 0},
    {"cos", (funcptr)cos,      TE_FUNCTION1 | TE_FLAG_PURE, 0},
    {"cosh", (funcptr)cosh,    TE_FUNCTION1 | TE_FLAG_PURE, 0},
    {"day", (funcptr)dom,    TE_FUNCTION1 | TE_FLAG_PURE, 0}, // new function (added by Marius Appel on Sept 10, 2019)
    {"dayofyear", (funcptr)doy,    TE_FUNCTION1 | TE_FLAG_PURE, 0}, // new function (added by Marius Appel on Sept 10, 2019)
    {"e", (funcptr)e,          TE_FUNCTION0 | TE_FLAG_PURE, 0},
    {"exp", (funcptr)exp,      TE_FUNCTION1 | TE_FLAG_PURE, 0},
    {"fac", (funcptr)fac,      TE_FUNCTION1 | TE_FLAG_PURE, 0},
    {"floor", (funcptr)floor,  TE_FUNCTION1 | TE_FLAG_PURE, 0},
    {"hours", (funcptr)hours,    TE_FUNCTION1 | TE_FLAG_PURE, 0}, // new function (added by Marius Appel on Sept 10, 2019)
    {"ifelse", (funcptr)iif,    TE_FUNCTION3 | TE_FLAG_PURE, 0}, // new function (added by Marius Appel on Jan 23, 2019)
    {"isfinite", (funcptr)is_finite,    TE_FUNCTION1 | TE_FLAG_PURE, 0}, // new function (added by Marius Appel on Jan 23, 2019)
    {"isnan", (funcptr)is_nan,    TE_FUNCTION1 | TE_FLAG_PURE, 0}, // new function (added by Marius Appel on Jan 23, 2019)
    {"year", (funcptr)year,    TE_FUNCTION1 | TE_FLAG_PURE, 0}, // new function (added by Marius Appel on Sept 10, 2019)
    {"ln", (funcptr)log,       TE_FUNCTION1 | TE_FLAG_PURE, 0},
#ifdef TE_NAT_LOG
    {"log", (funcptr)log,      TE_FUNCTION1 | TE_FLAG_PURE, 0},
#else
    {"log", (funcptr)log10,    TE_FUNCTION1 | TE_FLAG_PURE, 0},
#endif
    {"log10", (funcptr)log10,  TE_FUNCTION1 | TE_FLAG_PURE, 0},
    {"minutes", (funcptr)minutes,    TE_FUNCTION1 | TE_FLAG_PURE, 0}, // new function (added by Marius Appel on Sept 10, 2019)
    {"month", (funcptr)month,    TE_FUNCTION1 | TE_FLAG_PURE, 0}, // new function (added by Marius Appel on Sept 10, 2019)
    {"ncr", (funcptr)ncr,      TE_FUNCTION2 | TE_FLAG_PURE, 0},
    {"npr", (funcptr)npr,      TE_FUNCTION2 | TE_FLAG_PURE, 0},
    {"pi", (funcptr)pi,        TE_FUNCTION0 | TE_FLAG_PURE, 0},
    {"pow", (funcptr)pow,      TE_FUNCTION2 | TE_FLAG_PURE, 0},
    {"seconds", (funcptr)seconds,    TE_FUNCTION1 | TE_FLAG_PURE, 0}, // new function (added by Marius Appel on Sept 10, 2019)
    {"sin", (funcptr)sin,      TE_FUNCTION1 | TE_FLAG_PURE, 0},
    {"sinh", (funcptr)sinh,    TE_FUNCTION1 | TE_FLAG_PURE, 0},
    {"sqrt", (funcptr)sqrt,    TE_FUNCTION1 | TE_FLAG_PURE, 0},
    {"tan", (funcptr)tan,      TE_FUNCTION1 | TE_FLAG_PURE, 0},
    {"tanh", (funcptr)tanh,    TE_FUNCTION1 | TE_FLAG_PURE, 0},
    {"weekday", (funcptr)dow,    TE_FUNCTION1 | TE_FLAG_PURE, 0}, // new function (added by Marius Appel on Sept 10, 2019)
    {"year", (funcptr)year,    TE_FUNCTION1 | TE_FLAG_PURE, 0}, // new function (added by Marius Appel on Sept 10, 2019)

    {0, 0, 0, 0}
};

static const te_function *find_builtin(const char *name, int len) {
    int imin = 0;
    int imax = sizeof(functions) / sizeof(te_function) - 2;

    /*Binary search.*/
    while (imax >= imin) {
        const int i = (imin + ((imax-imin)/2));
        int c = strncmp(name, functions[i].name, len);
        if (!c) c = '\0' - functions[i].name[len];
        if (c == 0) {
            return functions + i;
        } else if (c > 0) {
            imin = i + 1;
        } else {
            imax = i - 1;
        }
    }

    return 0;
}

static const te_variable *find_lookup(const state *s, const char *name, int len) {
    int iters;
    const te_variable *var;
    if (!s->lookup) return 0;

    for (var = s->lookup, iters = s->lookup_len; iters; ++var, --iters) {
        if (strncmp(name, var->name, len) == 0 && var->name[len] == '\0') {
            return var;
        }
    }
    return 0;
}



static double add(double a, double b) {return a + b;}
static double sub(double a, double b) {return a - b;}
static double mul(double a, double b) {return a * b;}
static double divide(double a, double b) {return a / b;}
static double negate(double a) {return -a;}
static double comma(double a, double b) {(void)a; return b;}


union te_symbol {
  const te_variable *var;
  const te_function *func;
};

void next_token(state *s) {
    s->type = TOK_NULL;

    do {

        if (!*s->next){
            s->type = TOK_END;
            return;
        }

        /* Try reading a number. */
        if ((s->next[0] >= '0' && s->next[0] <= '9') || s->next[0] == '.') {
            s->binding.value = strtod(s->next, (char**)&s->next);
            s->type = TOK_NUMBER;
        } else {
            /* Look for a variable or builtin function call. */
            if (s->next[0] >= 'a' && s->next[0] <= 'z') {
                const char *start;
                start = s->next;
                while ((s->next[0] >= 'a' && s->next[0] <= 'z') || (s->next[0] >= '0' && s->next[0] <= '9') || (s->next[0] == '_')) s->next++;

                union te_symbol var;
                var.var = find_lookup(s, start, s->next - start);
                if (!var.var) var.func = find_builtin(start, s->next - start);
                if (!var.func) {
                    s->type = TOK_ERROR;
                } else {
                    switch(TYPE_MASK(var.func->type))
                    {
                        case TE_VARIABLE:
                            s->type = TOK_VARIABLE;
                            s->binding.bound = var.var->address; // does this work?
                            break;

                        case TE_CLOSURE0: case TE_CLOSURE1: case TE_CLOSURE2: case TE_CLOSURE3:         /* Falls through. */
                        case TE_CLOSURE4: case TE_CLOSURE5: case TE_CLOSURE6: case TE_CLOSURE7:         /* Falls through. */
                            s->context = var.func->context;                                                  /* Falls through. */

                        case TE_FUNCTION0: case TE_FUNCTION1: case TE_FUNCTION2: case TE_FUNCTION3:     /* Falls through. */
                        case TE_FUNCTION4: case TE_FUNCTION5: case TE_FUNCTION6: case TE_FUNCTION7:     /* Falls through. */
                            s->type = var.func->type;
                            s->binding.function = var.func->address;
                            break;
                    }
                }

            } else {

                if (strncmp(s->next, "||", 2) == 0) {
                    s->type = TOK_INFIX; s->binding.function = (funcptr)lor;
                    s->next += 2;
                }
                else if (strncmp(s->next, "&&", 2) == 0) {
                    s->type = TOK_INFIX; s->binding.function = (funcptr)land;
                    s->next += 2;
                }
                else if (strncmp(s->next, "<=", 2) == 0) {
                    s->type = TOK_INFIX; s->binding.function = (funcptr)lte;
                    s->next += 2;
                }
                else if (strncmp(s->next, ">=", 2) == 0) {
                    s->type = TOK_INFIX; s->binding.function = (funcptr)gte;
                    s->next += 2;
                }
                else if (strncmp(s->next, "==", 2) == 0) {
                    s->type = TOK_INFIX; s->binding.function = (funcptr)eq;
                    s->next += 2;
                }
                else if (strncmp(s->next, "!=", 2) == 0) {
                    s->type = TOK_INFIX; s->binding.function = (funcptr)neq;
                    s->next += 2;
                }
                else if (strncmp(s->next, "<<", 2) == 0) {
                    s->type = TOK_INFIX; s->binding.function = (funcptr)shl;
                    s->next += 2;
                }
                else if (strncmp(s->next, ">>", 2) == 0) {
                    s->type = TOK_INFIX; s->binding.function = (funcptr)shr;
                    s->next += 2;
                }
                else if (strncmp(s->next, "+", 1) == 0) {
                    s->type = TOK_INFIX; s->binding.function = (funcptr)add;
                    s->next += 1;
                }
                else if (strncmp(s->next, "-", 1) == 0) {
                    s->type = TOK_INFIX; s->binding.function = (funcptr)sub;
                    s->next += 1;
                }
                else if (strncmp(s->next, "*", 1) == 0) {
                    s->type = TOK_INFIX; s->binding.function = (funcptr)mul;
                    s->next += 1;
                }
                else if (strncmp(s->next, "/", 1) == 0) {
                    s->type = TOK_INFIX; s->binding.function = (funcptr)divide;
                    s->next += 1;
                }
                else if (strncmp(s->next, "^", 1) == 0) {
                    s->type = TOK_INFIX; s->binding.function = (funcptr)pow;
                    s->next += 1;
                }
                else if (strncmp(s->next, "%", 1) == 0) {
                    s->type = TOK_INFIX; s->binding.function = (funcptr)fmod;
                    s->next += 1;
                }
                else if (strncmp(s->next, "(", 1) == 0) {
                    s->type = TOK_OPEN;
                    s->next += 1;
                }
                else if (strncmp(s->next, ")", 1) == 0) {
                    s->type = TOK_CLOSE;
                    s->next += 1;
                }
                else if (strncmp(s->next, ",", 1) == 0) {
                    s->type = TOK_SEP;
                    s->next += 1;
                }
                else if (s->next[0] == ' ' || s->next[0] == '\t' || s->next[0] == '\n' || s->next[0] == '\r' ) {
                    s->next += 1;
                }
                else if (strncmp(s->next, "<", 1) == 0) {
                    s->type = TOK_INFIX; s->binding.function = (funcptr)lt;
                    s->next += 1;
                }
                else if (strncmp(s->next, ">", 1) == 0) {
                    s->type = TOK_INFIX; s->binding.function = (funcptr)gt;
                    s->next += 1;
                }
                else if (strncmp(s->next, "|", 1) == 0) {
                    s->type = TOK_INFIX; s->binding.function = (funcptr)bor;
                    s->next += 1;
                }
                else if (strncmp(s->next, "&", 1) == 0) {
                    s->type = TOK_INFIX; s->binding.function = (funcptr)band;
                    s->next += 1;
                }
                else if (strncmp(s->next, "~", 1) == 0) {
                    s->type = TOK_INFIX; s->binding.function = (funcptr)bnot;
                    s->next += 1;
                }
                else if (strncmp(s->next, "!", 1) == 0) {
                    s->type = TOK_INFIX; s->binding.function = (funcptr)lnot;
                    s->next += 1;
                }



                else {
                    s->type = TOK_ERROR;
                }
            }
        }
    } while (s->type == TOK_NULL);
}


static te_expr *list(state *s);
static te_expr *expr(state *s);
static te_expr *power(state *s);

static te_expr *base(state *s) {
    /* <base>      =    <constant> | <variable> | <function-0> {"(" ")"} | <function-1> <power> | <function-X> "(" <expr> {"," <expr>} ")" | "(" <list> ")" */
    te_expr *ret;
    int arity;

    switch (TYPE_MASK(s->type)) {
        case TOK_NUMBER:
            ret = new_expr(TE_CONSTANT, 0);
            ret->binding.value = s->binding.value;
            next_token(s);
            break;

        case TOK_VARIABLE:
            ret = new_expr(TE_VARIABLE, 0);
            ret->binding.bound = s->binding.bound;
            next_token(s);
            break;

        case TE_FUNCTION0:
        case TE_CLOSURE0:
            ret = new_expr(s->type, 0);
            ret->binding.function = s->binding.function;
            if (IS_CLOSURE(s->type)) ret->parameters[0] = s->context;
            next_token(s);
            if (s->type == TOK_OPEN) {
                next_token(s);
                if (s->type != TOK_CLOSE) {
                    s->type = TOK_ERROR;
                } else {
                    next_token(s);
                }
            }
            break;

        case TE_FUNCTION1:
        case TE_CLOSURE1:
            ret = new_expr(s->type, 0);
            ret->binding.function = s->binding.function;
            if (IS_CLOSURE(s->type)) ret->parameters[1] = s->context;
            next_token(s);
            ret->parameters[0] = power(s);
            break;

        case TE_FUNCTION2: case TE_FUNCTION3: case TE_FUNCTION4:
        case TE_FUNCTION5: case TE_FUNCTION6: case TE_FUNCTION7:
        case TE_CLOSURE2: case TE_CLOSURE3: case TE_CLOSURE4:
        case TE_CLOSURE5: case TE_CLOSURE6: case TE_CLOSURE7:
            arity = ARITY(s->type);

            ret = new_expr(s->type, 0);
            ret->binding.function = s->binding.function;
            if (IS_CLOSURE(s->type)) ret->parameters[arity] = s->context;
            next_token(s);

            if (s->type != TOK_OPEN) {
                s->type = TOK_ERROR;
            } else {
                int i;
                for(i = 0; i < arity; i++) {
                    next_token(s);
                    ret->parameters[i] = expr(s);
                    if(s->type != TOK_SEP) {
                        break;
                    }
                }
                if(s->type != TOK_CLOSE || i != arity - 1) {
                    s->type = TOK_ERROR;
                } else {
                    next_token(s);
                }
            }

            break;

        case TOK_OPEN:
            next_token(s);
            ret = list(s);
            if (s->type != TOK_CLOSE) {
                s->type = TOK_ERROR;
            } else {
                next_token(s);
            }
            break;

        default:
            ret = new_expr(0, 0);
            s->type = TOK_ERROR;
            ret->binding.value = NAN;
            break;
    }

    return ret;
}


static te_expr *power(state *s) {
    /* <power>     =    {("-" | "+")} <base> */
    int sign = 1;
    while (s->type == TOK_INFIX && (s->binding.function == (funcptr)add || s->binding.function == (funcptr)sub || s->binding.function == (funcptr)lnot || s->binding.function == (funcptr)bnot)) {
        if (s->binding.function == (funcptr)sub) sign = -sign;
        next_token(s);
    }

    te_expr *ret;

    if (sign == 1) {
        ret = base(s);
    } else {
        ret = NEW_EXPR(TE_FUNCTION1 | TE_FLAG_PURE, base(s));
        ret->binding.function = (funcptr)negate;
    }

    return ret;
}

#ifdef TE_POW_FROM_RIGHT
static te_expr *factor(state *s) {
    /* <factor>    =    <power> {"^" <power>} */
    te_expr *ret = power(s);

    int neg = 0;
    te_expr *insertion = 0;

    if (ret->type == (TE_FUNCTION1 | TE_FLAG_PURE) && ret->function == negate) {
        te_expr *se = ret->parameters[0];
        free(ret);
        ret = se;
        neg = 1;
    }

    while (s->type == TOK_INFIX && (s->function == pow)) {
        te_fun2 t = s->function;
        next_token(s);

        if (insertion) {
            /* Make exponentiation go right-to-left. */
            te_expr *insert = NEW_EXPR(TE_FUNCTION2 | TE_FLAG_PURE, insertion->parameters[1], power(s));
            insert->function = t;
            insertion->parameters[1] = insert;
            insertion = insert;
        } else {
            ret = NEW_EXPR(TE_FUNCTION2 | TE_FLAG_PURE, ret, power(s));
            ret->function = t;
            insertion = ret;
        }
    }

    if (neg) {
        ret = NEW_EXPR(TE_FUNCTION1 | TE_FLAG_PURE, ret);
        ret->function = negate;
    }

    return ret;
}
#else
static te_expr *factor(state *s) {
    /* <factor>    =    <power> {"^" <power>} */
    te_expr *ret = power(s);

    while (s->type == TOK_INFIX && (s->binding.function == (funcptr)pow)) {
        te_fun2 t = (te_fun2)(s->binding.function);
        next_token(s);
        ret = NEW_EXPR(TE_FUNCTION2 | TE_FLAG_PURE, ret, power(s));
        ret->binding.function = (funcptr)t;
    }

    return ret;
}
#endif



static te_expr *term(state *s) {
    /* <term>      =    <factor> {("*" | "/" | "%") <factor>} */
    te_expr *ret = factor(s);

    while (s->type == TOK_INFIX && (s->binding.function == (funcptr)mul || s->binding.function == (funcptr)divide || s->binding.function == (funcptr)fmod)) {
        te_fun2 t = (te_fun2)(s->binding.function);
        next_token(s);
        ret = NEW_EXPR(TE_FUNCTION2 | TE_FLAG_PURE, ret, factor(s));
        ret->binding.function = (funcptr)t;
    }

    return ret;
}


static te_expr *expr(state *s) {
    /* <expr>      =    <term> {("+" | "-") <term>} */
    te_expr *ret = term(s);

    while (s->type == TOK_INFIX && (s->binding.function == (funcptr)add ||
            s->binding.function == (funcptr)sub ||
            s->binding.function == (funcptr)lt ||
            s->binding.function == (funcptr)lte ||
            s->binding.function == (funcptr)gt ||
            s->binding.function == (funcptr)gte ||
            s->binding.function == (funcptr)eq ||
            s->binding.function == (funcptr)neq ||
            s->binding.function == (funcptr)lor ||
            s->binding.function == (funcptr)land ||
            s->binding.function == (funcptr)band ||
            s->binding.function == (funcptr)bor ||
            s->binding.function == (funcptr)shr ||
            s->binding.function == (funcptr)shl)) {
        te_fun2 t = (te_fun2)s->binding.function;
        next_token(s);
        ret = NEW_EXPR(TE_FUNCTION2 | TE_FLAG_PURE, ret, term(s));
        ret->binding.function = (funcptr)t;
    }

    return ret;
}


static te_expr *list(state *s) {
    /* <list>      =    <expr> {"," <expr>} */
    te_expr *ret = expr(s);

    while (s->type == TOK_SEP) {
        next_token(s);
        ret = NEW_EXPR(TE_FUNCTION2 | TE_FLAG_PURE, ret, expr(s));
        ret->binding.function = (funcptr)comma;
    }

    return ret;
}


#define TE_FUN(...) ((double(*)(__VA_ARGS__))n->binding.function)
#define M(e) te_eval(n->parameters[e])


double te_eval(const te_expr *n) {
    if (!n) return NAN;

    switch(TYPE_MASK(n->type)) {
        case TE_CONSTANT: return n->binding.value;
        case TE_VARIABLE: return *n->binding.bound;

        case TE_FUNCTION0: case TE_FUNCTION1: case TE_FUNCTION2: case TE_FUNCTION3:
        case TE_FUNCTION4: case TE_FUNCTION5: case TE_FUNCTION6: case TE_FUNCTION7:
            switch(ARITY(n->type)) {
                case 0: return TE_FUN(void)();
                case 1: return TE_FUN(double)(M(0));
                case 2: return TE_FUN(double, double)(M(0), M(1));
                case 3: return TE_FUN(double, double, double)(M(0), M(1), M(2));
                case 4: return TE_FUN(double, double, double, double)(M(0), M(1), M(2), M(3));
                case 5: return TE_FUN(double, double, double, double, double)(M(0), M(1), M(2), M(3), M(4));
                case 6: return TE_FUN(double, double, double, double, double, double)(M(0), M(1), M(2), M(3), M(4), M(5));
                case 7: return TE_FUN(double, double, double, double, double, double, double)(M(0), M(1), M(2), M(3), M(4), M(5), M(6));
                default: return NAN;
            }

        case TE_CLOSURE0: case TE_CLOSURE1: case TE_CLOSURE2: case TE_CLOSURE3:
        case TE_CLOSURE4: case TE_CLOSURE5: case TE_CLOSURE6: case TE_CLOSURE7:
            switch(ARITY(n->type)) {
                case 0: return TE_FUN(void*)(n->parameters[0]);
                case 1: return TE_FUN(void*, double)(n->parameters[1], M(0));
                case 2: return TE_FUN(void*, double, double)(n->parameters[2], M(0), M(1));
                case 3: return TE_FUN(void*, double, double, double)(n->parameters[3], M(0), M(1), M(2));
                case 4: return TE_FUN(void*, double, double, double, double)(n->parameters[4], M(0), M(1), M(2), M(3));
                case 5: return TE_FUN(void*, double, double, double, double, double)(n->parameters[5], M(0), M(1), M(2), M(3), M(4));
                case 6: return TE_FUN(void*, double, double, double, double, double, double)(n->parameters[6], M(0), M(1), M(2), M(3), M(4), M(5));
                case 7: return TE_FUN(void*, double, double, double, double, double, double, double)(n->parameters[7], M(0), M(1), M(2), M(3), M(4), M(5), M(6));
                default: return NAN;
            }

        default: return NAN;
    }

}

#undef TE_FUN
#undef M

static void optimize(te_expr *n) {
    /* Evaluates as much as possible. */
    if (n->type == TE_CONSTANT) return;
    if (n->type == TE_VARIABLE) return;

    /* Only optimize out functions flagged as pure. */
    if (IS_PURE(n->type)) {
        const int arity = ARITY(n->type);
        int known = 1;
        int i;
        for (i = 0; i < arity; ++i) {
            optimize(n->parameters[i]);
            if (((te_expr*)(n->parameters[i]))->type != TE_CONSTANT) {
                known = 0;
            }
        }
        if (known) {
            const double value = te_eval(n);
            te_free_parameters(n);
            n->type = TE_CONSTANT;
            n->binding.value = value;
        }
    }
}


te_expr *te_compile(const char *expression, const te_variable *variables, int var_count, int *error) {
    state s;
    s.start = s.next = expression;
    s.lookup = variables;
    s.lookup_len = var_count;

    next_token(&s);
    te_expr *root = list(&s);

    if (s.type != TOK_END) {
        te_free(root);
        if (error) {
            *error = (s.next - s.start);
            if (*error == 0) *error = 1;
        }
        return 0;
    } else {
        optimize(root);
        if (error) *error = 0;
        return root;
    }
}


double te_interp(const char *expression, int *error) {
    te_expr *n = te_compile(expression, 0, 0, error);
    double ret;
    if (n) {
        ret = te_eval(n);
        te_free(n);
    } else {
        ret = NAN;
    }
    return ret;
}



//static void pn (const te_expr *n, int depth) {
//    int i, arity;
//    printf("%*s", depth, "");
//
//    switch(TYPE_MASK(n->type)) {
//    case TE_CONSTANT: printf("%f\n", n->value); break;
//    case TE_VARIABLE: printf("bound %p\n", n->bound); break;
//
//    case TE_FUNCTION0: case TE_FUNCTION1: case TE_FUNCTION2: case TE_FUNCTION3:
//    case TE_FUNCTION4: case TE_FUNCTION5: case TE_FUNCTION6: case TE_FUNCTION7:
//    case TE_CLOSURE0: case TE_CLOSURE1: case TE_CLOSURE2: case TE_CLOSURE3:
//    case TE_CLOSURE4: case TE_CLOSURE5: case TE_CLOSURE6: case TE_CLOSURE7:
//         arity = ARITY(n->type);
//         printf("f%d", arity);
//         for(i = 0; i < arity; i++) {
//             printf(" %p", n->parameters[i]);
//         }
//         printf("\n");
//         for(i = 0; i < arity; i++) {
//             pn(n->parameters[i], depth + 1);
//         }
//         break;
//    }
//}


//void te_print(const te_expr *n) {
//    pn(n, 0);
//}
