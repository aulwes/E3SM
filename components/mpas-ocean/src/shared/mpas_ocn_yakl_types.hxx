#pragma once

#include "YAKL.h"

typedef double real;

typedef yakl::Array<real,1,yakl::memHost,yakl::styleFortran>  h_double_1d_t;
typedef yakl::Array<real,1,yakl::memDefault,yakl::styleFortran>  d_double_1d_t;
typedef yakl::Array<int,1,yakl::memHost,yakl::styleFortran>  h_int_1d_t;
typedef yakl::Array<int,1,yakl::memDefault,yakl::styleFortran>  d_int_1d_t;

typedef yakl::Array<real,2,yakl::memHost,yakl::styleFortran>  h_double_2d_t;
typedef yakl::Array<real,2,yakl::memDefault,yakl::styleFortran>  d_double_2d_t;
typedef yakl::Array<int,2,yakl::memDefault,yakl::styleFortran>  d_int_2d_t;
typedef yakl::Array<int,2,yakl::memHost,yakl::styleFortran>  h_int_2d_t;

typedef yakl::Array<real,3,yakl::memHost,yakl::styleFortran>  h_double_3d_t;
typedef yakl::Array<real,3,yakl::memDefault,yakl::styleFortran>  d_double_3d_t;

