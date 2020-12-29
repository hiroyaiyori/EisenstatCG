//
// Created by 伊従寛哉 on 2020/12/29.
//

#ifndef CG_METHOD_PRINT_UTIL_H
#define CG_METHOD_PRINT_UTIL_H
#include "matrix_setting.h"


void print_double_vector(double *x);

void print_int_vector(ivector x);


void print_matrix_image(crscol col, crsrow row);

void print_crs_matrix(crsdata data, crscol col, crsrow row);

#endif //CG_METHOD_PRINT_UTIL_H