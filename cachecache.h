//
// Created by 伊従寛哉 on 2020/12/29.
//

#ifndef CG_METHOD_CACHECACHE_H
#define CG_METHOD_CACHECACHE_H

#include "matrix_setting.h"

// CRS方式
typedef double crsdata[NON_ZERO];
typedef int crscol[NON_ZERO];
typedef int crsrow[N+1];

typedef double matrix[N][N];    // 行列
typedef double vector[N];       // ベクトル
typedef int ivector[N];

#define M    3
#define TMAX 100
#define EPS  (1.0e-6)
#define EPS_D (2.0e-6)
#define OMEGA  1.5


int print_vector(vector x);


int print_int_vector(ivector x);


void print_matrix_image(crscol col, crsrow row);


void copy_vector(vector y, const vector x);


void lhat_res_split(crsdata data, crsdata resdata, crscol col, crscol rescol, crsrow row, crsrow resrow);


double crs_dotproduct(const crsdata data, const crscol col, const int si, const int ei, const vector y);


double dotproduct(const vector x, const vector y);


double norm(vector x);

void multiply_mv(vector Ab, crsdata data, crsrow row, crscol col, vector b);


void multiply_sv(vector ab, double a, const vector b);

void forward_substitution(const crsdata L, const crsrow row, const crscol col, vector x, const ivector Di);

void backward_substitution(const crsdata U, const crsrow row, const crscol col, vector x, const ivector Di);

void cache_cache(const int maxt, crsdata data, crsrow row, crscol col, const vector b, vector x, const double eps,
                 const double eps_d, const double omega, const int m, int *itern);
#endif //CG_METHOD_CACHECACHE_H