//
// Created by 伊従寛哉 on 2020/12/27.
//

#ifndef CG_METHOD_GPS_OPT_H
#define CG_METHOD_GPS_OPT_H

#define N 48
#define NON_ZERO	224 * 2 - N
#define TMAX 100
#define EPS  (1.0e-6)
#define EPS_D (2.0e-6)
#define OMEGA  1.5

// CRS方式
typedef double crsdata[NON_ZERO];
typedef int crscol[NON_ZERO];
typedef int crsrow[N+1];
typedef int level_index_structure[N+1];

typedef double matrix[N][N];    // matrix
typedef double vector[N];       // vector
typedef int ivector[N];

int print_int_vector(ivector x);

void print_double_vector(double *x);

void print_crs_matrix(crsdata data, crscol col, crsrow row);

void print_matrix_image(crscol col, crsrow row);

void minimal_digree_all(crsrow row, int *mdi);

int minimal_digree(crsrow row, ivector level, const int si, const int ei);

int level_width(level_index_structure level_index, int level_depth);

int generate_level_structure(ivector level, level_index_structure level_index, ivector x_level, crscol col, crsrow row, int v);

void degree_increasing_sort(ivector i_arr, const int len, crsrow row);

void pseudo_diameter(int *level_index_v, int *level_index_u,
                     int *x_level_v, int *x_level_u, int *level_v, int *level_u, crsrow row, crscol col, int *u, int *v, int *level_depth);

void minimizing_level_width(crsrow row, crscol col, ivector x_level, level_index_structure level_index, level_index_structure level_index_v,
                            level_index_structure level_index_u, ivector level_v, ivector level_u, int *u, int *v, int *level_depth, int *chosed_element);

void gps(crsrow row, crscol col, ivector numbered);

void opt_mat(ivector numbered, crsdata data, crsdata optdata, crscol col, crscol optcol, crsrow row, crsrow optrow);

#endif //CG_METHOD_GPS_OPT_H
