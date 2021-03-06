//
// Created by 伊従寛哉 on 2020/12/28.
//

#ifndef CG_METHOD_MATRIX_SETTING_H
#define CG_METHOD_MATRIX_SETTING_H
/*
 * Define here the length and the element number of Non Zeros.
 */
#define N 2910
#define NON_ZERO 88603*2-N
#define PARALLEL_N  8


// CRS方式
typedef double crsdata[NON_ZERO];
typedef int crscol[NON_ZERO];
typedef int crsrow[N+1];

typedef double matrix[N][N];    // 行列
typedef double vector[N];       // ベクトル
typedef int ivector[N];

#endif //CG_METHOD_MATRIX_SETTING_H
