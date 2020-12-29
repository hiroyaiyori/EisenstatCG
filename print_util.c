//
// Created by 伊従寛哉 on 2020/12/29.
//
#include <stdio.h>
#include "matrix_setting.h"


void print_double_vector(double *x){
    int i;
    for(i = 0; i < 10; i++){
        printf("vec[%d] = %3.1lf\n", i, x[i]);
    }
    return;
}

void print_int_vector(ivector x)
{
    int i;
    for(i = 0; i < 10; i++){
        printf("vec[%d] = %d\n", i, x[i]);
    }
    return;
}

void print_matrix_image(crscol col, crsrow row) {
    int i, j, k, flg;
    flg = 0;
    for (i = 0; i < N; i++) {
        for (j = 0; j < N; j++) {
            for (k = row[i]; k < row[i + 1]; k++) {
                if (col[k] == j) {
                    printf("◉ ");
                    flg = 1;
                    break;
                }
            }
            if (!flg) {
                printf("◯ ");
            }
            flg = 0;
        }
        printf("\n");
    }
}


void print_crs_matrix(crsdata data, crscol col, crsrow row){
    int i, j, k, flg;
    flg = 0;
    for (i = 0; i < N; i++){
        for (j = 0; j < N; j++){
            for (k = row[i]; k < row[i + 1]; k++){
                if (col[k] == j){
                    fprintf(stdout, "%3.1lf\t", data[k]);
                    flg = 1;
                    break;
                }
            }
            if (!flg){
                fprintf(stdout, "%3.1lf\t", 0.0);
            }
            flg = 0;
        }
        printf("\n");
    }
}

