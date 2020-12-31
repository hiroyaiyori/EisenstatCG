//
// Created by 伊従寛哉 on 2020/12/29.
//
#include <stdio.h>
#include <time.h>

#include "matrix_setting.h"
#include "cachecache.h"

int main() {
    int i,j;
    // 連立方程式 Ax = b
    // 行列Aは正定値対象行列
    matrix a = {{5.0, 2.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0},
                 {2.0, 5.0, 2.0, 1.0, 1.0, 1.0, 1.0, 2.0, 1.0, 1.0},
                 {1.0, 2.0, 5.0, 2.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0},
                 {1.0, 1.0, 2.0, 5.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0},
                 {1.0, 1.0, 1.0, 2.0, 5.0, 2.0, 1.0, 1.0, 0.0, 0.0},
                 {1.0, 1.0, 1.0, 0.0, 2.0, 5.0, 2.0, 1.0, 0.0, 0.0},
                 {1.0, 1.0, 1.0, 0.0, 1.0, 2.0, 5.0, -2.0, 0.0, 0.0},
                 {1.0, 2.0, 1.0, 0.0, 1.0, 1.0, -2.0, 5.0, 2.0, 0.0},
                 {1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 2.0, 5.0, 2.0},
                 {1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 5.0}};

//    matrix a = { {5.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
//                 {2.0, 5.0, 2.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0},
//                 {0.0, 2.0, 5.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
//                 {0.0, 0.0, 2.0, 5.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0},
//                 {0.0, 0.0, 0.0, 2.0, 5.0, 2.0, 0.0, 0.0, 0.0, 0.0},
//                 {0.0, 0.0, 0.0, 0.0, 2.0, 5.0, 2.0, 0.0, 0.0, 0.0},
//                 {0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 5.0, -2.0, 0.0, 0.0},
//                 {0.0, 2.0, 0.0, 0.0, 0.0, 0.0, -2.0, 5.0, 2.0, 0.0},
//                 {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 5.0, 2.0},
//                 {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 5.0} };

    int nz = 0;
    for (int i = 0; i < N; i++){
        for (int j = 0; j < N; j++){
            if (a[i][j] != 0.0)
                nz++;
        }
    }
    printf("nz = %d\n", nz);

    crsdata data, resdata;
    crscol col, rescol;
    crsrow row, resrow;
    row[0] = 0;
    resrow[0] = 0;
    int data_i = 0;
    for (i=0; i<N; i++){
        for (j=0; j<N; j++){
            if (a[i][j] != 0.0){
                data[data_i] = a[i][j];
                col[data_i] = j;
                data_i++;
            }
        }
        row[i+1] = data_i;
    }
//    printf("\noriginal matrix\n");
//    print_matrix_image(data, col, row);
//    lhat_res_split(data, resdata, col, rescol, row, resrow);
//    printf("\nhat matrix\n");
//    print_matrix_image(col, row);
//    printf("\nresidual matrix\n");
//    print_matrix_image(rescol, resrow);



    vector b = {3.0, 1.0, 4.0, 0.0, 5.0, -1.0, 6.0, -2.0, 7.0, -15.0};
    // 初期値を設定
    vector x = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};

    int itern;
    long cpu_time_0, cpu_time_1, cg_cpu_time;
    cpu_time_0 = clock();

    // TRI preconditioning CG法でAx=bを解く
    cache_cache(TMAX, data, row, col, b, x, EPS, EPS_D, OMEGA, M, &itern);

    cpu_time_1 = clock();
    cg_cpu_time = cpu_time_1 - cpu_time_0;

    printf("--------------- calc done -------------\n");
    printf("CG Method CPU time： %ld\n", cg_cpu_time);
    vector ans_x = {1.22896, 1.43526, 2.06826, -2.3582, 2.49526, -0.427202, -1.94958, -4.6877, 4.65143, -5.80701};
    for (i = 0; i < N; i++){
        printf("%2g, ",x[i]);
    }
    printf("\n");
    double err_sum;
    err_sum = 0;
    for(i = 0; i < N; i++){
        err_sum += (x[i] - ans_x[i]) * (x[i] - ans_x[i]);
        printf("x[%d] = %2g\n", i, x[i]);
    }
    printf("error %2g\n", err_sum);
    return 0;
}

