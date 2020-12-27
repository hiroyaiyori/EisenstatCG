//
// Created by 伊従寛哉 on 2020/12/27.
//
#include <stdio.h>
#include <time.h>
#include "gps_opt.h"

void main(){
    // 連立方程式 Ax = b
    // 行列Aは正定値対象行列
    matrix a = { {5.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0},
                 {0.0, 5.0, 2.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0},
                 {0.0, 2.0, 5.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
                 {1.0, 0.0, 0.0, 5.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0},
                 {0.0, 0.0, 0.0, 0.0, 5.0, 0.0, 0.0, 0.0, 0.0, 1.0},
                 {0.0, 0.0, 0.0, 0.0, 0.0, 5.0, 2.0, 0.0, 0.0, 0.0},
                 {0.0, 0.0, 0.0, 1.0, 0.0, 2.0, 5.0, 0.0, 0.0, 1.0},
                 {1.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 5.0, 2.0, 0.0},
                 {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 5.0, 2.0},
                 {0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 1.0, 0.0, 2.0, 5.0} };
    crsdata data;
    crscol col;
    crsrow row;
    int data_i = 0;
    for (int i=0; i<N; i++){
        for (int j=0; j<N; j++){
            if (a[i][j] != 0.0){
                data[data_i] = a[i][j];
                col[data_i] = j;
                data_i++;
            }
        }
        row[i+1] = data_i;
    }

    int v;
    long cpu_time_0, cpu_time_1, cg_cpu_time;
    ivector numbered;
    cpu_time_0 = clock();

    // TRI preconditioning CG法でAx=bを解く
    minimal_digree_all(row, &v);

//    printf("--------------- minimal_degree_all done -------------\n");
//
//    printf("--------------- generate_level_structure start -------------\n");
//    ivector level, x_level;
//    level_index_structure level_index;
//    generate_level_structure(level,  level_index, x_level, col, row, v);
//
//    printf("--------------- x_level -------------\n");
//    print_int_vector(x_level);
//    printf("--------------- generate_level_structure done -------------\n");
//
//
//    printf("--------------- degree_increasing_sort -------------\n");
//    ivector i_arr = {1,5,3,7,8,2,0,0,0,0};
//    int len = 6;
//    degree_increasing_sort(i_arr, len, row);
//    printf("=================== i_arr =========================\n");
//    print_int_vector(i_arr);
//    printf("--------------- degree_increasing_sort end-------------\n");
//
//
//    printf("--------------- pseudo_diameter start-------------\n");
//    level_index_structure level_index_v, level_index_u;
//    ivector x_level_v, x_level_u, level_v, level_u;
//    int u, level_depth;
//    pseudo_diameter(level_index_v, level_index_u, x_level_v, x_level_u, level_v, level_u, row, col, &u, &v, &level_depth);
//    printf("===============x_level_v===========\n");
//    print_int_vector(x_level_v);
//    printf("===============x_level_u=============\n");
//    print_int_vector(x_level_u);
//    printf("--------------- pseudo_diameter done-------------\n");
//
//
//    printf("--------------- minimizing_level_width start-------------\n");
//    int chosed_element;
//    minimizing_level_width(row, col, x_level,level_index, level_index_v, level_index_u,level_v, level_u, &u, &v, &level_depth, &chosed_element);
//    printf("x_level\n");
//    print_int_vector(x_level);
//    printf("--------------- minimizing_level_width done-------------\n");

    printf("--------------- gps start -------------\n");
    gps(row, col, numbered);
    printf("numbered\n");
    print_int_vector(numbered);
    printf("--------------- gps done -------------\n");
    crsdata optdata;
    crscol optcol;
    crsrow optrow = {0};
    opt_mat(numbered, data, optdata, col, optcol, row, optrow);

    printf("&&&       original matrix        &&&\n");
    print_crs_matrix(data, col, row);
    printf("&&&       original matrix        &&&\n");
    printf("&&&       bandwidth minimized matrix        &&&\n");
    print_crs_matrix(optdata, optcol, optrow);
    printf("&&&       bandwidth minimized matrix        &&&\n");

    cpu_time_1 = clock();
    cg_cpu_time = cpu_time_1 - cpu_time_0;
    printf("--------------- calc done -------------\n");
    printf("GPS ALGO CPU time： %ld\n", cg_cpu_time);

    return;
}