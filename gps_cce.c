//
// Created by 伊従寛哉 on 2020/12/29.
//
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "mmio.h"
#include "gps_opt.h"
#include "cachecache.h"

int main(int argc, char *argv[]){
    int ret_code;
    MM_typecode matcode;
    FILE *f;
    int m, n, nz;
    int i, j, k, l, c;
    double init_val;
    int coltmp, rowtmp;
    double valtmp;
    crscol col, r_row;
    crsdata val;
    crsrow row = {0};

    if (argc < 2)
    {
        fprintf(stderr, "Usage: %s [martix-market-filename]\n", argv[0]);
        exit(1);
    }
    else
    {
        if ((f = fopen(argv[1], "r")) == NULL){
            printf("file open failed");

            exit(1);
        }
    }

    if (mm_read_banner(f, &matcode) != 0)
    {
        printf("Could not process Matrix Market banner.\n");
        exit(1);
    }


    /*  This is how one can screen matrix types if their application */
    /*  only supports a subset of the Matrix Market data types.      */

    if (mm_is_complex(matcode) && mm_is_matrix(matcode) &&
        mm_is_sparse(matcode) )
    {
        printf("Sorry, this application does not support ");
        printf("Market Market type: [%s]\n", mm_typecode_to_str(matcode));
        exit(1);
    }

    /* find out size of sparse matrix .... */

    if ((ret_code = mm_read_mtx_crd_size(f, &m, &n, &nz)) !=0){
        printf("mm_read_mtx_crd_size failed");
        exit(1);
    }

    if((m != N) || (n != N) || (2 * nz - N != NON_ZERO)){
        printf("m = %d, n = %d, nz = %d. Invalid value", m, n, nz);
        exit(1);
    }

    /* NOTE: when reading in doubles, ANSI C requires the use of the "l"  */
    /*   specifier as in "%lg", "%lf", "%le", otherwise errors will occur */
    /*  (ANSI C X3.159-1989, Sec. 4.9.6.2, p. 136 lines 13-15)            */
    j = 1;
    l = 0;
    for (i=0; i<nz; i++)
    {
//        fscanf(f, "%d %d\n", &col[l], &r_row[l]);
//        val[l] = 1;
        fscanf(f, "%d %d %lg\n", &col[l], &r_row[l], &val[l]);

        if( i == 0)
            init_val = val[0];
        val[l] /= init_val;
        col[l]--;  /* adjust from 1-based to 0-based */
        r_row[l]--;
        if ((l != 0) && (r_row[l] != r_row[l - 1])){
            if (r_row[l] == r_row[l - 1] + 1){
                coltmp = col[l];
                rowtmp = r_row[l];
                valtmp = val[l];
                row[j] = l;
                for (k = 0; k < j; k++){
                    for (c = row[k]; c < row[k + 1]; c++){
                        if (col[c] == j){
                            col[l] = k;
                            r_row[l] = j;
                            val[l] = val[c];
                            l++;
                        }
                    }
                }
                col[l] = coltmp;
                r_row[l] = rowtmp;
                val[l] = valtmp;
                l++;
                j++;

            }else{
                printf("matrix is not in correct order%d, %d", r_row[l], r_row[l - 1]);
                exit(1);
            }
        }else{
            l++;
        }
    }


    if (l != NON_ZERO){
        printf("something wrong with making csr");
        exit(1);
    }

    row[j] = NON_ZERO;
    if (j != N){
        printf("the number of row is invalid\n");
    }

    if (f !=stdin) fclose(f);

    /************************/
    /* now write out matrix */
    /************************/


//    for (i=0; i<NON_ZERO; i++)
//        fprintf(stdout, "%d %d %2.0lf\n", col[i]+1, r_row[i]+1, val[i]);


    /************************/
    /* now write out matrix on original format*/
    /************************/
//    fprintf(stdout, "@@@@@@@@@@@@@@@ original matrix @@@@@@@@@@@@@@@@\n\n\n");
//    print_matrix_image(col, row);
    long gps_cpu_0, gps_cpu_1, gps_cpu;
    gps_cpu_0 = clock();
    ivector numbered;
    gps(row, col, numbered);
    crsdata optdata;
    crscol optcol;
    crsrow optrow;
    opt_mat(numbered, val, optdata, col, optcol, row, optrow);
    gps_cpu_1 = clock();
    gps_cpu = gps_cpu_1 - gps_cpu_0;
//    double md=0;
//    for (i=0; i<NON_ZERO; i++){
//        if (val[i]>md)
//            md = val[i];
//    }
//    fprintf(stdout, "\n\n\n\n\n---------@@@@@@@@@@@@@@@@@ minimized matrix @@@@@@@@@@@@@@@@@-----------\n\n\n\n\n");
//    print_matrix_image(optcol, optrow);



//// PRINT OPTED MATRIX
//
//    mm_write_banner(stdout, matcode);
//    mm_write_mtx_crd_size(stdout, N, N, NON_ZERO);
//    int J;
//    J = 0;
//    printf("%d\n",row[N]);
//    for (i=0; i<NON_ZERO; i++){
//        if (optrow[J+1] <= i)
//            J++;
//        fprintf(stdout, "%d %d %20.19g\n", optcol[i]+1, J+1, optdata[i]);
//    }
//
//// PRINT OPTED MATRIX





    int mb, mbopt;
    int mbs, mbopts;
    mb = 0;
    mbopt = 0;
    mbs = 0;
    mbopts = 0;
    for (i = 0; i < N; i++){
        if (mb < col[row[i + 1] - 1] - i){
            mb = col[row[i + 1] - 1] - i;
        }
        if (mbopt < optcol[optrow[i + 1] - 1] - i){
            mbopt = optcol[optrow[i + 1] - 1] - i;
        }
        mbs +=  col[row[i + 1] - 1] - i;
        mbopts += optcol[optrow[i + 1] - 1] - i;

    }
    printf("N:%d,P:%d,omega:%2g\n",N,PARALLEL_N,OMEGA);
    printf("bandwidth decreased %d to %d\n", mb, mbopt);
    printf("bandwidth sum decreased %d to %d\n", mbs, mbopts);

//    print_matrix_image(col, row);
//    printf("\n$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n");
//    print_matrix_image(optcol, optrow);
//    printf("\n$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n");

    vector x = {};
    vector opt_x = {};
//    for (i = 0; i < N; i++){
//        x[i] = 1.0;
//        opt_x[i] = 1.0;
//    }
    vector b = {};
    for (i = 0; i < N; i++){
        for (j = row[i]; j < row[i + 1]; j++){
            b[i] += val[j];
        }
    }
    vector optb = {};
    for (i = 0; i < N; i++){
        for (j = optrow[i]; j < optrow[i + 1]; j++){
            optb[i] += optdata[j];
        }
    }
    int itern, opt_itern;
    long cpu_time_0, cpu_time_1, cpu_time_2;
    cpu_time_0 = clock();
    // TRI preconditioning CG法でAx=bを解く
    cache_cache(TMAX, val, row, col, b, x, EPS, EPS_D, OMEGA, M, &itern);
    cpu_time_1 = clock();
    cache_cache(TMAX, optdata, optrow, optcol, optb, opt_x, EPS, EPS_D, OMEGA, M, &opt_itern);
    cpu_time_2 = clock();
    printf("original cachecache LOOP N: %d\n", itern);
    printf("GPS cachecache LOOP N: %d\n", opt_itern);


//    printf("--------------- calc done -------------\n");
    printf("ORIGINAL MATRIX CPU time： %2g\n", (double)(cpu_time_1 - cpu_time_0)/CLOCKS_PER_SEC);
    printf("BW MATRIX CPU time： %2g, GPS CPU time： %2g\n", (double)(cpu_time_2 - cpu_time_1)/CLOCKS_PER_SEC, (double)gps_cpu/CLOCKS_PER_SEC);
//    for (i = 0; i < N; i++){
//        printf("%2g, ",x[i]);
//    }
    double err_sum, opt_err_sum;
    err_sum = 0;
    opt_err_sum = 0;
    for(i = 0; i < N; i++){
        err_sum += (opt_x[i] - 1) * (opt_x[i] - 1);
        opt_err_sum += (x[i] - 1) * (x[i] - 1);
    }
    printf("error %2g\n", err_sum/N);
    printf("opt_error %2g\n", opt_err_sum/N);
    return 0;
}