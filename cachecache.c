//
// Created by 伊従寛哉 on 2020/12/28.
//

#include <stdio.h>
#include <math.h>
#include <time.h>
#include "matrix_setting.h"
#include "cachecache.h"
#include "print_util.h"



void copy_vector(vector y, const vector x)
{
    int i;
    for (i = 0; i < N; i++)
        y[i] = x[i];
}

void lhat_res_split(crsdata data, crsdata resdata, crscol col, crscol rescol, crsrow row, crsrow resrow){
    int i, j, k, blcsep, rtmprown, htmprown;
    int spl_p[PARALLEL_N + 1] = {0};
    int rnp = (N - 1) % PARALLEL_N;
    int quo = (N - 1) / PARALLEL_N;
    int hn = 0;
    int resn = 0;
    blcsep = 0;
    htmprown = 0;
    rtmprown = 0;
    for (i = 0; i < PARALLEL_N; i++){
        if (i < rnp)
            spl_p[i + 1] = spl_p[i] + quo + 1;
        else
            spl_p[i + 1] = spl_p[i] + quo;
    }
    if (spl_p[PARALLEL_N] != N - 1){
        printf("\n++++++++++++++++++++++++++++++++ERROR!!++++++++++++++++++++++++++%d\n", spl_p[PARALLEL_N]);
    }
    for (i = 0; i < PARALLEL_N; i++){
        j =  spl_p[i];
        for (k = blcsep; k < row[j + 1]; k++){
            if (col[k] <  1 + spl_p[i + 1]){
                data[hn] = data[k];
                col[hn] = col[k];
                hn++;
            }else {
                resdata[resn] = data[k];
                rescol[resn] = col[k];
                resn++;
            }
        }
        row[j] = htmprown;
        resrow[j] = rtmprown;
        htmprown = hn;
        rtmprown = resn;
        for (j = spl_p[i] + 1; j < spl_p[i + 1]; j++){
            for (k = row[j]; k < row[j + 1]; k++){
                if (col[k] <  spl_p[i]){
                    resdata[resn] = data[k];
                    rescol[resn] = col[k];
                    resn ++;
                }else if (col[k] > spl_p[i + 1]){
                    resdata[resn] = data[k];
                    rescol[resn] = col[k];
                    resn++;
                }else {
                    data[hn] = data[k];
                    col[hn] = col[k];
                    hn++;
                }
            }
            row[j] = htmprown;
            resrow[j] = rtmprown;
            htmprown = hn;
            rtmprown = resn;
        }
        j = spl_p[i + 1];
        for (k = row[j]; k < row[j + 1]; k++){
            if (col[k] <  spl_p[i]){
                resdata[resn] = data[k];
                rescol[resn] = col[k];
                resn++;
            }else {
                data[hn] = data[k];
                col[hn] = col[k];
                hn++;
                if (col[k] == j){
                    blcsep = k + 1;
                    break;
                }
            }
        }
    }
    row[N - 1] = htmprown;
    resrow[N - 1] = rtmprown;
    row[N] = hn;
    resrow[N] = resn;

    return;
}

double crs_dotproduct(const crsdata data, const crscol col, const int si, const int ei, const vector y)
{
    int i;
    double sum;
    sum = 0;
    for (i = si; i < ei; i++)
        sum += data[i] * y[col[i]];
    return sum;
}

double dotproduct(const vector x, const vector y) {
    int i;
    double sum;
    sum = 0;
    for (i = 0; i < N; i++)
        sum += x[i] * y[i];
    return sum;
}

double norm(vector x)
{
    return sqrt(dotproduct(x, x));
}

void multiply_mv(vector Ab, crsdata data, crsrow row, crscol col, vector b)
{
    int i;
    for (i = 0; i < N; i++){
        Ab[i] = crs_dotproduct(data, col, row[i], row[i+1], b);
    }
}

void multiply_sv(vector ab, double a, const vector b)
{
    int i;
    for (i = 0; i < N; i++){
        ab[i] = a * b[i];
    }
}

void forward_substitution(const crsdata L, const crsrow row, const crscol col, vector x, const ivector Di){
    int i, j;
    for(i = 0; i < N; i++){
        for(j = row[i]; j < Di[i]; j++){
            x[i] -= L[j] * x[col[j]];
        }
        x[i] /= L[Di[i]];
    }
}

void backward_substitution(const crsdata U, const crsrow row, const crscol col, vector x, const ivector Di){
    int i, j;
    for(i = N - 1; i >= 0; i--){
        for(j = row[i+1] - 1; j > Di[i]; j--){
            x[i] -= U[j] * x[col[j]];
        }
        x[i] /= U[Di[i]];
    }
}

void cache_cache(const int maxt, crsdata data, crsrow row, crscol col, const vector b, vector x, const double eps,
            const double eps_d, const double omega, const int m, int *itern) {
    int i, j, k;
    crsdata resdata = {};
    crscol rescol = {};
    crsrow resrow = {};
    double eps_tilde_rxr0, eps_rxr0, pAp, alpha, beta, rxr, pre_rxr;
    vector r = {};
    vector p = {};
    vector q = {};
    vector Ap = {};
    vector D = {};
    ivector Di;
    vector rs = {};
    vector s = {};
    vector rt = {};
    double omega_i = 1 - 2 / omega;

/* r := A x */
    multiply_mv(r, data, row, col, x);
/* L := L^hat*/
    lhat_res_split(data, resdata, col, rescol, row, resrow);

/* r := b - A x */
/* A := L + D/ω + L^T*/
    for (i = 0; i < N; i++){
        r[i] = b[i] - r[i];
        for (j = row[i]; j < row[i+1]; j++){
            if (col[j] == i){
                D[i] = data[j];
                Di[i] = j;
                data[j] /= omega;
                break;
            }
        }
    }
    /* eps_d_rxr0 := (εd * (r,r)) */
    eps_rxr0 = eps_d * dotproduct(r,r);

/* r_tilde := (L^T + D/ω)^-1 * r */
    backward_substitution(data, row, col, r, Di);
/* rxr := (r, r) */
    rxr = dotproduct(r,r);

    /* eps_tilde_rxr0 := (ε * (r_tilde,r_tilde)) */
    eps_tilde_rxr0 = eps * rxr;


/* p = r */
    copy_vector(p, r);
/* 反復 */
    for (k = 0; k < maxt; k++) {
        pre_rxr = rxr;
/* CACHE-CACHE ELEMENTS  */

/* s := (L + D/ω)^-1 * p  */
        copy_vector(s, p);
        forward_substitution(data, row, col, s, Di);
/* q := p + (1 - 2/ω)D(L + D/ω)^-1 * p  */
        for (i = 0; i < N; i++){
            q[i] = p[i] + omega_i * D[i] * s[i];
        }
/* q := q + R * s */
        multiply_mv(rs, resdata, resrow, rescol, s);
        for (i = 0; i < N; i++){
            q[i] += rs[i];
        }
/* q := (L^T + D/ω)^-1 * q  */
        backward_substitution(data, row, col, q, Di);
/* Ap := s + q  */
        for (i = 0; i < N; i++){
            Ap[i] = s[i] + q[i];
        }

/* EISENSTAT TRICK DONE */

/* pAp := (p, Ap) */
        pAp = dotproduct(p, Ap);
/* α = (r,r)/(p,Ap) */
        alpha = pre_rxr / pAp;
/* x, r の更新 */
        for (i = 0; i < N; i++) {
/* x = x + α*p */
            x[i] += alpha * s[i];
/* r = r - α*A*p */
            r[i] -= alpha * Ap[i];
        }
/* calc new rxr */
        rxr = dotproduct(r,r);


        printf("LOOP : %d\t Error = %g\n", k, rxr);
/* if (rt, rt)/(r0, r0) < eps then stop */
        if ((rxr < eps_tilde_rxr0) && (k % m == 0)){
            /* rt = (L^T + D/ω) * r */
            for (i = 0; i < N; i++){
                rt[i] = 0;
                for (j = Di[i]; j < row[i+1]; j++){
                    rt[i] += data[j] * r[col[j]];
                }
            }
            if (dotproduct(rt, rt) < eps_rxr0)
                break;
        }

        /* β := rxr / pre_rxr */
        beta = rxr / pre_rxr;
        /* p := r + β*p */
        for (i = 0; i < N; i++){
            p[i] = beta * p[i] + r[i];
        }
    }
    *itern = k;
}
