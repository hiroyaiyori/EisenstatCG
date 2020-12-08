#include <stdio.h>
#include <math.h>
#include <time.h>

#define N    10
#define M    3
#define NON_ZERO	30	
#define TMAX 100
#define EPS  (1.0e-6)
#define EPS_D (2.0e-6)
#define OMEGA  1.5

// CRS方式
typedef double crsdata[NON_ZERO];
typedef int crscol[NON_ZERO];
typedef int crsrow[N+1];

typedef double matrix[N][N];    // 行列
typedef double vector[N];       // ベクトル
typedef int ivector[N];



int print_vector(vector x)
{
    int i;
    for(i = 0; i < N; i++){
        printf("vec[%d] = %2g\n", i, x[i]);
    }
    return 0;
}

int print_int_vector(ivector x)
{
    int i;
    for(i = 0; i < 10; i++){
        printf("vec[%d] = %d\n", i, x[i]);
    }
    return 0;
}

void copy_vector(vector y, const vector x)
{
    int i;
    for (i = 0; i < N; i++)
        y[i] = x[i];
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

void tri_cg(const int maxt, crsdata data, crsrow row, crscol col, const vector b, vector x, const double eps,
        const double eps_d, const double omega, const int m, int *itern) {
    int i, j, k;
    double eps_tilde_rxr0, eps_rxr0, pAp, alpha, beta, rxr, pre_rxr;
    static vector r;
    static vector p;
    static vector q;
    static vector Ap;
    static vector D;
    static ivector Di;
    static vector s;
    static vector rt;
    double omega_i = 1 - 2 / omega;

/* r := A x */
    multiply_mv(r, data, row, col, x);

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
/* EISENSTAT TRICK  */

/* s := (L + D/ω)^-1 * p  */
        copy_vector(s, p);
        forward_substitution(data, row, col, s, Di);
/* q := p + (1 - 2/ω)D(L + D/ω)^-1 * p  */
        for (i = 0; i < N; i++){
            q[i] = p[i] + omega_i * D[i] * s[i];
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
int main(void) {
    // 連立方程式 Ax = b
    // 行列Aは正定値対象行列
    matrix a = { {5.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
                 {2.0, 5.0, 2.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0},
                 {0.0, 2.0, 5.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
                 {0.0, 0.0, 2.0, 5.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0},
                 {0.0, 0.0, 0.0, 2.0, 5.0, 2.0, 0.0, 0.0, 0.0, 0.0},
                 {0.0, 0.0, 0.0, 0.0, 2.0, 5.0, 2.0, 0.0, 0.0, 0.0},
                 {0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 5.0, -2.0, 0.0, 0.0},
                 {0.0, 2.0, 0.0, 0.0, 0.0, 0.0, -2.0, 5.0, 2.0, 0.0},
                 {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 5.0, 2.0},
                 {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 5.0} };
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

    vector b = {3.0, 1.0, 4.0, 0.0, 5.0, -1.0, 6.0, -2.0, 7.0, -15.0};
    // 初期値を設定
    vector x = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};

    int    i;
    int itern;
	long cpu_time_0, cpu_time_1, cg_cpu_time;
	cpu_time_0 = clock();

    // TRI preconditioning CG法でAx=bを解く
    tri_cg(TMAX, data, row, col, b, x, EPS, EPS_D, OMEGA, M, &itern);

	cpu_time_1 = clock();
	cg_cpu_time = cpu_time_1 - cpu_time_0;

    vector ans_x = {0.637207, -0.093018, 1.42197, -1.4619, 2.23279, -1.62006, 1.31737, -1.32663, 3.72697, -4.49079};

    printf("--------------- calc done -------------\n");
	printf("CG Method CPU time： %ld\n", cg_cpu_time);
    double err_sum = 0.0;
    for(i = 0; i < N; i++){
        err_sum += (x[i] - ans_x[i]) * (x[i] - ans_x[i]);
        printf("x[%d] = %2g\n", i, x[i]);
    }
    printf("error %2g\n", err_sum);
    return 0;
}
