#include <stdio.h>
#include <math.h>

#define N    10
#define TMAX 100
#define EPS  (1.0e-6)
#define OMEGA  1.5

typedef double vector[N];       // ベクトル
typedef double matrix[N][N];    // 行列



int print_vector(vector x)
{
    int i;
    for(i = 0; i < N; i++){
        printf("vec[%d] = %2g\n", i, x[i]);
    }
    return 0;
}

void copy_vector(vector y, const vector x)
{
    int i;
    for (i = 0; i < N; i++)
        y[i] = x[i];
}

double dotproduct(const vector x, const vector y)
{
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

void multiply_mv(vector Ab, matrix A, vector b)
{
    int i;
    for (i = 0; i < N; i++){
        Ab[i] = dotproduct(A[i], b);
    }
}

void multiply_sv(vector ab, double a, const vector b)
{
   int i;
   for (i = 0; i < N; i++){
       ab[i] = a * b[i];
   }
}

void forward_substitution(matrix L, vector x){
    int i, j;
    for(i = 0; i < N; i++){
        for(j = 0; j < i; j++){
            x[i] -= L[j][i] * x[j];
        }
        x[i] = x[i] / L[i][i];
    }
}

void backward_substitution(matrix U, vector x){
    int i, j;
    for(i = N - 1; i >= 0; i--){
        for(j = i + 1; j < N; j++){
            x[i] = x[i] - U[i][j] * x[j];
        }
        x[i] = x[i] / U[i][i];
    }
}

void eisenstat_cg(const int maxt, matrix A, const vector b, vector x, const double eps, const double omega, int *itern) {
    int i, j, k;
    double eps_rxr0, pAp, alpha, beta, rxr, pre_rxr;
    static vector r;
    static vector p;
    static vector q;
    static vector Ap;
    static vector D;
    static vector s;
    static vector rt;
    double omega_i = 1 - 2 / omega;

/* r := A x */
    multiply_mv(r, A, x);

/* r := b - A x */
/* A := L + D/ω + L^T*/
    for (i = 0; i < N; i++){
        r[i] = b[i] - r[i];
        D[i] = A[i][i];
        A[i][i] /= omega;
    }
    /* eps_b2 := (ε ||b||)^2 */
    eps_rxr0 = eps * dotproduct(r,r);


/* r_tilde := (L^T + D/ω)^-1 * r */
    backward_substitution(A, r);
/* rxr := (r, r) */
    rxr = dotproduct(r,r);


/* p = r */
    copy_vector(p, r);
/* 反復 */
    for (k = 0; k < maxt; k++) {
        pre_rxr = rxr;
/* EISENSTAT TRICK  */

/* s := (L + D/ω)^-1 * p  */
        copy_vector(s, p);
        forward_substitution(A, s);
/* q := p + (1 - 2/ω)D(L + D/ω)^-1 * p  */
        for (i = 0; i < N; i++){
            q[i] = p[i] + omega_i * D[i] * s[i];
        }
/* q := (L^T + D/ω)^-1 * q  */
        backward_substitution(A, q);
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

/* rt = (L^T + D/ω) * r */
        for (i = 0; i < N; i++){
            rt[i] = 0;
            for (j = i; j < N; j++){
                rt[i] += A[i][j] * r[j];
            }
        }
        printf("LOOP : %d\t Error = %g\n", k, rxr);
/* if (rt, rt)/(r0, r0) < eps then stop */
        if (dotproduct(rt, rt) < eps_rxr0)
            break;
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
    vector b = {3.0, 1.0, 4.0, 0.0, 5.0, -1.0, 6.0, -2.0, 7.0, -15.0};
    // 初期値を設定
    vector x = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};

    int    i;
    int itern;

    // EISENSTAT CG法でAx=bを解く
    eisenstat_cg(TMAX, a, b, x, EPS, OMEGA, &itern);

    vector ans_x = {0.637207, -0.093018, 1.42197, -1.4619, 2.23279, -1.62006, 1.31737, -1.32663, 3.72697, -4.49079};

    printf("--------------- calc done -------------\n");
    double err_sum = 0.0;
    for(i = 0; i < N; i++){
        err_sum += (x[i] - ans_x[i]) * (x[i] - ans_x[i]);
        printf("x[%d] = %2g\n", i, x[i]);
    }
    printf("error %2g\n", err_sum);
    return 0;
}