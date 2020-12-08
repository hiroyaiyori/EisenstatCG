#include <stdio.h>
#include <math.h>
#include <time.h>

#define N    10
#define TMAX 100
#define EPS  (1.0e-6)

typedef double vector[N];       // ベクトル
typedef double matrix[N][N];    // 行列



int print_vector(vector x)
{
    int i;
    for(i = 0; i < N; i++){
        printf("x[%d] = %2g\n", i, x[i]);
    }
    return 0;
}
//void print_matrix(int m, int n, matrix a)
//{
//    int i, j;
//    for (i = 0; i < m; i++) {
//        for (j = 0; j < n; j++) {
//            cout << a[i][j];
//            if (j % 5 == 4)
//                cout << endl;
//        }
//        if (j % 5 != 0)
//            cout << endl;
//    }
//}
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

void copy_vector(int n, vector y, const vector x)
{
    int i;
    for (i = 0; i < n; i++)
        y[i] = x[i];
}

void cg(int maxt, matrix A, vector b, vector x, double eps, int *maxiter) {
    int i, k;
    double eps_b, eps_b2, pAp, alpha, beta, rxr;
    static vector r, p, Ap;
/* eps_b2 := (ε ||b||)^2 */
    eps_b = eps * norm(b);
    eps_b2 = eps_b * eps_b;
/* r := A x */
    multiply_mv(r, A, x);

/* r := b - A x */
    for (i = 0; i < N; i++)
        r[i] = b[i] - r[i];
/* β := 1 / (r,r); p := β r */
    beta = 1.0 / dotproduct(r, r);
    multiply_sv(p, beta, r);
/* 反復 */
    for (k = 0; k < maxt; k++) {
/* Ap := A * p */
        multiply_mv(Ap, A, p);
/* pAp := (p, Ap) */
        pAp = dotproduct(p, Ap);
/* α = 1/(p,Ap) */
        alpha = 1.0 / pAp;
/* x, r の更新 */
        for (i = 0; i < N; i++) {
/* x = x + α*p */
            x[i] += alpha * p[i];
/* r = r - α*A*p */
            r[i] -= alpha * Ap[i];
        }
/* ||r||<ε||b|| ならば反復を終了 */
        rxr = dotproduct(r,r);
        printf("LOOP : %d\t Error = %g\n", k, rxr);
        if (rxr < eps_b2)
            break;
/* β := 1 / (r, r) */
        beta = 1.0 / rxr;
/* p := r + β*p */
        for (i = 0; i < N; i++){
            p[i] += beta * r[i];
        }
    }
    if (maxiter != NULL)
        *maxiter = k;
}

int main(void)
{
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
//    matrix a = {{5.0/18.0, 0.0}, {0.0, 0.5}};
//    vector b = {1.0/6.0, 0.5};
//    vector x = {3.0,3.0};
    int    i;
    int maxiter;
	long cpu_time_0, cpu_time_1, cg_cpu_time;
	
	cpu_time_0 = clock();
	// printf("initial setting CPU time： %ld\n", cpu_time_0);
    // CG法でAx=bを解く
    cg(TMAX, a, b, x, EPS, &maxiter);
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
