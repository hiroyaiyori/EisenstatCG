// CG法テストプログラム
/*==================================================*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define N    10
#define TMAX 100
#define EPS  (1.0e-6)

typedef double Vector[N];       // ベクトル
typedef double Matrix[N][N];    // 行列

// ベクトル初期化
void init_vector(Vector x)
{
  int i;
  for(i = 0; i < N; i++){
    x[i] = 0;
  }
}

// ベクトルに行列を作用 y = Ax
void vector_x_matrix(Vector y, Matrix a, Vector x)
{
  int    i, j;
  double vxm;
  for(i = 0; i < N; i++){
    vxm = 0;
    for(j = 0; j < N; j++){
      vxm += a[i][j]*x[j];
    }
    y[i] = vxm;
  }
}

// 内積を計算
double dot_product(Vector x, Vector y)
{
  int    i;
  double dot_p = 0;
  for(i = 0; i < N; i++){
    dot_p += x[i]*y[i];
  }
  return dot_p;
}

// ベクトルノルムを計算
// ベクトルノルム := sgm(0〜N-1)|x[i]|
double vector_norm(Vector x)
{
  int    i;
  double norm = 0;
  for(i = 0; i < N; i++){
    norm += fabs(x[i]);
  }
  return norm;
}

// CG法
void cg_method(Matrix a, Vector x, Vector b)
{
  static Vector p;   // 探索方向ベクトル
  static Vector r;   // 残差ベクトル
  static Vector ax;
  static Vector ap;
  static double pre_rxr;
  int           i, iter;

  // Axを計算
  vector_x_matrix(ax, a, x);

  // pとrを計算 p = r := b - Ax
  for(i = 0; i < N; i++){
    p[i] = b[i] - ax[i];
    r[i] = p[i];
  }

  // 反復計算
  for(iter = 1; iter < TMAX; iter++){
    double alpha, beta, err = 0;

    // alphaを計算
    vector_x_matrix(ap, a, p);
	pre_rxr = dot_product(r, r);
    alpha = pre_rxr/dot_product(p, ap);
    
    for(i = 0; i < N; i++){
      x[i] += +alpha*p[i];
      r[i] += -alpha*ap[i];
    }
    
    err = vector_norm(r);   // 誤差を計算
    printf("LOOP : %d\t Error = %g\t beta = %g\n", iter, err, beta);
    if(EPS > err) break;

    // EPS < err ならbetaとpを計算してループを継続
    beta = dot_product(r, r)/pre_rxr;
    for(i = 0; i < N; i++){
      p[i] = r[i] + beta*p[i];
    } 
  }
}

int main(void)
{
  // 連立方程式 Ax = b
  // 行列Aは正定値対象行列
    Matrix a = { {5.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
                 {2.0, 5.0, 2.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0},
                 {0.0, 2.0, 5.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
                 {0.0, 0.0, 2.0, 5.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0},
                 {0.0, 0.0, 0.0, 2.0, 5.0, 2.0, 0.0, 0.0, 0.0, 0.0},
                 {0.0, 0.0, 0.0, 0.0, 2.0, 5.0, 2.0, 0.0, 0.0, 0.0},
                 {0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 5.0, -2.0, 0.0, 0.0},
                 {0.0, 2.0, 0.0, 0.0, 0.0, 0.0, -2.0, 5.0, 2.0, 0.0},
                 {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 5.0, 2.0},
                 {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 5.0} };
    Vector b = {3.0, 1.0, 4.0, 0.0, 5.0, -1.0, 6.0, -2.0, 7.0, -15.0};
    // 初期値を設定
    Vector x = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
    int    i;

  // CG法でAx=bを解く
  cg_method(a, x, b);
  Vector ans_x = {0.637207, -0.093018, 1.42197, -1.4619, 2.23279, -1.62006, 1.31737, -1.32663, 3.72697, -4.49079};
  
  printf("--------------- calc done -------------\n");
  double err_sum = 0.0;
  for(i = 0; i < N; i++){
    err_sum += (x[i] - ans_x[i]) * (x[i] - ans_x[i]);
    printf("x[%d] = %2g\n", i, x[i]);
  }
  printf("error %2g\n", err_sum);

    return 0;
}
