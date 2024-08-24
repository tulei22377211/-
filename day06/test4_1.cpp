//求解超定方程组Ax=b
//其中 计算时用到了列主元高斯消元法

#include <stdio.h>
#include <stdlib.h>
#include <math.h>


void gauss_elimination(int n, double **A_tA, double *A_tb, double *x) {
    for (int i = 0; i < n; i++) {
        int max_index = i;
        for (int j = i + 1; j < n; j++) {//找列最大值
            if (fabs(A_tA[j][i]) > fabs(A_tA[max_index][i])) {
                max_index = j;
            }
        }
        if (A_tA[max_index][i] == 0) {//如果最大值为0，则矩阵为奇异矩阵
            printf("The matrix is singular.\n");
            exit(0);
        }
        if (max_index != i) {//交换两行
            for(int j = i; j < n; j++) {
                double temp = A_tA[i][j];
                A_tA[i][j] = A_tA[max_index][j];
                A_tA[max_index][j] = temp;
            }
            double temp_b = A_tb[i];
            A_tb[i] = A_tb[max_index];
            A_tb[max_index] = temp_b;
        }
        for (int j = i + 1; j < n; j++) {//消去法
            double factor = A_tA[j][i] / A_tA[i][i];
            A_tA[j][i] = 0;
            for (int k = i + 1; k < n; k++) {
                A_tA[j][k] -= factor * A_tA[i][k];
            }
            A_tb[j] -= factor * A_tb[i];
        }
    }
                

    // 解线性方程组
    for (int i = n - 1; i >= 0; i--) {
        double ad = A_tb[i];
        for (int j = n-1; j >= i+1; j--) {
            ad -= A_tA[i][j] * x[j];
        }
        
        x[i] = ad/A_tA[i][i];
    }
}

void result(int n, int num, double **A, double *b, double *x) {
    /* // 输出矩阵A和b
    printf("A = \n");
    for (int i = 0; i < num; i++) {
        for (int j = 0; j < n; j++) {
            printf("%lf ", A[i][j]);
        }
        printf("\n");
    }

    printf("b = \n");
    for (int i = 0; i < num; i++) {
        printf("%lf\n", b[i]);
    }
    printf("\n"); */

    //分配内存
    double **A_tA = (double **)malloc(n * sizeof(double *));
    if (A_tA == NULL) {
        printf("Memory allocation failed.\n");
        exit(0);
    }
    for (int i = 0; i < n; i++) {
        A_tA[i] = (double *)malloc(n * sizeof(double));
        if (A_tA[i] == NULL) {
            printf("Memory allocation failed.\n");
            exit(0);
        }
    }

    double *A_tb = (double *)malloc(n * sizeof(double));
    if (A_tb == NULL) {
        printf("Memory allocation failed.\n");
        exit(0);
    }

    //计算A_tA和A_tb
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            A_tA[i][j] = 0;
        }
    }

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            for (int k = 0; k < num; k++) {
                A_tA[i][j] += A[k][i] * A[k][j];
            }
        }
    }

    for (int i = 0; i < n; i++) {
        A_tb[i] = 0;
    }

    for (int i = 0; i < n; i++) {
        for (int k = 0; k < num; k++) {
            A_tb[i] += A[k][i] * b[k];
        }
    }

    /* //输出A_tA和A_tb
    printf("A_tA = \n");
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            printf("%lf ", A_tA[i][j]);
        }
        printf("\n");
    }
    printf("\n");
    
    printf("A_tb = \n");
    for (int i = 0; i < n; i++) {
        printf("%lf\n", A_tb[i]);
    }
    printf("\n"); */


    //调用高斯消元法求解线性方程组Ax=b
    gauss_elimination(n, A_tA, A_tb, x);
    

    //释放内存
    for (int i = 0; i < n; i++) {
        free(A_tA[i]);
    }
    free(A_tA);
    free(A_tb);
}

double S2(double x,double y,double x1,double y1) {
    return sqrt((x-x1)*(x-x1) + (y-y1)*(y-y1));
}


int main(void) {
    
    int n = 3;//列数//维数+1
    int num = 4; //行数/点数

    double eps = 1e-5; //误差

    double k = 1;//系数

    int max_iter = 100; //最大迭代次数
    int iter = 0; //迭代次数

    //分配内存
    double ** A = (double **)malloc(num * sizeof(double *));
    if (A == NULL) {
        printf("Memory allocation failed.\n");
        exit(0);
    }
    for (int i = 0; i < num; i++) {
        A[i] = (double *)malloc(n * sizeof(double));
        if (A[i] == NULL) {
            printf("Memory allocation failed.\n");
            exit(0);
        }
    }

    double * b = (double *)malloc(num * sizeof(double));
    if (b == NULL) {
        printf("Memory allocation failed.\n");
        exit(0);
    }

    double* x = (double*)malloc(num * sizeof(double));
    if (x == NULL) {
        printf("Memory allocation failed.\n");
        exit(0);
    }

    double* y = (double*)malloc(num * sizeof(double));
    if (y == NULL) {
        printf("Memory allocation failed.\n");
        exit(0);
    }

    double* r = (double*)malloc(num * sizeof(double));
    if (r == NULL) {
        printf("Memory allocation failed.\n");
        exit(0);
    }

    double* x_result = (double*)malloc(n * sizeof(double));
    if (x_result == NULL) {
        printf("Memory allocation failed.\n");
        exit(0);
    }

    x[0] = -1; x[1] = 1; x[2] = 1; x[3] = 0;
    y[0] = 0; y[1] = 0.5; y[2] = -0.5; y[3] = 1;
    r[0] = 1; r[1] = 0.5; r[2] = 0.5; r[3] = 0.5;


    // 计算A和b
    double* s = (double*)malloc(num * sizeof(double));
    if (s == NULL) {
        printf("Memory allocation failed.\n");
        exit(0);
    }

    x_result[0] = 0; x_result[1] = 0; x_result[2] = 0;

    
    double result_1 = sqrt((x_result[0]-x[0])*(x_result[0]-x[0]) + (x_result[1]-y[0])*(x_result[1]-y[0])) + sqrt((x_result[0]-x[1])*(x_result[0]-x[1]) + (x_result[1]-y[1])*(x_result[1]-y[1])) + sqrt((x_result[0]-x[2])*(x_result[0]-x[2]) + (x_result[1]-y[2])*(x_result[1]-y[2])) + sqrt((x_result[0]-x[3])*(x_result[0]-x[3]) + (x_result[1]-y[3])*(x_result[1]-y[3]))-r[0]-r[1]-r[2]-r[3]-4*x_result[2];

    printf("\n");
    
    printf("----------------------------------------------------\n");
    printf("iter\tx\t\ty\t\tK\t\terror\n");
    printf("----------------------------------------------------\n");
    printf("%d\t%lf\t%lf\t%lf\t%lf\n", iter, x_result[0], x_result[1], x_result[2], result_1);
    iter++;

    
    while (iter < max_iter && fabs(result_1) > eps) {    
        for (int i = 0; i < num; i++) {
            //printf("xi = %lf, yi = %lf, ri = %lf\n", x[i], y[i], r[i]);
            s[i] = S2(x_result[0],x_result[1],x[i],y[i]);
        }
        //printf("s1 = %lf, s2 = %lf, s3 = %lf, s4 = %lf \n", s[0], s[1], s[2], s[3]);

        for (int i = 0; i < num; i++) {
            for (int j = 0; j < n; j++) {
                if (j == 0) {
                    A[i][j] = (x_result[j] - x[i])/s[i];
                } else if (j == 1) {
                    A[i][j] = (x_result[j] - y[i])/s[i];
                } else {
                    A[i][j] = -1*k;
                }
            }
        }
        
        for (int i = 0; i < num; i++) {
            b[i] = -1*(s[i] - (r[i] + x_result[2]));
        }

        double* v = (double*)malloc(n * sizeof(double));
        if (v == NULL) {
            printf("Memory allocation failed.\n");
            exit(0);
        }

        result(n, num, A, b, v);
        
        for (int i = 0; i < n; i++) {
            x_result[i] += v[i];
        }

        result_1 = sqrt((x_result[0]-x[0])*(x_result[0]-x[0]) + (x_result[1]-y[0])*(x_result[1]-y[0])) + sqrt((x_result[0]-x[1])*(x_result[0]-x[1]) + (x_result[1]-y[1])*(x_result[1]-y[1])) + sqrt((x_result[0]-x[2])*(x_result[0]-x[2]) + (x_result[1]-y[2])*(x_result[1]-y[2])) + sqrt((x_result[0]-x[3])*(x_result[0]-x[3]) + (x_result[1]-y[3])*(x_result[1]-y[3]))-r[0]-r[1]-r[2]-r[3]-4*x_result[2];

        printf("%d\t%lf\t%lf\t%lf\t%lf\n", iter, x_result[0], x_result[1], x_result[2], result_1);
        iter++;
    }

    // 输出解
    printf("\n");
    printf("The solution is:\n");
    for(int i = 0; i < n-1; i++) {
        printf("x%d = %lf\n",i+1, x_result[i]);
    }
    printf("K = %lf\n", x_result[n-1]);
    printf("\n");

    // 释放内存
    for (int i = 0; i < num; i++) {
        free(A[i]);
    }
    free(A);
    free(b);
    free(x);
    free(y);
    free(r);
    free(x_result);
    free(s);

    return 0;
}
