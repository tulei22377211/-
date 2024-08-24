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
    // 输出矩阵A和b
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
    printf("\n");

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

    //输出A_tA和A_tb
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
    printf("\n");


    //调用高斯消元法求解线性方程组Ax=b
    gauss_elimination(n, A_tA, A_tb, x);

    //释放内存
    for (int i = 0; i < n; i++) {
        free(A_tA[i]);
    }
    free(A_tA);
    free(A_tb);
}

int main(void) {
    
    int n = 3;//列数
    int num = 4; //行数

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

    A[0][0] = 1; A[0][1] = 2; A[0][2] = 4;
    A[1][0] = 2; A[1][1] = 1; A[1][2] = 1;
    A[2][0] = 1; A[2][1] = 1; A[2][2] = 2;
    A[3][0] = 1; A[3][1] = -1; A[3][2] = -2;
    
    b[0] = -1; 
    b[1] = 4; 
    b[2] = 2; 
    b[3] = 1;

    double* x = (double*)malloc(n * sizeof(double));
    if (x == NULL) {
        printf("Memory allocation failed.\n");
        exit(0);
    }

    result(n, num, A, b, x);

    // 输出解
    printf("The solution is:\n");
    for(int i = 0; i < n; i++) {
        printf("%lf\n", x[i]);
    }
    printf("\n");

    // 释放内存
    for (int i = 0; i < num; i++) {
        free(A[i]);
    }
    free(A);
    free(b);
    free(x);    

    return 0;
}
