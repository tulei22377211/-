#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "error.hpp"

double max(double a, double b, double c) {
    return a > b? (a > c? a : c) : (b > c? b : c);
}
    
    
int main(void) {

    // 输入矩阵大小n
    int n;
    printf("Enter the size of the matrix: ");
    scanf("%d", &n);

    /* //输入最大迭代次数
    int max_iter;
    printf("Enter the maximum number of iterations: ");
    scanf("%d", &max_iter); */
    
    // 输入A
    double **A = NULL;
    A = (double **)malloc(n * sizeof(double *));
    for (int i = 0; i < n; i++) {
        A[i] = (double *)malloc(n * sizeof(double));
    }
    if(A == NULL) {
        printf("Memory allocation failed.\n");
        return 1;
    }
    printf("Enter the elements of the matrix:\n");
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            scanf("%lf", &A[i][j]);
        }
        printf("\n");
    }

    // 输入b
    printf("Enter the number of iterations:\n");
    double* b = (double*)malloc(n * sizeof(double));
    if(b == NULL) {
        printf("Memory allocation failed.\n");
        return 1;
    }
    for (int i = 0; i < n; i++) {
        scanf("%lf", &b[i]);
    }

    // 进行Jacobi迭代

    //计算D矩阵
    double** D = (double**)malloc(n * sizeof(double*));
    for (int i = 0; i < n; i++) {
        D[i] = (double*)malloc(n * sizeof(double));
    }
    if(D == NULL) {
        printf("Memory allocation failed.\n");
        return 1;
    }

    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++) {
            if (i == j) {
                D[i][j] = A[i][j];
            } else {
                D[i][j] = 0.0;
            }
        }
    }

    /* //计算L矩阵
    double** L = (double**)malloc(n * sizeof(double*));
    for (int i = 0; i < n; i++) {
        L[i] = (double*)malloc(n * sizeof(double));
    }
    if(L == NULL) {
        printf("Memory allocation failed.\n");
        return 1;
    }

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (i == j) {
                L[i][j] = 0.0;
            } else if (j > i) {
                L[i][j] = -A[i][j];
            } else {
                L[i][j] = 0.0;
            }
        }
    }
    
    //计算U矩阵
    double** U = (double**)malloc(n * sizeof(double*));
    for (int i = 0; i < n; i++) {
        U[i] = (double*)malloc(n * sizeof(double));
    }
    if(U == NULL) {
        printf("Memory allocation failed.\n");
        return 1;
    }

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (i == j) {
                U[i][j] = 0.0;
            } else if (j < i) {
                U[i][j] = -A[i][j];
            } else {
                U[i][j] = 0.0;
            }
        }
    } */

    /* // 输出矩阵A
    printf("The matrix A is:\n");
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            printf("%lf\t", A[i][j]);
        }
        printf("\n");
    }
    
    // 输出矩阵D
    printf("The matrix D is:\n");
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            printf("%lf\t", D[i][j]);
        }
        printf("\n");
    }

    // 输出矩阵b
    printf("The vector b is:\n");
    for (int i = 0; i < n; i++) {
        printf("%lf\t", b[i]);
    }
    printf("\n");
 */
    //设置数列来存储n维向量的n个分量迭代结果
    
    // 定义最大迭代次数
    int max_iter = 50;

    double** x = (double**)malloc((max_iter+1) * sizeof(double*));
    for (int i = 0; i < (max_iter+1); i++) {
        x[i] = (double*)malloc(n * sizeof(double));
    }
    if(x == NULL) {
        printf("Memory allocation failed.\n");
        return 1;
    }

    // 初始化n维向量的n个分量为0'
    for (int i = 0; i < n; i++) {
        x[0][i] = 0.0;
    }

    //存储L2范数的数组
    double* diff = (double*)malloc((max_iter+1) * sizeof(double));
    if(diff == NULL) {
        printf("Memory allocation failed.\n");
        return 1;
    }
    // 定义初始的L2范数
    diff[0] = 9999;

    // 开始迭代
    int iter = 1;
    double eps = 1e-5;
    
    // 开始迭代
    // 定义结果判断是否收敛
    double omiga = 1.46;
    int result = 0;
    while (iter <= max_iter) {
        for (int i = 0; i < n; i++) {
            double sum = 0.0;
            for (int j = 0; j < n; j++) {
                if (j != i) {
                    if (j < i) {
                        sum += A[i][j] * x[iter][j];
                    } else {
                        sum += A[i][j] * x[iter-1][j];
                    }
                } 
            }
            x[iter][i] = (((b[i]-sum) / D[i][i] ) - x[iter-1][i]) * omiga + x[iter-1][i];
        }
        
        // L2范数计算
        diff[iter] = 0.0;
        for (int i = 0; i < n; i++) {
            diff[iter] += (x[iter][i] - x[iter-1][i])*(x[iter][i] - x[iter-1][i]);
        }
        diff[iter] = sqrt(diff[iter]);
        
        
        // 输出当前迭代结果
        printf("Iteration %d:\n", iter);
        printf("x = [");
        for (int i = 0; i < n; i++) {
            printf("%lf ", x[iter][i]);
        }
        printf("]\n");
        printf("L2_error = %lf\n", diff[iter]);
        // 判断是否收敛
        if (diff[iter] < eps) {
            result = 1;
            break;
        } 
        iter++;
    }
    if(iter > max_iter) {
        iter--;
    }
    printf("----------\n");
    
    double* x_real = (double*)malloc(n * sizeof(double));
    if(x_real == NULL) {
        printf("Memory allocation failed.\n");
        return 1;
    }
    x_real[0] = 2;
    x_real[1] = 1;
    x_real[2] = -1;
    
    
    // 输出迭代过程
    printf("--------------------------------------------\n");
    printf("inter\tx1\t\tx2\t\tx3\t\tL_infinity_error\t\tL_2_error\n");
    printf("--------------------------------------------\n");
    printf("%d\t%lf\t%lf\t%lf\t%lf\t%s\n", 0, x[0][0], x[0][1], x[0][2],max(fabs(x[0][0]-x_real[0]),fabs(x[0][1]-x_real[1]),fabs(x[0][2])-x_real[2]), "initial error");
    for (int i = 1; i <= iter; i++){
        printf("%d\t%lf\t%lf\t%lf\t%lf\t%lf\n", i, x[i][0], x[i][1], x[i][2],max(fabs(x[i][0]-1),fabs(x[i][1]-1),fabs(x[i][2])-1), diff[i]);
    }
    printf("--------------------------------------------\n");
    
// 输出结果
    if (result == 1) {
        printf("The result was found within %d iterations.\n", iter);
    } else {
        //告诉用户结果未收敛
        printf("The result was not found within %d iterations.\n", max_iter);
        
        // 输出绝对误差(需要知道真正结果)
        printf("The absolute error is (%lf,%lf,%lf)'.\n",absolute_error(x[iter][0],1),absolute_error(x[iter][1],1),absolute_error(x[iter][2],1));

        // 输出x带入方程计算的L2范数误差（不需要知道真正结果）
        double *x_dairu = (double*)malloc(n * sizeof(double));
        if(x_dairu == NULL) {
            printf("Memory allocation failed.\n");
            return 1;
        }
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                x_dairu[i] += A[i][j] * x[max_iter][j];
            }
            x_dairu[i] -= b[i];
        }
        double error_norm = 0.0;
        for (int i = 0; i < n; i++) {
            error_norm += (x_dairu[i] * x_dairu[i]);
        }
        error_norm = sqrt(error_norm);
        printf("The L2 norm error when calculate result brought into the equation is %lf.\n", error_norm);
    }
       

    // 释放内存
    for (int i = 0; i < n; i++) {
        free(A[i]);
        free(D[i]);
        free(x[i]);
    }
    free(A);
    free(D);
    free(x);
    free(b);
    free(diff);

    return 0;
}