#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "error.hpp"
#include <time.h>

double max(double a, double b, double c) {
    return a > b? (a > c? a : c) : (b > c? b : c);
}
    
    
int main(void) {
    clock_t start, end;
    start = clock();
    double time_used = 0;



    // 输入矩阵大小n
    int n = 160;

    // 输入矩阵A和b
    double pi = 3.14159265358979323846;

    /* int n;
    printf("Enter the size of the matrix: ");
    scanf("%d", &n); */
    
    double h = (double)1.0 / n;
    n=n+1;
    // 输入矩阵A和b
    double **A = NULL;
    
    
    A = (double **)malloc(n * sizeof(double *));
    for (int i = 0; i < n; i++) {
        A[i] = (double *)malloc(n * sizeof(double));
    }
    if(A == NULL) {
        printf("Memory allocation failed.\n");
        exit(1);
    }
    
    
    
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            A[i][j] = 0;
        }
    }
    A[0][0]=1;
    A[n-1][n-1] = 1;
    for (int i = 1; i < n-1; i++) {
        A[i][i-1] = (double)(2.0) + h*(1.0+i*h)*(1.0+i*h);
        A[i][i] =  -1*(double)4.0 - 2*h*h * exp(-1*i*h);
        A[i][i+1] = 2.0 - h*(1+i*h)*(1+i*h);
    }

    double* b = (double*)malloc(n * sizeof(double));
    if(b == NULL) {
        printf("Memory allocation failed.\n");
        exit(1);
    }

    for (int i = 0; i < n; i++) {
        b[i] = 2 *h*h*(  (1-(1+i*h)*(1+i*h))*exp(i*h) - pi*pi*cos(pi*i*h) + pi*(1+i*h)*(1+i*h)*sin(pi*i*h)  - 1 -exp(-1*i*h)*cos(pi*i*h)   );
    }
    b[0] = 2;
    b[n-1] = exp(1)-1;

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


    //设置数列来存储n维向量的n个分量迭代结果
    
    // 定义最大迭代次数
    int max_iter = 999999;

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
    double eps = 1e-7;
    
    // 开始迭代
    // 定义结果判断是否收敛
    double omiga = 1.959;
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


    n=n-1;
    // 计算误差
    double * u = (double*)malloc((n+1) * sizeof(double));
    if(u == NULL) {
        printf("Memory allocation failed.\n");
        exit(1);
    }

    for (int i = 0; i < n+1; i++) {
        u[i] = exp(i*h) + cos(pi*i*h);
    }
    
    double err_1 = fabs(x[iter][1] - u[1]);
    for (int i = 2; i < n; i++) {
        err_1 = err_1 > fabs(x[iter][i] - u[i]) ? err_1 : fabs(x[iter][i] - u[i]);
    }

    double err_2 = 0;
    for (int i = 1; i < n; i++) {
        err_2 += h*(x[iter][i] - u[i])*(x[iter][i] - u[i]);
    }
    err_2 = sqrt(err_2);
    
    
    printf("The result is:\n");
    printf("error_1 = %lf\n", err_1);
    printf("error_2 = %lf\n", err_2);
    // 输出运行时间
    end = clock();
    time_used = (double)(end - start) / CLOCKS_PER_SEC;
    printf("The time used is %lf seconds.\n", time_used);
    
    
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