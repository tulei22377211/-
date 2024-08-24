#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

void jisuan_zhijie(int n, double * er_1, double * er_2, double * time_used) {
    clock_t start, end;
    
    start = clock();
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
    

    //列主元Gauss消去法
    for (int i = 0; i < n; i++) {
        int max_index = i;
        for (int j = i + 1; j < n; j++) {//找列最大值
            if (fabs(A[j][i]) > fabs(A[max_index][i])) {
                max_index = j;
            }
        }
        if (A[max_index][i] == 0) {//如果最大值为0，则矩阵为奇异矩阵
            printf("The matrix is singular.\n");
           
        }
        if (max_index != i) {//交换两行
            for(int j = i; j < n; j++) {
                double temp = A[i][j];
                A[i][j] = A[max_index][j];
                A[max_index][j] = temp;
            }
            double temp_b = b[i];
            b[i] = b[max_index];
            b[max_index] = temp_b;
        }
        for (int j = i + 1; j < n; j++) {//消去法
            double factor = A[j][i] / A[i][i];
            A[j][i] = 0;
            for (int k = i + 1; k < n; k++) {
                A[j][k] -= factor * A[i][k];
            }
            b[j] -= factor * b[i];
        }
    }
                

    /* / 输出变换后的A和b
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            printf("%lf ", A[i][j]);
        }
        printf("\n");
    }
    printf("\n");
    for(int i = 0; i < n; i++) {
        printf("%lf ", b[i]);
    }
    printf("\n"); */

    double* x = (double*)malloc(n * sizeof(double));
    if(x == NULL) {
        printf("Memory allocation failed.\n");
        exit(1);
    }

    // 解线性方程组
    for (int i = n - 1; i >= 0; i--) {
        double ad = b[i];
        for (int j = n-1; j >= i+1; j--) {
            ad -= A[i][j] * x[j];
        }
        x[i] = ad/A[i][i];
    }

    /* // 输出解
    printf("The solution is:\n");
    for(int i = 0; i < n; i++) {
        printf("%lf\n", x[i]);
    }
    printf("\n"); */
    
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

    double err_1 = fabs(x[1] - u[1]);
    for (int i = 2; i < n; i++) {
        err_1 = err_1 > fabs(x[i] - u[i]) ? err_1 : fabs(x[i] - u[i]);
    }

    double err_2 = 0;
    for (int i = 1; i < n; i++) {
        err_2 += h*(x[i] - u[i])*(x[i] - u[i]);
    }
    err_2 = sqrt(err_2);

    *er_2 = err_2;
    *er_1 = err_1;

    // 计算运行时间
    end = clock();
    *time_used = (double)(end - start) / CLOCKS_PER_SEC;

    // 释放内存
    for (int i = 0; i < n; i++) {
        free(A[i]);
    }
    free(A);
    free(b);
    free(x);
    free(u);
}

void jisuan_diedai(int n, double * er_1, double * er_2, double * time_used) {
    clock_t start, end;
    start = clock();

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

    // 迭代法求解


    //设置数列来存储n维向量的n个分量迭代结果
    
    // 定义最大迭代次数
    int max_iter = 999999;

    double** x = (double**)malloc((max_iter+1) * sizeof(double*));
    for (int i = 0; i < (max_iter+1); i++) {
        x[i] = (double*)malloc(n * sizeof(double));
    }
    if(x == NULL) {
        printf("Memory allocation failed.\n");
        exit(1);
    }

    // 初始化n维向量的n个分量为0'
    for (int i = 0; i < n; i++) {
        x[0][i] = 0.0;
    }

    //存储L2范数的数组
    double* diff = (double*)malloc((max_iter+1) * sizeof(double));
    if(diff == NULL) {
        printf("Memory allocation failed.\n");
        exit(1);
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
            x[iter][i] = (((b[i]-sum) / A[i][i] ) - x[iter-1][i]) * omiga + x[iter-1][i];
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

    *er_2 = err_2;
    *er_1 = err_1;

    // 计算运行时间
    end = clock();
    *time_used = (double)(end - start) / CLOCKS_PER_SEC;

    // 释放内存
    for (int i = 0; i < (max_iter+1); i++) {
        free(x[i]);
    }
    free(x);
    free(diff);
    for (int i = 0; i < n; i++) {
        free(A[i]);
    }
    free(A);
    free(b);
    free(u);
}



int main(void) {
    int n1 = 10;
    int n2 = 20;
    int n3 = 40;
    int n4 = 80;
    int n5 = 160;

    double err_1_1 = 0;
    double err_2_1 = 0;
    double err_1_2 = 0;
    double err_2_2 = 0;
    double err_1_3 = 0;
    double err_2_3 = 0;
    double err_1_4 = 0;
    double err_2_4 = 0;
    double err_1_5 = 0;
    double err_2_5 = 0;

    double time_used_1 = 0;
    double time_used_2 = 0;
    double time_used_3 = 0;
    double time_used_4 = 0;
    double time_used_5 = 0;
    
    jisuan_zhijie(n1, &err_1_1, &err_2_1, &time_used_1);
    jisuan_zhijie(n2, &err_1_2, &err_2_2, &time_used_2);
    jisuan_zhijie(n3, &err_1_3, &err_2_3, &time_used_3);
    jisuan_zhijie(n4, &err_1_4, &err_2_4, &time_used_4);
    jisuan_zhijie(n5, &err_1_5, &err_2_5, &time_used_5);

    printf("Direct solution method\n");
    printf("----------\n");
    printf("index\tn\terror_1\t\terror_2\t\ttime_used\n");
    printf("--------------------------\n");
   
    printf("n1\t10\t%lf\t%lf\t%lf\n", err_1_1, err_2_1, time_used_1);
    printf("n2\t20\t%lf\t%lf\t%lf\n", err_1_2, err_2_2, time_used_2);
    printf("n3\t40\t%lf\t%lf\t%lf\n", err_1_3, err_2_3, time_used_3);
    printf("n4\t80\t%lf\t%lf\t%lf\n", err_1_4, err_2_4, time_used_4);
    printf("n5\t160\t%lf\t%lf\t%lf\n", err_1_5, err_2_5, time_used_5);
    printf("--------------------------\n");

    printf("\n");
    
    jisuan_diedai(n1, &err_1_1, &err_2_1, &time_used_1);
    jisuan_diedai(n2, &err_1_2, &err_2_2, &time_used_2);
    jisuan_diedai(n3, &err_1_3, &err_2_3, &time_used_3);  
    jisuan_diedai(n4, &err_1_4, &err_2_4, &time_used_4);
    jisuan_diedai(n5, &err_1_5, &err_2_5, &time_used_5);

    printf("Iterative solution method\n");
    printf("----------\n");
    printf("index\tn\terror_1\t\terror_2\t\ttime_used\n");
    printf("--------------------------\n");
    
    printf("n1\t10\t%lf\t%lf\t%lf\n", err_1_1, err_2_1, time_used_1);
    printf("n2\t20\t%lf\t%lf\t%lf\n", err_1_2, err_2_2, time_used_2);
    printf("n3\t40\t%lf\t%lf\t%lf\n", err_1_3, err_2_3, time_used_3);
    printf("n4\t80\t%lf\t%lf\t%lf\n", err_1_4, err_2_4, time_used_4);
    printf("n5\t160\t%lf\t%lf\t%lf\n", err_1_5, err_2_5, time_used_5);
    printf("----------\n");
    return 0;
}