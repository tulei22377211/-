// 线性拟合五个点 输出y = ax + b中的a,b

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void result(double *x_i, double *y_i, int num, int n, double *a, double *b) {

    double ** A = (double **)malloc(num * sizeof(double *));
    //分配内存
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

    double * B = (double *)malloc(num * sizeof(double));
    if (B == NULL) {
        printf("Memory allocation failed.\n");
        exit(0);
    }

    // 构造矩阵A和b
    for (int i = 0; i < num; i++) {
        for (int j = 0; j < n; j++) {
            if (j == 0) {
                A[i][j] = 1;
            } else if (j == 1) {
                A[i][j] = x_i[i];
            } else {
                A[i][j] = 0;
            }
        }
        B[i] = y_i[i];
    }

    printf("A = \n");
    for (int i = 0; i < num; i++) {
        for (int j = 0; j < n; j++) {
            printf("%lf ", A[i][j]);
        }
        printf("\n");
    }

    printf("b = \n");
    for (int i = 0; i < num; i++) {
        printf("%lf\n", B[i]);
    }
    printf("\n");

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

    printf("A_tA = \n");
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            printf("%lf ", A_tA[i][j]);
        }
        printf("\n");
    }

    double *A_tb = (double *)malloc(n * sizeof(double));
    if (A_tb == NULL) {
        printf("Memory allocation failed.\n");
        exit(0);
    }

    for (int i = 0; i < n; i++) {
        A_tb[i] = 0;
        
    }


    for (int i = 0; i < n; i++) {
        
        for (int k = 0; k < num; k++) {
            A_tb[i] += A[k][i] * B[k];
        }
    }
    
    printf("A_tb = \n");
    for (int i = 0; i < n; i++) {
        printf("%lf\n", A_tb[i]);
    }
    printf("\n");

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
                

    

    double* x = (double*)malloc(n * sizeof(double));
    // 解线性方程组
    for (int i = n - 1; i >= 0; i--) {
        double ad = A_tb[i];
        for (int j = n-1; j >= i+1; j--) {
            ad -= A_tA[i][j] * x[j];
        }
        x[i] = ad/A_tA[i][i];
    }

    // 输出解
    printf("The solution is:\n");
    for(int i = 0; i < n; i++) {
        printf("%lf\n", x[i]);
    }
    printf("\n");
    
    *a = x[1];
    *b = x[0];
    
    // free memory
    for (int i = 0; i < num; i++) {
        free(A[i]);
    }
    free(b);
    for (int i = 0; i < n; i++) {
        free(A_tA[i]);
    }
    free(A_tA);
    free(A_tb);
    free(x);

}


// 线性拟合五个点 输出y = ax + b中的a,b
int main(void) {
    //列数
    int n = 2;
    //点数
    int num = 5;

    double * x_i = (double *)malloc(num * sizeof(double));
    if (x_i == NULL) {
        printf("Memory allocation failed.\n");
        return -1;
    }
    
    double * y_i = (double *)malloc(num * sizeof(double));
    if (y_i == NULL) {
        printf("Memory allocation failed.\n");
        return -1;
    }

    x_i[0] = 25; x_i[1] = 27; x_i[2] = 31; x_i[3] = 33; x_i[4] = 35;
    y_i[0] = 110; y_i[1] = 115; y_i[2] = 155; y_i[3] = 160; y_i[4] = 180;

    double a = 0, b = 0;

    result(x_i, y_i, num, n, &a, &b);

    printf("a = %lf\n", a);
    printf("b = %lf\n", b);

    // free memory
    free(x_i);
    free(y_i);

    return 0;
}
