#include <stdio.h>
#include <stdlib.h>
#include <math.h>
//可以对顺序gauss消去法进行改进，但也不保证误差一定很小

int main(void) {

    // 输入矩阵A和b
    double **A = NULL;
    int n;
    printf("Enter the size of the matrix: ");
    scanf("%d", &n);
    
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
    

    
    printf("Enter the number of iterations:\n ");
    double* b = (double*)malloc(n * sizeof(double));
    for (int i = 0; i < n; i++) {
        scanf("%lf", &b[i]);
    }

    /* // 顺序gauss消去法
    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {    
            double factor = A[j][i] / A[i][i];   
            A[j][i] = 0;    
            for (int k = i+1; k < n; k++) {
                A[j][k] -= factor * A[i][k];
            }
            b[j] -= factor * b[i];
        }
    } */

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
            return 1;
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
                

    // 输出变换后的A和b
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
    printf("\n");

    double* x = (double*)malloc(n * sizeof(double));
    // 解线性方程组
    for (int i = n - 1; i >= 0; i--) {
        double ad = b[i];
        for (int j = n-1; j >= i+1; j--) {
            ad -= A[i][j] * x[j];
        }
        x[i] = ad/A[i][i];
    }

    // 输出解
    printf("The solution is:\n");
    for(int i = 0; i < n; i++) {
        printf("%lf\n", x[i]);
    }
    printf("\n");

    for (int i = 0; i < n; i++) {
        free(A[i]);
    }
    free(A);
    free(b);
    return 0;
}