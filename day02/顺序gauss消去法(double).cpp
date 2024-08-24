#include <stdio.h>
#include <stdlib.h>
//适用于对角线严格占优的矩阵


int main(void) {
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
    

    
    printf("Enter the number of iterations:\n");
    double* b = (double*)malloc(n * sizeof(double));
    for (int i = 0; i < n; i++) {
        scanf("%lf", &b[i]);
    }

    // 顺序gauss消去法
    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {    
            double factor = A[j][i] / A[i][i];   
            A[j][i] = 0;    
            for (int k = i+1; k < n; k++) {
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