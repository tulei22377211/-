#include <stdio.h>
#include <stdlib.h>

int main(void) {
    float **A = NULL;
    int n = 3;

    A = (float **)malloc(n * sizeof(float *));
    if(A == NULL) {
        printf("Memory allocation failed.\n");
        return 1;
    }

    for (int i = 0; i < n; i++) {
        A[i] = (float *)malloc(n * sizeof(float));
    }

    A[0][0] = 0.012;A[0][1] = 0.01;A[0][2] = 0.167;
    A[1][0] = 1;A[1][1] = 0.8334;A[1][2] = 5.91;
    A[2][0] = 3200;A[2][1] = 1200;A[2][2] = 4.2;


    float* b = (float*)malloc(n * sizeof(float));
    b[0] = 0.6781;
    b[1] = 12.1;
    b[2] = 981;

    /* printf("Enter the size of the matrix: ");
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
    } */

    // 顺序gauss消去法
    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {    
            float factor = A[j][i] / A[i][i];   
            A[j][i] = 0;    
            for (int k = i+1; k < n; k++) {
                A[j][k] -= factor * A[i][k];
            }
            b[j] -= factor * b[i];
        }
    }

    // 输出变换后的A和b
    printf("The transformed matrix is:\n");
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            printf("%f ", A[i][j]);
        }
        printf("\n");
    }
    printf("The transformed vector is:\n");
    printf("\n");
    for(int i = 0; i < n; i++) {
        printf("%f ", b[i]);
    }
    printf("\n");

    float* x = (float*)malloc(n * sizeof(float));
    // 解线性方程组
    for (int i = n - 1; i >= 0; i--) {
        float ad = b[i];
        for (int j = n-1; j >= i+1; j--) {
            ad -= A[i][j] * x[j];
        }
        x[i] = ad/A[i][i];
    }

    // 输出解
    printf("The solution is:\n");
    for(int i = 0; i < n; i++) {
        printf("%f\n", x[i]);
    }
    printf("\n");

    for (int i = 0; i < n; i++) {
        free(A[i]);
    }
    free(A);
    free(b);
    return 0;
}