#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main(void) {
    float **A = NULL;
    int n = 2;

    A = (float **)malloc(n * sizeof(float *));
    if(A == NULL) {
        printf("Memory allocation failed.\n");
        return 1;
    }

    for (int i = 0; i < n; i++) {
        A[i] = (float *)malloc(n * sizeof(float));
    }

    A[0][0] = 0.000000010;A[0][1] = 1.00;
    A[1][0] = 1.00;A[1][1] = 1.00;

    float* b = (float*)malloc(n * sizeof(float));
    b[0] = 1.00;
    b[1] = 2.00;

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

    /* // 顺序gauss消去法
    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {    
            float factor = A[j][i] / A[i][i];   
            A[j][i] = 0;    
            for (int k = i+1; k < n; k++) {
                A[j][k] -= factor * A[i][k];
            }
            b[j] -= factor * b[i];
        }
    } */

    //列主元高斯消去法
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
                float temp = A[i][j];
                A[i][j] = A[max_index][j];
                A[max_index][j] = temp;
            }
            float temp_b = b[i];
            b[i] = b[max_index];
            b[max_index] = temp_b;
        }
        for (int j = i + 1; j < n; j++) {//消去法
            float factor = A[j][i] / A[i][i];
            A[j][i] = 0;
            for (int k = i + 1; k < n; k++) {
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