#include <stdio.h>
#include <stdlib.h>

double newton_1(double x, double * x_now, double * y_now, int n) {
    //假设x是递增的,并且x不与已知的x_now重复,且x在x_now的范围内
    
    double ** A = (double**)malloc((n+1) * sizeof(double*));
    if (A == NULL) {
        printf("Memory allocation failed.\n");
        exit(1);
    }
    for(int i=0;i<=n;i++) {
        A[i] = (double*)malloc((n+2) * sizeof(double));
        if (A[i] == NULL) {
            printf("Memory allocation failed.\n");
            exit(1);
        }
    }


    for(int i=0;i<=n;i++) {
        for(int j=0;j<=n+1;j++) {
            A[i][j] = 0;
        }
    }

    for (int i = 0; i <= n; i++) {
        A[i][0] = x_now[i];
        A[i][1] = y_now[i];
    }

    for (int j = 2; j <= n+1; j++) {
        for (int i = j-1; i <= n; i++) {
            A[i][j] = (A[i][j-1] - A[i-1][j-1]) / (A[i][0] - A[i-(j-1)][0]);
        }
    }

    double result = A[0][1] + A[1][2] *(x - x_now[0]);   
    //释放内存
    for(int i=0;i<=n;i++) {
        free(A[i]);
    }
    free(A);


    return result;
}

double newton_2(double x, double * x_now, double * y_now, int n) {
    //假设x是递增的,并且x不与已知的x_now重复,且x在x_now的范围内

    double ** A = (double**)malloc((n+1) * sizeof(double*));
    if (A == NULL) {
        printf("Memory allocation failed.\n");
        exit(1);
    }
    for(int i=0;i<=n;i++) {
        A[i] = (double*)malloc((n+2) * sizeof(double));
        if (A[i] == NULL) {
            printf("Memory allocation failed.\n");
            exit(1);
        }
    }

    for(int i=0;i<=n;i++) {
        for(int j=0;j<=n+1;j++) {
            A[i][j] = 0;
        }
    }

    for (int i = 0; i <= n; i++) {
        A[i][0] = x_now[i];
        A[i][1] = y_now[i];
    }

    for (int j = 2; j <= n+1; j++) {
        for (int i = j-1; i <= n; i++) {
            A[i][j] = (A[i][j-1] - A[i-1][j-1]) / (A[i][0] - A[i-(j-1)][0]);
        }
    }
    
    double result = A[0][1] + A[1][2] *(x - x_now[0]) + A[2][3] *(x - x_now[0])*(x - x_now[1]);    
    //释放内存
    for(int i=0;i<=n;i++) {
        free(A[i]);
    }
    free(A);
    
    return result;
}


int main() {

    int n = 4;
    double * x_now;
    double * y_now;
    x_now = (double*)malloc(n * sizeof(double));
    y_now = (double*)malloc(n * sizeof(double));

    x_now[0] = -2;
    x_now[1] = -1;
    x_now[2] = 1;
    x_now[3] = 2;

    y_now[0] = 5;
    y_now[1] = 3;
    y_now[2] = 17;
    y_now[3] = 21;

    double x = 1.5;

    double result_1 = newton_1(x, x_now, y_now, n);

    double result_2 = newton_2(x, x_now, y_now, n);

    printf("linear\t : %lf\n", result_1);
    printf("quadratic: %lf\n", result_2);

    return 0;
}