#include <stdio.h>
#include <stdlib.h>

double f(double x) {

    int n = 7;
    double * x_now;
    double * y_now;
    x_now = (double*)malloc(n * sizeof(double));
    y_now = (double*)malloc(n * sizeof(double));

    x_now[0] = -0.1;
    x_now[1] = 0;
    x_now[2] = 1;
    x_now[3] = 4;
    x_now[4] = 5;
    x_now[5] = 5;
    x_now[6] = 6;

    y_now[0] = 0;
    y_now[1] = -8;
    y_now[2] = 0;
    y_now[3] = 6;
    y_now[4] = 1;
    y_now[5] = 1;
    y_now[6] = 4;

    
    double ** A = (double**)malloc(n * sizeof(double*));
    if (A == NULL) {
        printf("Memory allocation failed.\n");
        exit(1);
    }
    for(int i=0;i<n;i++) {
        A[i] = (double*)malloc((n+1) * sizeof(double));
        if (A[i] == NULL) {
            printf("Memory allocation failed.\n");
            exit(1);
        }
    }


    for(int i=0;i<n;i++) {
        for(int j=0;j<n+1;j++) {
            A[i][j] = 0;
        }
    }

    

    for (int i = 0; i < n; i++) {
        A[i][0] = x_now[i];
        A[i][1] = y_now[i];
    }

    for (int j = 2; j < n+1; j++) {
        for (int i = j-1; i < n; i++) {
            
            A[i][j] = (A[i][j-1] - A[i-1][j-1]) / (A[i][0] - A[i-(j-1)][0]);
            
        }
    }

    
    
    double result = A[0][1];
    for (int i = 1; i < n; i++) {
        double sum = 1;
        for (int k = 0; k < i;k++){
            sum = sum * (x - A[k][0]);
        }
        result = result + A[i][i+1]*sum; 
    }


    //释放内存
    for(int i=0;i<n;i++) {  
        free(A[i]);
    }
    free(A);
    free(x_now);
    free(y_now);

    return result;
}

double df(double x) {
    double eps = 0.000001;
    return (f(x+eps) - f(x-eps))/(2*eps);
}

double ddf(double x) {
    double eps = 0.000001;
    return (df(x+eps) - df(x-eps))/(2*eps);
}


int main() {
    double x1 = -0.1;
    double x2 = 0;
    double x3 = 1;
    double x4 = 4;
    double x5 = 5;
    double x6 = 5;
    double x7 = 6;

    

    double y1 = f(x1);
    double y2 = f(x2);
    double y3 = f(x3);  
    double y4 = f(x4);
    double y5 = f(x5);
    double y6 = f(x6);
    double y7 = f(x7);


    printf("x1 = %lf\ty1 = %lf\n", x1, y1);
    printf("x2 = %lf\ty2 = %lf\n", x2, y2);
    printf("x3 = %lf\ty3 = %lf\n", x3, y3);  
    printf("x4 = %lf\ty4 = %lf\n", x4, y4);
    printf("x5 = %lf\ty5 = %lf\n", x5, y5);
    printf("x6 = %lf\ty6 = %lf\n", x6, y6);
    printf("x7 = %lf\ty7 = %lf\n", x7, y7); 

    printf("\n");

    double eps = 0.000001;
    double x = 5;
    double dy = (f(x+eps) - f(x-eps))/(2*eps);
    printf("dy = %lf\n", dy);

    return 0;
}