#include <stdio.h>
#include <stdlib.h>
#include <math.h>


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
            if (i == 5 && j == 2){
                A[i][j] = 0;
            } else{
                A[i][j] = (A[i][j-1] - A[i-1][j-1]) / (A[i][0] - A[i-(j-1)][0]);
            }
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

double x_step(double a, double b){
    return b - f(b)*(b-a)/(f(b)-f(a));
}

int main() {

    double a = 0;
    double b = 6;
    double min_x = a;
    double max_x = b;
    double eps_x = 1e-5;
    double eps_y = 1e-4;

    int count = 1;
    int MAX_COUNT = 100;
    // dichotomy method
    double x = (min_x + max_x) / 2;
    double y = f(x);
    printf("dichotomy\n");
    printf("--------------------------\n");
    printf("iter\tx\t\tf(x)\t\tx-f(x)\n");
    printf("%d\t%lf\t%lf\t%lf\n", count, x, y, x-f(x));
    while ((fabs(y) >= eps_y) && (fabs(x-min_x) >= eps_x) && (count <= MAX_COUNT )) {
        if (f(x) * f(min_x) < 0) {
            max_x = x;
        } else {
            min_x = x;
        }
        x = (min_x + max_x) / 2;
        y = f(x);        
        printf("%d\t%lf\t%lf\t%lf\n", count+1, x, y, x-f(x));
        count++;
    }
    printf("--------------------------\n");
    if((count > MAX_COUNT) && (fabs(y) >= eps_y) && (fabs(x-min_x) >= eps_x)) {
        printf("The result was not found within %d iterations.\n", MAX_COUNT);
    } else if(isnan(x)){
        printf("The result was not found within %d iterations.\n", MAX_COUNT);
    } else if(isnan(f(x))){
        printf("The result was not found within %d iterations.\n", MAX_COUNT);
    } else if(x <=a || x >= b){
        printf("The result was not found within %d iterations.\n", MAX_COUNT);
    } else {
        printf("root is %lf\n", x);
    }
    printf("\n");

    //trial value method
    printf("trial value method\n");
    min_x = a;
    max_x = b;
    count = 1;
    x = x_step(min_x, max_x);
    y = f(x);
    printf("dichotomy\n");
    printf("--------------------------\n");
    printf("iter\tx\t\tf(x)\t\tx-f(x)\n");
    printf("%d\t%lf\t%lf\t%lf\n", count, x, y, x-f(x));
    int last_x = min_x - max_x;//保证初始时，上一个x不在区间内
    while ((fabs(y) >= eps_y) && (fabs(max_x-min_x) >= eps_x) && (count <= MAX_COUNT ) && (fabs((x_step(min_x, max_x)-min_x)*(x_step(min_x, max_x)-max_x))>=eps_x)) {
        if (f(x) * f(min_x) < 0) {
            max_x = x;
        } else {
            min_x = x;
        }
        last_x = x;
        x = x_step(min_x, max_x);
        y = f(x);        
        printf("%d\t%lf\t%lf\t%lf\n", count+1, x, y, x-f(x));
        count++;
    }
    printf("--------------------------\n");
    if((count > MAX_COUNT) && (fabs(y) >= eps_y) && (fabs(max_x-min_x) >= eps_x)  && (fabs((x_step(min_x, max_x)-min_x)*(x_step(min_x, max_x)-max_x))>=eps_x)) {
        printf("The result was not found within %d iterations.\n", MAX_COUNT);
    } else if(isnan(x)){
        printf("The result was not found within %d iterations.\n", MAX_COUNT);
    } else if(isnan(f(x))){
        printf("The result was not found within %d iterations.\n", MAX_COUNT);
    } else if(x <=a || x >= b){
        printf("The result was not found within %d iterations.\n", MAX_COUNT);
    } else {
        printf("root is %lf\n", x);
    }
    printf("\n");
    

    return 0;
}