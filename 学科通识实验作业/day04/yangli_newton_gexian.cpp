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

double newton(double x0) {
    return x0 - f(x0)/df(x0);
}

double gexian(double x0, double x1) {
    return x1 - f(x1)*(x1-x0)/(f(x1)-f(x0));
}

int main(void) {

    double a = -0.2;
    double b = 6.1;
    
    double eps_x = 1e-5;
    double eps_y = 1e-4;
    int iter = 0;
    int max_iter = 100;
    

    // 牛顿迭代法
    double x0 = 4;
    double x = x0;
    double y = f(x0);
    
    printf("\n");
    printf("Newton's iteration method:\n");
    printf("--------------------------\n");
    printf("iter\tx\t\tf(x)\t\tx-f(x)\n");
    printf("%d\t%lf\t%lf\t%lf\n", iter, x, y, x-f(x));
    while ((fabs(y) >= eps_y) && (fabs(newton(x)-x) >= eps_x) && (fabs(df(x))>=eps_y) && (iter <= max_iter ) ){
        
        x = newton(x);
        y = f(x);
        printf("%d\t%lf\t%lf\t%lf\n", iter+1, x, y, x-f(x));
        iter++;
    }
    printf("--------------------------\n");
    if((iter > max_iter) && (fabs(y) >= eps_y) && (fabs(newton(x)-x) >= eps_x) && (fabs(df(x))>=eps_y)) {
        printf("The result was not found within %d iterations.\n", max_iter);
    } else if(isnan(x)){
        printf("The result was not found within %d iterations.\n", max_iter);
    } else if(isnan(f(x))){
        printf("The result was not found within %d iterations.\n", max_iter);
    }  else if(x <=a || x >= b){
        printf("The result was not found within %d iterations.\n", max_iter);
    }  else {
        printf("root is %lf\n", x);
    }
    printf("\n");

    //割线法
    x0 = 4;
    double x1 = 4.5;
    x = x1;
    y = f(x1);
    iter = 0;
    
    printf("\n");
    printf("secant method:\n");
    printf("--------------------------\n");
    printf("iter\tx\t\tf(x)\t\tx-f(x)\n");
    printf("%d\t%lf\t%lf\t%lf\n", iter, x, y, x-f(x));
    while ((fabs(y) >= eps_y) && (fabs(gexian(x0,x1)-x) >= eps_x) && (fabs((gexian(x0,x1)-x1)*(gexian(x0,x1)-x0))>=eps_x) && (iter <= max_iter ) ){
        
        x = newton(x);
        y = f(x);
        printf("%d\t%lf\t%lf\t%lf\n", iter+1, x, y, x-f(x));
        iter++;
    }
    printf("--------------------------\n");
    if((iter > max_iter) && (fabs(y) >= eps_y) && (fabs(gexian(x0,x1)-x) >= eps_x) && (fabs((gexian(x0,x1)-x1)*(gexian(x0,x1)-x0))>=eps_x)) {
        printf("The result was not found within %d iterations.\n", max_iter);
    } else if(isnan(x)){
        printf("The result was not found within %d iterations.\n", max_iter);
    } else if(isnan(f(x))){
        printf("The result was not found within %d iterations.\n", max_iter);
    }  else if(x <=a || x >= b){
        printf("The result was not found within %d iterations.\n", max_iter);
    }  else {
        printf("root is %lf\n", x);
    }
    printf("\n");

    
    return 0;
}