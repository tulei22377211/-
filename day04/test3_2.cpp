#include <stdio.h>
#include <math.h>

double f(double x) {
    return x*x*x - x - 1;
}

double df(double x) {
    return 3*x*x - 1;
}

double newton(double x0) {
    return x0 - f(x0)/df(x0);
}

double gexian(double x0, double x1) {
    return x1 - f(x1)*(x1-x0)/(f(x1)-f(x0));
}

int main() {
    double a = 1;
    double b = 2;
    
    double eps_x = 1e-5;
    double eps_y = 1e-4;
    int iter = 0;
    int max_iter = 100;
    

    // 牛顿迭代法
    double x0 = 1.5;
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
    x0 = 1;
    double x1 = 2;
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