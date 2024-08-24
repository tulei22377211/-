#include <stdio.h>
#include <math.h>

double f(double x) {
    return x*x*x - x - 1;
}

double diedai_1(double x0){
    return x0*x0*x0 - 1;
}

double diedai_2(double x0){
    return diedai_1(x0)*diedai_1(x0)*diedai_1(x0) - 1;
}

double jiasu(double x0 ){
    return x0 - (diedai_1(x0) - x0)*(diedai_1(x0) - x0)/(x0 - 2*diedai_1(x0) +diedai_2(x0));
}

int main(void) {
    double a = 1;
    double b = 2;
    double x0 = 1.5;
    double eps_x = 1e-5;
    double eps_y = 1e-4;
    int max_iter = 1000;

    
    //简单迭代法
    double last_x = 999;
    int iter = 0;
    double x = x0;
    double y = f(x);
    printf("\n");
    printf("simple iteration method:\n");
    printf("--------------------------\n");
    printf("iter\tx\t\tf(x)\t\tx-f(x)\n");
    printf("%d\t%lf\t%lf\t%lf\n", iter, x, y, x-f(x));
    while ((fabs(y) >= eps_y) && (fabs(diedai_1(x)-x) >= eps_y) && (iter <= max_iter ) && (fabs(diedai_1(x)-x) >= eps_x)){
        last_x = x;
        x = diedai_1(x);
        y = f(x);
        printf("%d\t%lf\t%lf\t%lf\n", iter+1, x, y, x-f(x));
        iter++;
    }
    printf("--------------------------\n");
    if((iter > max_iter) && (fabs(y) >= eps_y) && (fabs(diedai_1(x)-x) >= eps_y) && (fabs(diedai_1(x)-x) >= eps_x)) {
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

    //加速迭代法
    last_x = 999;
    iter = 0;
    x = x0;
    y = f(x);
    printf("\n");
    printf("accelerated iteration method:\n");
    printf("--------------------------\n");
    printf("iter\tx\t\tf(x)\t\tx-f(x)\n");
    printf("%d\t%lf\t%lf\t%lf\n", iter, x, y, x-f(x));
    while ((fabs(diedai_1(x)-x) >= eps_x) && (fabs(jiasu(x)-x) >= eps_y) && (iter <= max_iter ) && (fabs(diedai_2(x)-2*diedai_1(x)+x) >= eps_x)){
        last_x = x;
        x = jiasu(x);
        y = f(x);
        printf("%d\t%lf\t%lf\t%lf\n", iter+1, x, y, x-f(x));
        iter++;
    }
    printf("--------------------------\n");
    if((iter > max_iter) && (fabs(diedai_1(x)-x) >= eps_x) && (fabs(jiasu(x)-x) >= eps_y) && (fabs(diedai_2(x)-2*diedai_1(x)+x) >= eps_x) ) {
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