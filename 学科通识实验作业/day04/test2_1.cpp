#include <stdio.h>
#include <math.h>

double f(double x) {
    return x*x*x + 4 * x*x - 10;
}

double diedai_1(double x0){
    return 0.5 * sqrt(10 - x0*x0*x0);
}

double diedai_2(double x0){
    return  sqrt(10/x0 - 4*x0);
}

double diedai_3(double x0){
    return  x0 - x0*x0*x0 - 4*x0*x0 + 10;
}

int main(void) {
    double a = 1;
    double b = 2;
    double x0 = 1.5;
    double eps_x = 1e-5;
    double eps_y = 1e-4;
    int MAX_ITER = 100;

    
    //迭代公式一
    double last_x = 999;
    int iter = 0;
    double x = x0;
    double y = f(x);
    printf("simple iteration method_1:\n");
    printf("--------------------------\n");
    printf("iter\tx\t\tf(x)\t\tx-f(x)\n");
    printf("%d\t%lf\t%lf\t%lf\n", iter, x, y, x-f(x));
    while ((fabs(y) >= eps_y) && (fabs(diedai_1(x)-x) >= eps_y) && (iter <= MAX_ITER ) && (fabs(diedai_1(x)-x) >= eps_x)){
        last_x = x;
        x = diedai_1(x);
        y = f(x);
        printf("%d\t%lf\t%lf\t%lf\n", iter+1, x, y, x-f(x));
        iter++;
    }
    printf("--------------------------\n");
    if((iter > MAX_ITER) && (fabs(y) >= eps_y) && (fabs(diedai_1(x)-x) >= eps_y) && (fabs(diedai_1(x)-x) >= eps_x)) {
        printf("The result was not found within %d iterations.\n", MAX_ITER);
    } else if(isnan(x)){
        printf("The result was not found within %d iterations.\n", MAX_ITER);
    } else if(isnan(f(x))){
        printf("The result was not found within %d iterations.\n", MAX_ITER);
    } else if(x <a || x > b){
        printf("The result was not found within %d iterations.\n", MAX_ITER);
    } else {
        printf("root is %lf\n", x);
    }
    printf("\n");

    //迭代公式二
    last_x = 999;
    iter = 0;
    x = x0;
    y = f(x);
    printf("simple iteration method_2:\n");
    printf("--------------------------\n");
    printf("iter\tx\t\tf(x)\t\tx-f(x)\n");
    printf("%d\t%lf\t%lf\t%lf\n", iter, x, y, x-f(x));
    while ((fabs(y) >= eps_y) && (fabs(diedai_2(x)-x) >= eps_y) && (iter <= MAX_ITER ) && (fabs(diedai_2(x)-x) >= eps_x)){
        last_x = x;
        x = diedai_2(x);
        y = f(x);
        printf("%d\t%lf\t%lf\t%lf\n", iter+1, x, y, x-f(x));
        iter++;
    }
    printf("--------------------------\n");
    if((iter > MAX_ITER) && (fabs(y) >= eps_y) && (fabs(diedai_2(x)-x) >= eps_y) && (fabs(diedai_2(x)-x) >= eps_x)) {
        printf("The result was not found within %d iterations.\n", MAX_ITER);
    } else if(isnan(x)){
        printf("The result was not found within %d iterations.\n", MAX_ITER);
    } else if(isnan(f(x))){
        printf("The result was not found within %d iterations.\n", MAX_ITER);
    } else if(x <=a || x >= b){
        printf("The result was not found within %d iterations.\n", MAX_ITER);
    } else {
        printf("root is %lf\n", x);
    }
    printf("\n");

    //迭代公式三
    last_x = 999;
    iter = 0;
    x = x0;
    y = f(x);
    printf("simple iteration method_3:\n");
    printf("--------------------------\n");
    printf("iter\tx\t\tf(x)\t\tx-f(x)\n");
    printf("%d\t%lf\t%lf\t%lf\n", iter, x, y, x-f(x));
    while ((fabs(y) >= eps_y) && (fabs(diedai_3(x)-x) >= eps_y) && (iter <= MAX_ITER ) && (fabs(diedai_3(x)-x) >= eps_x)){
        last_x = x;
        x = diedai_3(x);
        y = f(x);
        printf("%d\t%lf\t%lf\t%lf\n", iter+1, x, y, x-f(x));
        iter++;
    }
    printf("--------------------------\n");
    if((iter > MAX_ITER) && (fabs(y) >= eps_y) && (fabs(diedai_3(x)-x) >= eps_y) && (fabs(diedai_3(x)-x) >= eps_x)) {
        printf("The result was not found within %d iterations.\n", MAX_ITER);
    } else if(isnan(x)){
        printf("The result was not found within %d iterations.\n", MAX_ITER);
    } else if(isnan(f(x))){
        printf("The result was not found within %d iterations.\n", MAX_ITER);
    } else if(x <=a || x >= b){
        printf("The result was not found within %d iterations.\n", MAX_ITER);
    } else if(fabs(f(x))> 0){
        printf("The result was not found within %d iterations.\n", MAX_ITER);
    } else {
        printf("root is %lf\n", x);
    }
    printf("\n");
    




    return 0;
}