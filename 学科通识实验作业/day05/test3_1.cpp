#include <stdio.h>
#include <math.h>

double f(double x) {
    return x*x -sin(x);
}

double df(double x) {
    return 2*x - cos(x);
}

double d2f(double x) {
    return 2 + sin(x);
}

double newton(double x0) {
    return x0 - df(x0)/d2f(x0);
}


    

double func(double x1, double x2, double x3) {
    //生成三点生成二次函数的极值点
    double y1 = f(x1);
    double y2 = f(x2);
    double y3 = f(x3);

    return 0.5*((x2*x2 - x3*x3)*y1 + (x3*x3 - x1*x1)*y2 + (x1*x1 - x2*x2)*y3)  /  ((x2 - x3)*y1 + (x3 - x1)*y2 + (x1 - x2)*y3);
}


int main() {
    double left = 0;
    double right = 1;
    double eps_x = 1e-5;
    double eps_dy = 1e-4;
    double eps_y = 1e-4;

    int count = 1;
    int MAX_COUNT = 100;

    // 牛顿迭代法
    double x0 = 1;
    double x = x0;
    double y = df(x0);
    
    printf("\n");
    printf("Newton's iteration method:\n");
    printf("--------------------------\n");
    printf("iter\tx\t\tf(x)\t\tdf(x)\n");
    printf("%d\t%lf\t%lf\t%lf\n", count, x, f(x), df(x));
    while ((fabs(y) >= eps_y) && (fabs(newton(x)-x) >= eps_x) && (fabs(df(x))>=eps_y) && (count <= MAX_COUNT ) ){
        
        x = newton(x);
        y = df(x);
        printf("%d\t%lf\t%lf\t%lf\n", count+1, x, f(x), df(x));
        count++;
    }
    printf("--------------------------\n");
    if((count > MAX_COUNT) && (fabs(y) >= eps_y) && (fabs(newton(x)-x) >= eps_x) && (fabs(df(x))>=eps_y)) {
        printf("The result was not found within %d iterations.\n", MAX_COUNT);
    } else if(isnan(x)){
        printf("The result was not found within %d iterations.\n", MAX_COUNT);
    } else if(isnan(f(x))){
        printf("The result was not found within %d iterations.\n", MAX_COUNT);
    }  else if(x <=left || x >= right){
        printf("The result was not found within %d iterations.\n", MAX_COUNT);
    }  else {
        printf("root is %lf\n", x);
    }
    printf("\n");




    return 0;
}