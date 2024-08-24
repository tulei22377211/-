#include <stdio.h>
#include <math.h>

double f(double x) {
    return x*sin(x) - 1;
}

double x_step(double a, double b){
    return b - f(b)*(b-a)/(f(b)-f(a));
}

int main(void) {
    double min_x = 0;
    double max_x = 2;
    double eps_x = 1e-5;
    double eps_y = 1e-4;

    int count = 1;
    int MAX_COUNT = 100;
    // trial value method
    double x = x_step(min_x, max_x);
    double y = f(x);
    printf("dichotomy\n");
    printf("--------------------------\n");
    printf("iter\tx\t\tf(x)\t\tx-f(x)\n");
    printf("%d\t%lf\t%lf\t%lf\n", count, x, y, x-f(x));
    int last_x = min_x - max_x;//保证初始时，上一个x不在区间内
    while ((fabs(y) >= eps_y) && (fabs(x-last_x) >= eps_x) && (count <= MAX_COUNT )) {
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
    if((count > MAX_COUNT) && (fabs(y) >= eps_y) && (fabs(x-min_x) >= eps_x)) {
        printf("The result was not found within %d iterations.\n", MAX_COUNT);
    } else {
        printf("root is %lf\n", x);
    }
    printf("\n");
    return 0;
}