#include <stdio.h>
#include <math.h>

double f(double x) {
    return x*sin(x) - 1;
}

int main(void) {
    double min_x = 0;
    double max_x = 2;
    double eps_x = 1e-5;
    double eps_y = 1e-4;

    int count = 1;
    int MAX_COUNT = 100;

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
    } else {
        printf("root is %lf\n", x);
    }
    return 0;
}