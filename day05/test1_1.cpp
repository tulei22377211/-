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

int main() {
    double a = 0;
    double b = 1;
    double eps_x = 1e-5;
    double eps_dy = 1e-4;
    double eps_y = 1e-4;

    int count = 1;
    int MAX_COUNT = 100;

    double min_x = a;
    double max_x = b;    

    double r = 0.6180340025899865;
    double r0 = 1 - r;

    double right_x = r* (max_x - min_x) + min_x;
    double left_x = r0* (max_x - min_x) + min_x;
    double result_x = right_x;

    

    printf("\n");
    printf("--------------------------\n");
    printf("iter\tx\t\tf(x)\t\tdf(x)\n");
    printf("%d\t%lf\t%lf\t%lf\n", count, result_x, f(result_x), df(result_x));
    while ((fabs(max_x - min_x) >= eps_x) && (fabs(df(max_x) - df(result_x)) + fabs(df(min_x) - df(result_x)) >= eps_y) && (count <= MAX_COUNT )) {
        if (f(right_x) > f(left_x)) {
            max_x = right_x;
            right_x = left_x;
            left_x = r0* (max_x - min_x) + min_x;
            result_x = left_x;
        } else {
            min_x = left_x;
            left_x = right_x;
            right_x = r* (max_x - min_x) + min_x;
            result_x = right_x;
        }
          
        printf("%d\t%lf\t%lf\t%lf\n", count+1, result_x, f(result_x), df(result_x));
        count++;
    }
    printf("--------------------------\n");
    if((count > MAX_COUNT) && (fabs(df(result_x)) >= eps_y) && (fabs(max_x - min_x) >= eps_x) && (count <= MAX_COUNT )) {
        printf("The result was not found within %d iterations.\n", MAX_COUNT);
    } else {
        printf("The result is %lf\n", result_x);
    }



    return 0;
}