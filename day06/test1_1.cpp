#include <stdio.h>
#include <stdlib.h>

double lagrange_1(double x, double * x_now, double * y_now, int n) {
    //假设x是递增的,并且x不与已知的x_now重复,且x在x_now的范围内
    int count = 0;

    for (int i = 0; i < n; i++) {
        if (x>x_now[i]) {
            count++;
        }
    }

    double left_x = x_now[count-1];
    double right_x = x_now[count];
    double left_y = y_now[count-1];
    double right_y = y_now[count];

    double result = (x - left_x) * right_y / (right_x - left_x) + (right_x - x) * left_y / (right_x - left_x);

    return result;
}

double lagrange_2(double x, double * x_now, double * y_now, int n) {
    //假设x是递增的,并且x不与已知的x_now重复,且x在x_now的范围内
    int count = 0;

    for (int i = 0; i < n; i++) {
        if (x>x_now[i]) {
            count++;
        }
    }

    double x_left = x_now[count-1];  
    double x_mid = x_now[count];
    double x_right = x_now[count+1];

    double y_left = y_now[count-1];  
    double y_mid = y_now[count];
    double y_right = y_now[count+1];

    double result = y_left * ((x - x_mid)* (x - x_right)) / ((x_left - x_mid) * (x_left - x_right)) + y_mid * ((x - x_left) * (x - x_right)) / ((x_mid - x_left) * (x_mid - x_right)) + y_right * ((x - x_left) * (x - x_mid)) / ((x_right - x_left) * (x_right - x_mid));

    return result;
}


int main() {

    int n = 3;
    double * x_now;
    double * y_now;
    x_now = (double*)malloc(n * sizeof(double));
    y_now = (double*)malloc(n * sizeof(double));

    x_now[0] = 100;
    x_now[1] = 121;
    x_now[2] = 144;

    y_now[0] = 10;
    y_now[1] = 11;
    y_now[2] = 12;

    double x = 115;

    double result_1 = lagrange_1(x, x_now, y_now, n);

    double result_2 = lagrange_2(x, x_now, y_now, n);

    printf("linear\t : %lf\n", result_1);
    printf("quadratic: %lf\n", result_2);

    return 0;
}