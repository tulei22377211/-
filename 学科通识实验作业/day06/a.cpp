#include <stdio.h>
#include <math.h>

int main() {

    double x1 = -1; double y1 = 0;
    double x2 = 1; double y2 = 0.5;
    double x3 = 1; double y3 = -0.5;
    double x4 = 0; double y4 = 1;

    double r1 = 1;
    double r2 = 0.5;
    double r3 = 0.5;
    double r4 = 0.5;

    double x = 5701178.766626;
    double y = 10549541.142106;
    double K = 11991506.766177;

    double result = sqrt((x-x1)*(x-x1) + (y-y1)*(y-y1)) + sqrt((x-x2)*(x-x2) + (y-y2)*(y-y2)) + sqrt((x-x3)*(x-x3) + (y-y3)*(y-y3)) + sqrt((x-x4)*(x-x4) + (y-y4)*(y-y4))-r1-r2-r3-r4-4*K;

    printf("%lf\n", result);
    return 0;
}