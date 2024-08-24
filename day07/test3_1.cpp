#include <iostream>
using namespace std;

#include "func.hpp"

double f(double t) {
    double x = 0.5 + 0.3*t + 3.9*t*t - 4.7*t*t*t;
    double y = 1.5 + 0.3*t + 0.9*t*t - 2.7*t*t*t;
    double dx = 0.3 + 3.9*2*t - 4.7*3*t*t;
    double dy = 0.3 + 0.9*2*t - 2.7*3*t*t;
    double result = dx*dx + dy*dy;
    return sqrt(result);
}

double f_x(double x) {
    return 0;
}
double f_y(double x) {
    return 0;
}

int main() {

    double a = 0.0, b = 1.0, n = 5;
    double length = B(a, b, n);
    cout << "The length of the curve is " << length << endl;

    return 0;
}