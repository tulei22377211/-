#include <iostream>
#include "func.hpp"

using namespace std;

double f(double x, double y) {
    return 0;
}

double x_t(double t, double x, double y) {
    return x + 2*y;
}

double y_t(double t, double x, double y) {
    return 3*x + 2*y;
}

int main() {

    double a = 0.0;
    double b = 0.2;
    double t0 = 0;
    double x0 = 6;
    double y0 = 4;
    double h = 0.02;

    cout << fixed;
    cout.precision(8);

    cout << "t\t\tx\t\ty" << endl;
    cout << "----------" << endl;
    for (int i = 0; i*h <= b; i++){
        double t = t0 + i*h;
        double x = RK_4_2_x(a, t, t0, x0, y0, h);
        double y = RK_4_2_y(a, t, t0, x0, y0, h);
        cout << t << "\t" << x << "\t" << y << endl;
    }

    return 0;
}