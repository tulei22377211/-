#include <iostream>
#include "func.hpp"


using namespace std;

double f(double x, double y){
    return (x-y)/2.0;
}

double x_t(double x, double y, double t){
    return 0;
} //定义微分方程的x项
double y_t(double x, double y, double t){
    return 0;
} //定义微分方程的y项

int main(){

    
    double t0 = 0.0;
    double y0 = 1.0;
    double h = 1/100.0;
    double a = 0;
    double b = 3;//求解y(b)

    cout << fixed;
    cout.precision(6);

    cout << "Runge-Kutta 4" << endl;
    cout << "---------------" << endl;
    cout << "h\t\ty(3)\t\treal_y(3)\terror" <<endl;
    cout << "---------------" << endl;
    
    for (double h = 1.0; h >= 1.0/128.0; h /= 2.0){
        cout << "1/" << (int)(1/h) << "\t\t" << RK_4(a, b, t0, y0, h) << "\t" << 3.0*exp(-3/2.0)-2+3.0<< "\t" << 3.0*exp(-3/2.0)-2+3.0 - RK_4(a, b, t0, y0, h) << endl;
    }
    cout << "---------------" << endl;
    
    cout << endl;

    cout << "When h = 1/64, the y(t) is: " << endl;
    h = 1.0/64;
    cout << "---------------" << endl;
    cout << "t\t\ty(t)\t\treal_y(t)\terror" << endl;
    cout << "---------------" << endl;
    for(int i = 0; i*h <= b; i++){
        double t = t0 + i*h;
        if (t == 0.0 || t == 0.125 || t == 0.25 || t == 0.375 || t == 0.5 || t == 0.75 ||t == 1.0 || t == 1.5 || t == 2.0 || t == 2.5 || t == 3.0){
            double y = RK_4(a, t, t0, y0, h);
        cout << t << "\t" << y << "\t" << 3.0*exp(-t/2.0)-2+t << "\t" << 3.0*exp(-t/2.0)-2+t - y << endl;
        }
    }
    
    return 0;
}