#include <iostream>
using namespace std;

#include "func.hpp"

double f(double x) {
    return 2 + sin(2*sqrt(x));
}

double f_x(double x) {
    return 0;
}
double f_y(double x) {
    return 0;
}

int main(void) {
    double real = 2*6 - sqrt(6)*cos(2*sqrt(6)) + (1.0/2) * sin(2*sqrt(6)) -(2*1 - sqrt(1)*cos(2 * sqrt(1)) + (1.0/2) * sin(2*sqrt(1)));
    double result_C_T = C_T(1,6,11);
    double result_C_S = C_S(1,6,11);

    cout << fixed;
    cout.precision(8);


    cout << "\t\treal_result\tC_T(f)\t\tC_S(f)" << endl;
    cout << "----------------" << endl;
    cout <<"result\t\t" << real << "\t" << result_C_T << "\t" << result_C_S << endl;
    cout <<"error\t\t"<< 0 << "\t\t" << err(real,result_C_T) << "\t" << err(real,result_C_S) << endl;
    cout << "----------------" << endl;
    return 0;
}