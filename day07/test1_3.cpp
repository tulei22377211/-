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
    double result_C_T_11 = C_T(1,6,11);
    double result_C_T_21 = C_T(1,6,21);
    double result_C_T_41 = C_T(1,6,41);
    double result_C_T_81 = C_T(1,6,81);
    double result_C_T_161 = C_T(1,6,161);

    double result_C_S_11 = C_S(1,6,11);
    double result_C_S_21 = C_S(1,6,21);
    double result_C_S_41 = C_S(1,6,41);
    double result_C_S_81 = C_S(1,6,81);
    double result_C_S_161 = C_S(1,6,161);

    cout << fixed;
    cout.precision(8);

    cout << "C_T(f)" << endl;
    cout << "----------------" << endl;
    cout << "n\tC_T(f)\t\terror" << endl;
    cout << "11\t" << result_C_T_11 << "\t" << err(real,result_C_T_11) << endl;
    cout << "21\t" << result_C_T_21 << "\t" << err(real,result_C_T_21) << endl;
    cout << "41\t" << result_C_T_41 << "\t" << err(real,result_C_T_41) << endl;
    cout << "81\t" << result_C_T_81 << "\t" << err(real,result_C_T_81) << endl;
    cout << "161\t" << result_C_T_161 << "\t" << err(real,result_C_T_161) << endl;
    cout << "----------------" << endl;

    cout << endl;

    cout << "C_S(f)" << endl;
    cout << "----------------" << endl;
    cout << "n\tC_S(f)\t\terror" << endl;
    cout << "11\t" << result_C_S_11 << "\t" << err(real,result_C_S_11) << endl;
    cout << "21\t" << result_C_S_21 << "\t" << err(real,result_C_S_21) << endl;
    cout << "41\t" << result_C_S_41 << "\t" << err(real,result_C_S_41) << endl;
    cout << "81\t" << result_C_S_81 << "\t" << err(real,result_C_S_81) << endl;
    cout << "161\t" << result_C_S_161 << "\t" << err(real,result_C_S_161) << endl;
    cout << "----------------" << endl;
    return 0;
}