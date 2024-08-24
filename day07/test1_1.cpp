#include <iostream>
using namespace std;

#include "func.hpp"

double f(double x) {
    return 1 + exp(-x) * sin(4*x);
}

double f_x(double x) {
    return 0;
}
double f_y(double x) {
    return 0;
}

int main(void) {
    double real = 1.3082506046426;
    double result_T = T(0,1,2);
    double result_S = S(0,1,3);
    double result_S_3_8 = S_3_8(0,1,4);
    double result_B = B(0,1,5);

    cout << fixed;
    cout.precision(8);

    cout << "\t\treal_result\tT(f)\t\tS(f)\t\tS_3_8(f)\t\tB(f)" << endl;
    cout << "----------------" << endl;
    cout <<"result\t\t" << real << "\t" << result_T << "\t" << result_S << "\t" << result_S_3_8 << "\t" << result_B << endl;
    cout <<"error\t\t"<< 0 << "\t\t" << err(real,result_T) << "\t" << err(real,result_S) << "\t" << err(real,result_S_3_8) << "\t" << err(real,result_B) << endl;
    cout << "----------------" << endl;
    return 0;
}