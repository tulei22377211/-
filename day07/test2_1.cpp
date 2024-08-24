#include <iostream>
using namespace std;

#include "func.hpp"

double f(double x) {
    return 1/(x+2);
}

double f_x(double x) {
    return 0;
}
double f_y(double x) {
    return 0;
}

int main(void) {
    double real = 1.09861228866810;
    double result_T = T(-1,1,2);
    double result_S = S(-1,1,3);
    double result_G_2 = G_2(-1,1);

    cout << fixed;
    cout.precision(8);

  
    cout << "\t\treal_result\tT(f)\t\tS(f)\t\tG_2(f)" << endl;
    cout << "----------------" << endl;
    cout <<"result\t\t" << real << "\t" << result_T << "\t" << result_S << "\t" << result_G_2 << endl;
    cout <<"error\t\t"<< 0 << "\t\t" << err(real,result_T) << "\t" << err(real,result_S) << "\t" << err(real,result_G_2) << endl;
    cout << "----------------" << endl;
    return 0;
}