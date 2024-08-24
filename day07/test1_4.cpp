#include <iostream>
using namespace std;

#include "func.hpp"

double f(double x) {
    return 1.0/x;
}

double f_x(double x) {
    return 0;
}
double f_y(double x) {
    return 0;
}

int main(void) {
    double real = 1.252762968495367995688120;
    
    cout << fixed;
    cout.precision(8);

    cout << "C_T(f)" << endl;
    cout << "n\tM\th\t\treal\t\tC_T(f)\t\terror" << endl;
    cout << "----------------" << endl;
    for (int i = 3; i <= 10000; i += 2) {
        double result_C_T = C_T(2, 7, i);
        if (abs(err(real,result_C_T)) > 5.0e-9){
            cout << i << "\t" << i-1 << "\t" << 1.0/(i-1) << "\t" << real << "\t" << result_C_T << "\t" << err(real,result_C_T) << endl;
        } else {
            cout << i << "\t" << i-1 << "\t" << 1.0/(i-1) << "\t" << real << "\t" << result_C_T << "\t" << err(real,result_C_T) << endl;
            break;
        }
        if (i > 9998) {
            cout << "Error: too many iterations" << endl;
            break;
        }
    }
    

    cout << endl;

    cout << "C_S(f)" << endl;
    cout << "n\tM\th\t\treal\t\tC_S(f)\t\terror" << endl;
    cout << "----------------" << endl;
    for (int i = 3; i <= 7000; i += 2) {
        double result_C_S = C_S(2, 7, i);
        if (abs(err(real,result_C_S)) > 5.0e-9){
            cout << i << "\t" << (1/2)*(i-1) << "\t" << 1.0/(i-1) << "\t" << real << "\t" << result_C_S << "\t" << err(real,result_C_S) << endl;
        } else {
            cout << i << "\t" << (1/2)*(i-1) << "\t" << 1.0/(i-1) << "\t" << real << "\t" << result_C_S << "\t" << err(real,result_C_S) << endl;
            break;
        }
        if (i >6998) {
            cout << "Error: too many iterations" << endl;
            break;
        }

    }
    
    return 0;
}

      
