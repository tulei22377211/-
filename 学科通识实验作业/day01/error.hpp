#include <stdio.h>

//计算相对误差
long double relative_error(double measure, double real) {
    if(real != 0) {
        return (long double)(real - measure) / (long double)real;
    } else {
        return 99999;
    }
}

//计算绝对误差
long double absolute_error(double measure, double real) {
    return (long double)(real - measure);
}