#include <stdio.h>
#include <math.h>

#include "convert.hpp"
#include "error.hpp"

//第一种迭代方法
double func1(double I0,int n){
    double I1 = 1 - I0 * n;
    return I1;
}
//第二种迭代方法
double func2(double I1, int n){
    double I0 = (double)(1/(double)n) * (double)(1 - I1);
    return I0;
}

int main(void) {
    //真实值
    //printf("真实值：\n\n");
    double e = exp(1);
    double I0 = 1 - 1 / e;
    double I1 = 1/e;
    double I2 = 1 - 2/e;
    double I3 = -2 + 6/e;
    double I4 = 9 - 24/e;
    double I5 = -44 + 120/e;
    double I6 = 265 - 720/e;
    double I7 = -1854 + 5040/e;
    double I8 = 14833 - 40320/e;
    double I9 = -133496 + 362880/e;

    double I_0[10] = {I0, I1, I2, I3, I4, I5, I6, I7, I8, I9};
    for(int i = 0; i < 10; i++) {

        printf("I%d = %lf\n", i, I_0[i]);
    }
    
    //printf("\n");

//第一种迭代方法
    //printf("第一种迭代方法：\n\n");
    double I_1[10];
    I_1[0] = I0;
    
    for(int i = 0; i < 10; i++) {
        I_1[i+1] = con_d_4(func1(I_1[i], i+1));
    }


    //for(int i = 0; i < 10; i++) {
    //    printf("I%d = %.5lf\n其与真值的绝对误差为%.16lf,相对误差为：%.16lf\n", i, I_1[i], absolute_error(I_1[i], I_0[i]), relative_error(I_1[i], I_0[i]));
    //}
    //printf("\n");


//第二种迭代方法
    //printf("第二种迭代方法：\n\n");
    double I_2[10];
    I_2[9] = I9;

    for(int i = 9; i > 0; i--) {
        I_2[i-1] = con_d_4(func2(I_2[i], i));
    }

    //for(int i = 0; i < 10; i++) {
    //    printf("I%d = %.5lf\n其与真值的绝对误差为%.16lf,相对误差为：%.16lf\n", i, I_2[i], absolute_error(I_2[i], I_0[i]), relative_error(I_2[i], I_0[i]));
    //}
    
    printf("--------------------\n");
    printf("Index\trealvalue\tobserve_1\tabsolute_error_1\trelative_error_1\t\tobserve_2\tabsolute_error_2\trelative_error_2\n");
    for(int i = 0; i < 10; i++) {
        printf("I%d\t%.4lf\t\t%.4lf\t\t%.16lf\t%.16lf\t\t%.4lf\t\t%.16lf\t%.16lf\n",i, I_0[i], I_1[i], absolute_error(I_1[i], I_0[i]), relative_error(I_1[i], I_0[i]), I_2[i], absolute_error(I_2[i], I_0[i]), relative_error(I_2[i], I_0[i]));
    }

    return 0;
}