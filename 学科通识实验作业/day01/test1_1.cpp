#include <stdio.h>

#include "convert.hpp"
#include "error.hpp"

int main(void) {

    double f1 = 0.12346845;
    double f2 = 0.012346845;
    double f3 = 0.0012346845;
    double f4 = 0.00012346845;

    double f1_4 = con_d_4(f1);
    double f2_4 = con_d_4(f2);
    double f3_4 = con_d_4(f3);
    double f4_4 = con_d_4(f4);

    printf("real\tobservation\tabsolute_error\trelative_error\n");
    printf("%f\t%f\t%.8f\t%.8f\n",f1, f1_4, absolute_error(f1_4, f1), relative_error(f1_4, f1));
    printf("%f\t%f\t%.8f\t%.8f\n",f2, f2_4, absolute_error(f2_4, f2), relative_error(f2_4, f2));
    printf("%f\t%f\t%.8f\t%.8f\n",f3, f3_4, absolute_error(f3_4, f3), relative_error(f3_4, f3));
    printf("%f\t%f\t%.8f\t%.8f\n",f4, f4_4, absolute_error(f4_4, f4), relative_error(f4_4, f4));

    //printf("%f 的四位有效数字近似值的绝对误差为 %.8f,相对误差为 %.8f\n", f1, absolute_error(f1_4, f1), relative_error(f1_4, f1));
    //printf("%f 的四位有效数字近似值的绝对误差为 %.8f,相对误差为 %.8f\n", f2, absolute_error(f2_4, f2), relative_error(f2_4, f2));
    //printf("%f 的四位有效数字近似值的绝对误差为 %.8f,相对误差为 %.8f\n", f3, absolute_error(f3_4, f3), relative_error(f3_4, f3));
    //printf("%f 的四位有效数字近似值的绝对误差为 %.8f,相对误差为 %.8f\n", f4, absolute_error(f4_4, f4), relative_error(f4_4, f4));
    return 0;
}