#include <stdio.h>

#include "convert.hpp"
#include "error.hpp"

int main(void) {

    double pi0 = 3.141592653589793;
    double pi1 = pi0 + 0.000123456789123;

    double ppi0 = con_d_8(pi0);
    double ppi1 = con_d_8(pi1);
    printf("Approximate to eight significant figures:\n");
    printf("real\tobservation\tabsolute_error\trelative_error\n");

    printf("%lf\t%lf\t%.16lf%\t%.16lf\n", pi0, ppi0, absolute_error(ppi0, pi0), relative_error(ppi0, pi0));
    printf("%lf\t%lf\t%.16lf%\t%.16lf\n", pi1, ppi1, absolute_error(ppi1, pi1), relative_error(ppi1, pi1));

    //printf("%lf 做八位有效数字近似的绝对误差为 %.16lf, 相对误差为 %.16lf\n", pi0, absolute_error(ppi0, pi0), relative_error(ppi0, pi0));
    //printf("%lf 做八位有效数字近似的绝对误差为 %.16lf, 相对误差为 %.16lf\n", pi1, absolute_error(ppi1, pi1), relative_error(ppi1, pi1));

    double mp = con_d_8(ppi0 - ppi1);
    printf("\n\n");
    printf("The difference between the two approximations\n");
    printf("real\tobservation\tabsolute_error\trelative_error\n");
    printf("%lf\t%lf\t%.16lf%\t%.16lf\n", 0.000123456789123, mp, absolute_error(mp, 0.000123456789123), relative_error(mp, 0.000123456789123));
    //printf("近似差值%lf 的绝对误差为 %.16lf, 相对误差为 %.16lf\n", mp, absolute_error(mp, 0.000123456789123), relative_error(mp, 0.000123456789123));
    return 0;
}