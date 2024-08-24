#include <stdio.h>

#include "convert.hpp"
#include "error.hpp"

int main(void) {
    float pi20 = 7654321.0f;
    float pi21 = 0.2f;

    float result1 = pi20;
    for(int i = 0; i < 5; i++){
        result1 += pi21;
    }

    float result2 = pi20;
    float cd = pi21;
    for (int i = 0; i < 4; i++){
        cd += pi21;
    }
    result2 += cd;

    //printf("The sequential sum of %f and five %fs results in %f\n", pi20, pi21, result1);
    //printf("The sum of %f and five %fs results in %f\n", pi20, pi21, result2);
    printf("We believe that after summing the five decimals first, then adding the large number is the true value\n\n");
    printf("real\tobservation\tabsolute_error\trelative_error\n");
    printf("%f\t%f\t\t%.8lf\t\t%.8lf\n", result2,result1,absolute_error(result1,result2),relative_error(result1,result2));
    //printf("绝对误差为%.8lf,相对误差为%.8lf\n", absolute_error(result1,result2), relative_error(result1,result2));

    return 0;
}