//割线法为 当不知函数的导数 而无法使用牛顿法时的平替
//此方程已知函数形式 且导数易得 故可使用牛顿法求根
//为防止牛顿法在极值点附近无限次徘徊的情况
//在外边包一层二分法

#include <stdio.h>
#include <math.h>

const double pi = 3.14159265358979323846264338327950288419716939937510;

double f(double x) {
    return x + 4.0 * pi * sin(pi*x);
}

double df(double x) {
    return 1.0 + 4.0 * pi * pi * cos(pi*x);
}

double newton(double x0) {
    return x0 - f(x0)/df(x0);
}

double mid(double x0, double x1) {
    return (x0 + x1) / 2.0;
}

int main(void) {

    double a = -3;
    double b = 4;
    
    double eps_x = 1e-5;
    double eps_y = 1e-5;
    int iter = 0;
    int max_iter_mid = 60;
    int max_iter_newton = 40;
    int max_iter = max_iter_mid + max_iter_newton;

    // 牛顿迭代法
    double x0 = 4;
    double x = x0;
    double y = f(x0);
    double min_x = a;
    double max_x = b;
    
    printf("\n");
    printf("The first 40 steps are Newton's method, and the last 60 steps are bisection method.:\n");
    printf("--------------------------\n");
    printf("iter\tx\t\tf(x)\t\tx-f(x)\n");
    printf("%d\t%lf\t%lf\t%lf\n", iter, x, y, x-f(x));

    
    int result = 0;
    
    while ((fabs(y) >= eps_y) && (fabs(newton(x)-x) >= eps_x) && (fabs(df(x))>=eps_y)  ){
        x = newton(x);
        y = f(x);
        printf("%d\t%lf\t%lf\t%lf\n", iter+1, x, y, x-f(x));
        iter++;
        if (iter+1 > max_iter_newton) {
            result = 0;
            break;
        } else if ((fabs(y) < eps_y) || (fabs(newton(x)-x) < eps_x) || (fabs(df(x))<eps_y) ){
            result = 1;
            if(x <=a || x >= b) result = 0;
            break;
        }
    }
    if (result == 0) printf("Use dichotomy from here\n");
    // 二分法
    if (result == 0) {
        max_x = b;
        min_x = a;
        x0 = -3;
        double x1 = 4;
        
        
        x = mid(x0, x1);
        y = f(x);
        
        printf("%d\t%lf\t%lf\t%lf\n", iter+1, x, y, x-f(x));
        while ((fabs(y) >= eps_y) && (fabs(x-min_x) >= eps_x) ) {
            if (f(x) * f(min_x) < 0) {
                max_x = x;
            } else {
                min_x = x;
            }
            x = (min_x + max_x) / 2;
            y = f(x);        
            printf("%d\t%lf\t%lf\t%lf\n", iter+1, x, y, x-f(x));
            iter++;
            if(iter-1 > max_iter_mid) {
                result = 0;
                break;
            } else if ((fabs(y) < eps_y) || (fabs(x-min_x) < eps_x) ){
                result = 1;
                break;
            }
        } 

    }
    
    

    
    

    printf("--------------------------\n");
    if((result == 0) ) {
        printf("The result was not found within %d iterations.\n", max_iter);
    } else if(isnan(x)){
        printf("The result was not found within %d iterations.\n", max_iter);
    } else if(isnan(f(x))){
        printf("The result was not found within %d iterations.\n", max_iter);
    }  else if(x <=a || x >= b){
        printf("The result was not found within %d iterations.\n", max_iter);
    }  else {
        printf("root is %lf\n", x);
        printf("x = %lf f(x) = %lf\n", x, f(x));
        printf("\n");
    }
    
    


    return 0;
}
