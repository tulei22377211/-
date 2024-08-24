#include <stdio.h>
#include <math.h>

double f(double x) {
    return x*x -sin(x);
}

double df(double x) {
    return 2*x - cos(x);
}

double d2f(double x) {
    return 2 + sin(x);
}


    

double func(double x1, double x2, double x3) {
    //生成三点生成二次函数的极值点
    double y1 = f(x1);
    double y2 = f(x2);
    double y3 = f(x3);

    return 0.5*((x2*x2 - x3*x3)*y1 + (x3*x3 - x1*x1)*y2 + (x1*x1 - x2*x2)*y3)  /  ((x2 - x3)*y1 + (x3 - x1)*y2 + (x1 - x2)*y3);
}


int main() {
    double left = 0;
    double right = 1;
    double eps_x = 1e-5;
    double eps_dy = 1e-4;
    double eps_y = 1e-4;

    int count = 1;
    int MAX_COUNT = 100;

    double min_x = left;
    double max_x = right;   
    double last_x = 0.5*(max_x - min_x) + min_x;

    double new_x = func(min_x,last_x,max_x);
    double result_x = new_x;

    

    printf("\n");
    printf("--------------------------\n");
    printf("iter\tx\t\tf(x)\t\tdf(x)\n");
    printf("%d\t%lf\t%lf\t%lf\n", count, result_x, f(result_x), df(result_x));
    while ((fabs(max_x - min_x) >= eps_x) && (fabs(df(max_x) - df(result_x)) + fabs(df(min_x) - df(result_x)) >= eps_y) && (count <= MAX_COUNT ) &&(fabs(new_x - last_x) >= eps_x)) {
        if (new_x > last_x) {
            if (f(new_x) > f(last_x)) {
                max_x = new_x;
                last_x = last_x;
                new_x = func(min_x,last_x,max_x);
                result_x = new_x;
            } else {
                min_x = last_x;
                last_x = new_x;
                new_x = func(min_x,last_x,max_x);
                result_x = new_x;
            }
            printf("%d\t%lf\t%lf\t%lf\n", count+1, result_x, f(result_x), df(result_x));
            count++;
        } else {
            if (f(last_x) > f(new_x)) {
                max_x = last_x;
                last_x = new_x;
                new_x = func(min_x,last_x,max_x);
                result_x = new_x;
            } else {
                min_x = new_x;
                last_x = last_x;
                new_x = func(min_x,last_x,max_x);
                result_x = new_x;
            }
            printf("%d\t%lf\t%lf\t%lf\n", count+1, result_x, f(result_x), df(result_x));
            count++;
        }
    }
        
          
        printf("%d\t%lf\t%lf\t%lf\n", count+1, result_x, f(result_x), df(result_x));
        count++;
    
    printf("--------------------------\n");
    if((count > MAX_COUNT) && (fabs(df(result_x)) >= eps_y) && (fabs(max_x - min_x) >= eps_x) && (count <= MAX_COUNT )&&(fabs(new_x - last_x) >= eps_x)) {
        printf("The result was not found within %d iterations.\n", MAX_COUNT);
    } else {
        printf("The result is %lf\n", result_x);
    }




    return 0;
}