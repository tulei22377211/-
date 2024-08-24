#include <stdio.h>
#include <iostream>
using namespace std;

double f(double x);
//待积分函数
double f_y(double x);
//待积分函数
double f_x(double x);
//待积分函数

double err(double real,double test){
    return real - test;
}


double x(double a, double b, int n, int k){
    //返回[a,b]n等分的第k个点  k从0 到 n-1;
    double h = (b-a)/(double)(n-1);
    double result = a + k*h;
    return result;
}


double T(double a, double b, int n){
    //梯形公式
    //a为左端点，b为右端点，n为采样点个数2，
    
    double h = (b-a)/(double)(n-1);
    
    double result = 0.5 * h * (f(a) + f(b));

    return result;

}


double S(double a, double b, int n){
    //辛普森公式
    //a为左端点，b为右端点，n为采样点个数3，
    double h = (b-a)/(double)(n-1);

    double result = (1.0/3) * h * ( f(a) + 4*f(x(a,b,n,1)) + f(b));

    return result;

}


double S_3_8(double a, double b, int n){
    //辛普森3/8公式
    //a为左端点，b为右端点，n为采样点个数4，

    double h = (b-a)/(double)(n-1);   

    double result = (3.0/8) * h * ( f(a) + 3*f(x(a,b,n,1)) + 3*f(x(a,b,n,2)) + f(b));
    
    return result; 
}

double S_3_8_x(double a, double b, int n){
    //辛普森3/8公式
    //a为左端点，b为右端点，n为采样点个数4，

    double h = (b-a)/(double)(n-1);   

    double result = (3.0/8) * h * ( f_x(a) + 3*f_x(a+1) + 3*f_x(a+2)) + f_x(a+3);
    
    return result; 
}
double S_3_8_y(double a, double b, int n){
    //辛普森3/8公式
    //a为左端点，b为右端点，n为采样点个数4，

    double h = (b-a)/(double)(n-1);   

    double result = (3.0/8) * h * ( f_y(a) + 3*f_y(a+1) + 3*f_y(a+2) + f_y(a+3));
    
    return result; 
}


double B(double a, double b, int n){
    //布尔公式
    //a为左端点，b为右端点，n为采样点个数5，

    double h = (b-a)/(double)(n-1);

    double result = (2.0/45) * h * ( 7*f(a) + 32*f(x(a,b,n,1)) + 12*f(x(a,b,n,2)) + 32*f(x(a,b,n,3)) + 7*f(b) );

    return result;
}


double C_T(double a, double b, int n){
    //组合梯形公式
    //a为左端点，b为右端点，n为采样点个数
    int M = n-1;
    double h = (b-a)/(double)(n-1);

    double result = (1.0/2) * h * (f(a) + f(b));
    for (int k = 1; k < M; k++) {
        result += h * f(x(a,b,n,k));
    }

    return result;
}

double C_S(double a, double b, int n){
    //组合辛普森公式
    //a为左端点，b为右端点，n为采样点个数(n为奇数)
    //int M = (n-1)/2;
    double h = (b-a)/(double)(n-1);

    double result = (1.0/3) * h * (f(a) + f(b));
    for(int i = 2; i < (n-1); i+=2) {
        result += (2.0/3) * h * f(x(a,b,n,i));
    }
    for(int j = 1; j < (n-1); j+=2) {
        result += (4.0/3) * h * f(x(a,b,n,j));
    }

    return result;
}

double G_transform_f(double a, double b, double x){
    //高斯-勒让德变换
    double k1 = (b-a)/2.0;
    double k2 = (b+a)/2.0;
    double result = f(k2 + k1 * x)*k1;

    return result;
}

double G_T_y(double a, double b,double x){
    double k1 = (b-a)/2.0;
    double k2 = (b+a)/2.0;
    double result = f_y(k2 + k1 * x)*k1;
    return result;
}

double G_T_x(double a, double b,double x){
    double k1 = (b-a)/2.0;
    double k2 = (b+a)/2.0;
    double result = f_x(k2 + k1 * x)*k1;
    return result;
}

double G_2(double a, double b){
    //两点高斯-勒让德
    
    double result = G_transform_f(a,b,-0.5773502692) + G_transform_f(a,b,0.5773502692);

    return result;
}

double G_3(double a, double b){
    //3点高斯-勒让德
    
    double result = (5.0/9) * G_transform_f(a,b,-0.77459666924) + (8.0/9) * G_transform_f(a,b,0) + (5.0/9) * G_transform_f(a,b,0.77459666924);

    return result;
}

double G_3_x(double a, double b){
    //3点高斯-勒让德
    double result = (5.0/9) * G_T_x(a,b,-0.77459666924) + (8.0/9) * G_T_x(a,b,0) + (5.0/9) * G_T_x(a,b,0.77459666924);

    return result;
}

double G_3_y(double a, double b){
    //3点高斯-勒让德
    double result = (5.0/9) * G_T_y(a,b,-0.77459666924) + (8.0/9) * G_T_y(a,b,0) + (5.0/9) * G_T_y(a,b,0.77459666924);

    return result;
}


