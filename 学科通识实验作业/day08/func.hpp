#include <iostream>
#include "func.cpp"

using namespace std;

double f(double x, double y);
//y(t,y)在（t,y）处的导数

//方程组x(t),y(t)的微分方程
double x_t(double t, double x, double y); //定义微分方程的x项
double y_t(double t, double x, double y); //定义微分方程的y项

double F_Euler(double a, double b, double t0, double y0, double h);
//用向前欧拉法求解y(t,y)以（t0,y0）为初值点的微分方程

double B_Euler(double a, double b, double t0, double y0, double h);
//用向后欧拉法求解y(t,y)以（t0,y0）为初值点的微分方程

double EC_Euler(double a, double b, double t0, double y0, double h);
//用预估修正向后欧拉法求解y(t,y)以（t0,y0）为初值点的微分方程

double RK_2(double a, double b, double t0, double y0, double h);
//用2阶Runge-Kutta法求解y(t,y)以（t0,y0）为初值点的微分方程

double RK_3(double a, double b, double t0, double y0, double h);
//用3阶Runge-Kutta法求解y(t,y)以（t0,y0）为初值点的微分方程

double RK_4(double a, double b, double t0, double y0, double h);
//用4阶Runge-Kutta法求解y(t,y)以（t0,y0）为初值点的微分方程

double RK_4_2_x(double a, double b, double t0, double x0,double y0, double h);
//用4阶Runge-Kutta法求解x(t),y(t)以（t0,x0,y0）为初值点的微分方程

double RK_4_2_y(double a, double b, double t0, double x0,double y0, double h);
//用4阶Runge-Kutta法求解x(t),y(t)以（t0,x0,y0）为初值点的微分方程

