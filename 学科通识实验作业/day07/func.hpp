#include "func.cpp"

double f(double x);
//待积分函数

double err(double real,double test);

double x(double a, double b, int n, int k);
//返回[a,b]n等分的第k个点  k从0 到 n-1;

double T();
//梯形公式
//a为左端点，b为右端点，n为采样点个数2，

double S();
//辛普森公式
//a为左端点，b为右端点，n为采样点个数3，

double S_3_8();
//辛普森3/8公式
//a为左端点，b为右端点，n为采样点个数4，

double B();
//布尔公式
//a为左端点，b为右端点，n为采样点个数5，

double C_T();
//组合梯形公式
//a为左端点，b为右端点，n为采样点个数

double C_S();
//组合辛普森公式
//a为左端点，b为右端点，n为采样点个数(n为奇数)

double G_transform_f();
//高斯-勒让德变换

double G_2();
//两点高斯-勒让德

double G_3();
//3点高斯-勒让德


/////////////////////////////////////////////





