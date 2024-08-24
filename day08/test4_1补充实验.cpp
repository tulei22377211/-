#include <iostream>

#define PI 3.14159265358979323846264338327950288419716939937510

using namespace std;

double x_t(double t, double x, double y,double p,double q) {
    double K = 1000.0;
    double l = 6.0;
    double a = 0.2;
    double omega = 2*PI*(38/60.0);
    double W = 80.0;
    double m = 2500.0;
    double d = 0.01;

    double result = q;
    return result;
}

double dx_t(double t, double x, double y, double p, double q) {
    double K = 1000.0;
    double l = 6.0;
    double a = 0.2;
    double omega = 2*PI*(38/60.0);
    double W = 80.0;
    double m = 2500.0;
    double d = 0.01;

  
    //double result = -q*d + 3*cos(x)*K*(1.0/(l*m*a)) * (exp(a*(y-l*sin(x))) - exp(a*(y)+sin(x)));
    //题目中所给公式 印刷错误，应该是：
    double result = -q*d + 3*cos(x)*K*(1.0/(l*m*a)) * (exp(a*(y-l*sin(x))) - exp(a*(y+l*sin(x))));

    return result;
}

double y_t(double t, double x, double y, double p, double q) {
    double K = 1000.0;
    double l = 6.0;
    double a = 0.2;
    double omega = 2*PI*(38/60.0);
    double W = 80.0;
    double m = 2500.0;
    double d = 0.01;

    double result = p;
    return result;
}

double dy_t(double t, double x, double y,double p, double q) {
    double K = 1000.0;
    double l = 6.0;
    double a = 0.2;
    double omega = 2*PI*(38/60.0);
    double W = 80.0;
    double m = 2500.0;
    double d = 0.01;

    double result = -d*p - K*(1.0/(m*a)) * (exp(a*(y-l*sin(x))) + exp(a*y+l*sin(x)) -2) + 0.2*W*sin(omega*t);
    return result;
}

double RK_4_4_x(double a, double b, double t0, double x0, double dx0, double y0,double dy0, double h) {
    //用4阶龙格库塔法求解x(t),dx(t),y(t),dy(t)以（t0,x0,dx0,y0,dy0）为初值点的微分方程
    double f1 = h * x_t(t0, x0, y0, dy0, dx0);
    double g1 = h * dx_t(t0, x0, y0, dy0, dx0);
    double h1 = h * y_t(t0, x0, y0, dy0, dx0);
    double k1 = h * dy_t(t0, x0, y0,dy0, dx0);
    double f2 = h * x_t(t0 + h / 2.0, x0 + f1 / 2.0, y0 + h1 / 2.0, dy0 + k1 / 2.0, dx0 + g1 / 2.0);
    double g2 = h * dx_t(t0 + h / 2.0, x0 + f1 / 2.0, y0 + h1 / 2.0, dy0 + k1 / 2.0, dx0 + g1 / 2.0);
    double h2 = h * y_t(t0 + h / 2.0, x0 + f1 / 2.0, y0 + h1 / 2.0, dy0 + k1 / 2.0, dx0 + g1 / 2.0);
    double k2 = h * dy_t(t0 + h / 2.0, x0 + f1 / 2.0, y0 + h1 / 2.0, dy0 + k1 / 2.0, dx0 + g1 / 2.0);
    double f3 = h * x_t(t0 + h / 2.0, x0 + f2 / 2.0, y0 + h2 / 2.0, dy0 + k2 / 2.0, dx0 + g2 / 2.0);
    double g3 = h * dx_t(t0 + h / 2.0, x0 + f2 / 2.0, y0 + h2 / 2.0, dy0 + k2 / 2.0, dx0 + g2 / 2.0);
    double h3 = h * y_t(t0 + h / 2.0, x0 + f2 / 2.0, y0 + h2 / 2.0, dy0 + k2 / 2.0, dx0 + g2 / 2.0);
    double k3 = h * dy_t(t0 + h / 2.0, x0 + f2 / 2.0, y0 + h2 / 2.0, dy0 + k2 / 2.0, dx0 + g2 / 2.0);
    double f4 = h * x_t(t0 + h, x0 + f3, y0 + h3, dy0 + k3, dx0 + g3);
    double g4 = h * dx_t(t0 + h, x0 + f3, y0 + h3, dy0 + k3, dx0 + g3);
    double h4 = h * y_t(t0 + h, x0 + f3, y0 + h3, dy0 + k3, dx0 + g3);
    double k4 = h * dy_t(t0 + h, x0 + f3, y0 + h3, dy0 + k3, dx0 + g3);
    double result_x = x0 ;
    double result_dx = dx0;
    double result_y = y0 ;
    double result_dy = dy0 ;
    for (double t = t0 + h; t <= b; t += h) {
        result_x += (f1 + 2.0 * f2 + 2.0 * f3 + f4) / 6.0;
        result_dx += (g1 + 2.0 * g2 + 2.0 * g3 + g4) / 6.0;
        result_y += (h1 + 2.0 * h2 + 2.0 * h3 + h4) / 6.0;
        result_dy += (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0;
        f1 = h * x_t(t, result_x, result_y, result_dy, result_dx);
        g1 = h * dx_t(t, result_x, result_y, result_dy, result_dx);
        h1 = h * y_t(t, result_x, result_y, result_dy, result_dx);
        k1 = h * dy_t(t, result_x, result_y, result_dy, result_dx);
        f2 = h * x_t(t + h / 2.0, result_x + f1 / 2.0, result_y + h1 / 2.0, result_dy + k1 / 2.0, result_dx + g1 / 2.0);      
        g2 = h * dx_t(t + h / 2.0, result_x + f1 / 2.0, result_y + h1 / 2.0, result_dy + k1 / 2.0, result_dx + g1 / 2.0);
        h2 = h * y_t(t + h / 2.0, result_x + f1 / 2.0, result_y + h1 / 2.0, result_dy + k1 / 2.0, result_dx + g1 / 2.0);
        k2 = h * dy_t(t + h / 2.0, result_x + f1 / 2.0, result_y + h1 / 2.0, result_dy + k1 / 2.0, result_dx + g1 / 2.0);
        f3 = h * x_t(t + h / 2.0, result_x + f2 / 2.0, result_y + h2 / 2.0, result_dy + k2 / 2.0, result_dx + g2 / 2.0);
        g3 = h * dx_t(t + h / 2.0, result_x + f2 / 2.0, result_y + h2 / 2.0, result_dy + k2 / 2.0, result_dx + g2 / 2.0);
        h3 = h * y_t(t + h / 2.0, result_x + f2 / 2.0, result_y + h2 / 2.0, result_dy + k2 / 2.0, result_dx + g2 / 2.0);
        k3 = h * dy_t(t + h / 2.0, result_x + f2 / 2.0, result_y + h2 / 2.0, result_dy + k2 / 2.0, result_dx + g2 / 2.0);
        f4 = h * x_t(t + h, result_x + f3, result_y + h3, result_dy + k3, result_dx + g3);
        g4 = h * dx_t(t + h, result_x + f3, result_y + h3, result_dy + k3, result_dx + g3);
        h4 = h * y_t(t + h, result_x + f3, result_y + h3, result_dy + k3, result_dx + g3);
        k4 = h * dy_t(t + h, result_x + f3, result_y + h3, result_dy + k3, result_dx + g3);
    }
    return result_x;
}




int main() {

    //x = 0;
    //y = y;

    double a = 0.0;
    double b = 600.0;
    double t0 = 0.0;
    double x0 = 0.00000001;//初始时只有一点点倾斜角
                            //若初始时没有倾斜角，x一直为0
    double dx0 = 0.0;
    double y0 = 0.0;
    double dy0 = 0.0;
    double h = 1.0;

    cout << fixed;
    cout.precision(8);

    cout << "t\t\tx" << endl;
    cout << "----------" << endl;
    for (int i = 0; i*h <= b; i++){
        double t = t0 + i*h;
        double x = RK_4_4_x(a, t, t0, x0,dx0, y0,dy0, h);
        cout << t << "\t" << x << endl;
    }

    return 0;
}