double f(double x, double y); //定义微分方程的右端项

double x_t(double t, double x, double y); //定义微分方程的x项
double y_t(double t, double x, double y); //定义微分方程的y项

double F_Euler(double a, double b, double t0, double y0, double h){
    //用向前欧拉法求解y(t,y)以（t0,y0）为初值点的微分方程
    //h为步长,[a,b]为时间区间,t0为初值点的时刻,y0为初值点的函数值
    //默认计算的为t0到b的积分
    double result = y0;
    for (double t = t0 + h; t <= b; t += h) {
        result += h * f(t-h, y0); 
        y0 = result; //更新y0
    }
    return result;
}


double B_Euler(double a, double b, double t0, double y0, double h){
    //用向后欧拉法求解y(t,y)以（t0,y0）为初值点的微分方程
    //h为步长
    y0 = y0 + h * f(t0, y0);
    double result = y0;
    for (double t = t0 + h; t < b; t += h) {
        result += h * f(t, y0);
        y0 = result; //更新y0
    }
    return result;
}
    



double EC_Euler(double a, double b, double t0, double y0, double h){
    //用预估修正向后欧拉法求解y(t,y)以（t0,y0）为初值点的微分方程
    //h为步长
    double y1 = y0 + h * f(t0, y0);
    double result = y0; 
    for (double t = t0 + h; t <= b; t += h) {
        result += 0.5 * h * (f(t-h, y0) + f(t, y1)); //用y(t,y)的近似值来近似y(t+h,y+h*f(t,y))
        y0 = result; //更新y0
        y1 = y0 + h * f(t, y0); //更新y1
    }
    return result;
}


double RK_2(double a, double b, double t0, double y0, double h){
    //用2阶龙格库塔法求解y(t,y)以（t0,y0）为初值点的微分方程
    //h为步长
    double k1 = h * f(t0, y0);
    double k2 = h * f(t0 + 0.5 * h, y0 + 0.5 * k1);
    double result = y0 ;
    for (double t = t0 + h; t <= b; t += h) {
        result += k2; //更新y0
        k1 = h * f(t, result);
        k2 = h * f(t + 0.5 * h, result + 0.5 * k1);
    }
    return result;
}


double RK_3(double a, double b, double t0, double y0, double h){
    //用4阶龙格库塔法求解y(t,y)以（t0,y0）为初值点的微分方程
    //h为步长
    double k1 = h * f(t0, y0);
    double k2 = h * f(t0 + 0.5 * h, y0 + 0.5 * k1);
    double k3 = h * f(t0 + h, y0 - k1 + 2 * k2);
    double result = y0 ;
    for (double t = t0 + h; t <= b; t += h) {
        result += (k1 + 4 * k2 + k3) / 6.0; //更新y0
        k1 = h * f(t, result);
        k2 = h * f(t + 0.5 * h, result + 0.5 * k1);
        k3 = h * f(t + h, result - k1 + 2 * k2);
    }
    return result;
}

double RK_4(double a, double b, double t0, double y0, double h){
    //用4阶龙格库塔法求解y(t,y)以（t0,y0）为初值点的微分方程
    //h为步长
    double k1 = h * f(t0, y0);
    double k2 = h * f(t0 + 0.5 * h, y0 + 0.5 * k1);
    double k3 = h * f(t0 + 0.5 * h, y0 + 0.5 * k2);
    double k4 = h * f(t0 + h, y0 + k3);
    double result = y0 ;
    for (double t = t0 + h; t <= b; t += h) {
        result += (k1 + 2 * k2 + 2 * k3 + k4) / 6.0; //更新y0
        k1 = h * f(t, result);
        k2 = h * f(t + 0.5 * h, result + 0.5 * k1);
        k3 = h * f(t + 0.5 * h, result + 0.5 * k2);
        k4 = h * f(t + h, result + k3);
        
    }
    return result;
}

double RK_4_2_x(double a, double b, double t0, double x0,double y0, double h){
    //用4阶龙格库塔法求解x(t),y(t)以（t0,x0,y0）为初值点的微分方程
    //h为步长
    double f1 = h * x_t(t0, x0, y0);
    double g1 = h * y_t(t0, x0, y0);
    double f2 = h * x_t(t0 + 0.5 * h, x0 + 0.5 * f1, y0 + 0.5 * g1);
    double g2 = h * y_t(t0 + 0.5 * h, x0 + 0.5 * f1, y0 + 0.5 * g1);
    double f3 = h * x_t(t0 + 0.5 * h, x0 + 0.5 * f2, y0 + 0.5 * g2);
    double g3 = h * y_t(t0 + 0.5 * h, x0 + 0.5 * f2, y0 + 0.5 * g2);
    double f4 = h * x_t(t0 + h, x0 + f3, y0 + g3);
    double g4 = h * y_t(t0 + h, x0 + f3, y0 + g3);
    double result_x = x0;
    double result_y = y0;
    for (double t = t0 + h; t <= b; t += h) {
        result_x +=  (f1 + 2 * f2 + 2 * f3 + f4) / 6.0; //更新x0
        result_y +=  (g1 + 2 * g2 + 2 * g3 + g4) / 6.0; //更新y0
        f1 = h * x_t(t, result_x, result_y);
        g1 = h * y_t(t, result_x, result_y);
        f2 = h * x_t(t + 0.5 * h, result_x + 0.5 * f1, result_y + 0.5 * g1);
        g2 = h * y_t(t + 0.5 * h, result_x + 0.5 * f1, result_y + 0.5 * g1);
        f3 = h * x_t(t + 0.5 * h, result_x + 0.5 * f2, result_y + 0.5 * g2);
        g3 = h * y_t(t + 0.5 * h, result_x + 0.5 * f2, result_y + 0.5 * g2);
        f4 = h * x_t(t + h, result_x + f3, result_y + g3);
        g4 = h * y_t(t + h, result_x + f3, result_y + g3);
    }
    return result_x;
}

double RK_4_2_y(double a, double b, double t0, double x0,double y0, double h){
    //用4阶龙格库塔法求解x(t),y(t)以（t0,x0,y0）为初值点的微分方程
    //h为步长
    double f1 = h * x_t(t0, x0, y0);
    double g1 = h * y_t(t0, x0, y0);
    double f2 = h * x_t(t0 + 0.5 * h, x0 + 0.5 * f1, y0 + 0.5 * g1);
    double g2 = h * y_t(t0 + 0.5 * h, x0 + 0.5 * f1, y0 + 0.5 * g1);
    double f3 = h * x_t(t0 + 0.5 * h, x0 + 0.5 * f2, y0 + 0.5 * g2);
    double g3 = h * y_t(t0 + 0.5 * h, x0 + 0.5 * f2, y0 + 0.5 * g2);
    double f4 = h * x_t(t0 + h, x0 + f3, y0 + g3);
    double g4 = h * y_t(t0 + h, x0 + f3, y0 + g3);
    double result_x = x0;
    double result_y = y0;
    for (double t = t0 + h; t <= b; t += h) {
        result_x +=  (f1 + 2 * f2 + 2 * f3 + f4) / 6.0; //更新x0
        result_y +=  (g1 + 2 * g2 + 2 * g3 + g4) / 6.0; //更新y0
        f1 = h * x_t(t, result_x, result_y);
        g1 = h * y_t(t, result_x, result_y);
        f2 = h * x_t(t + 0.5 * h, result_x + 0.5 * f1, result_y + 0.5 * g1);
        g2 = h * y_t(t + 0.5 * h, result_x + 0.5 * f1, result_y + 0.5 * g1);
        f3 = h * x_t(t + 0.5 * h, result_x + 0.5 * f2, result_y + 0.5 * g2);
        g3 = h * y_t(t + 0.5 * h, result_x + 0.5 * f2, result_y + 0.5 * g2);
        f4 = h * x_t(t + h, result_x + f3, result_y + g3);
        g4 = h * y_t(t + h, result_x + f3, result_y + g3);
    }
    return result_y;
}
