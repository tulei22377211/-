#include <iostream>
#include <sstream> //字符串流函数的调用准备
#include <iomanip> //格式化操作的调用准备
#include "func6.hpp"
#include "cmath"

template<size_t n>
double norm_2(double X[n])
{
    double sum = 0;
    for (int k=0; k<n; k++)
        sum = sum + pow(X[k], 2);
    return sqrt(sum);
}

int main()
{
    std::cout.precision(12);
    double A[3][4] = {15600, 18760, 17610, 19170, 
                     7540, 2750, 14630, 610,
                     20140, 18610, 13480, 18390};
    double T[4] = {0.07074, 0.07220, 0.07690, 0.07242};
    double X[4] = {0, 0, 6370, 0};

    int cnt = 0; 
    while (cnt < 100)
    {    
        double P[4][4], Q[4];
        for (int i=0; i<4; i++)
        {
            double s = sqrt(pow(X[0]-A[0][i],2)+pow(X[1]- A[1][i],2)
                           +pow(X[2]-A[2][i],2));
            for (int j=0; j<3; j++)
                P[i][j] = (X[j] - A[j][i]) / s;
            P[i][3] = 299792.50;  //  注意index, 注意单位, 注意取精确值
            Q[i] = s - 299792.50 * (T[i] - X[3]);
        }
        double* value = Linear_Least_Squares<4, 4>(P, Q);
        if (norm_2<4>(value) < 1e-8)
        {
            std::cout << "---Convergence!---" << std::endl;
            return 0;
        }
        for (int k=0; k<4; k++)
            X[k] = X[k] - value[k];
        std::cout << "cnt="<< cnt << "\tx=" << X[0] << "\ty=" << X[1]
                  << "\tz=" << X[2] << "\td=" << X[3] << std::endl;
        cnt ++;
    }
    std::cout << "---Divergence!---" << std::endl;
    return 0;
}
