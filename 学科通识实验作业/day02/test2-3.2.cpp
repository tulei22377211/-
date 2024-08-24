
#include <iostream>
#include <sstream> //字符串流函数的调用准备
#include <iomanip> //格式化操作的调用准备
#define isetRow 3
#define jsetColumn 3

// #include "convert.hpp"

int main(int argc, char *argv[])
{
    std::cout.precision(8);

    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << std::endl;
    typedef double arry2d[isetRow][jsetColumn];
    typedef float  arry2f[isetRow][jsetColumn];
    typedef double arry1d[isetRow];
    typedef float  arry1f[isetRow];

    arry2f A = {2.0f,4.0f,-2.0f,1.0f,-3.0f,-3.0f,4.0f,2.0f,2.0f};
    arry1f b = {2.0f,-1.0f,3.0f};


    std::cout << "A = " << std::endl;
    for(int i = 0 ; i < isetRow;i++)
    {
        for(int j= 0;j <jsetColumn;j++)
        {
            
            std::cout <<std::setw(13)<< A[i][j]  << ",     ";
        }
        std::cout << std::endl;
    }
        std::cout << std::endl;
    std::cout << "b = " << std::endl;
    for(int i = 0 ; i < isetRow;i++)
    {

        std::cout <<std::setw(13)<<  b[i] << std::endl;
    }

    arry2f L={0.0};
    arry2f U={0.0};
    arry1f y={0.0};
    arry1f x={0.0};
    // for(int j = 0; j<isetRow;j++){
    //         for(int i = 0; i<isetRow;i++){

    //         }

    // }
    
    for(int j = 0; j<isetRow;j++){
        U[0][j] = A[0][j];
        L[j][j] = 1.0;
    }
    for(int i = 1; i<isetRow;i++){
        L[i][0] = A[i][0]/U[0][0];
    }
    std::cout << "U = " << std::endl;
    for(int i = 0 ; i < isetRow;i++)
    {
        for(int j= 0;j <jsetColumn;j++)
        {
            
            std::cout <<std::setw(13)<<  U[i][j]  << ",     ";
        }
        std::cout << std::endl;
    }

        std::cout << std::endl;
    std::cout << "L = " << std::endl;
    for(int i = 0 ; i < isetRow;i++)
    {
        for(int j= 0;j <jsetColumn;j++)
        {
            
            std::cout <<std::setw(13)<<  L[i][j]  << ",     ";
        }
        std::cout << std::endl;
    }

        std::cout << std::endl;



    double usumstep = 0.0;
    for(int k = 1; k<isetRow;k++){
        for(int j = k; j<isetRow;j++){
            usumstep = 0.0;
            for(int m =0;m<k;m++){
                usumstep = usumstep + L[k][m]*U[m][j];
            }
            U[k][j] = A[k][j] - usumstep;
        }
        for(int i = k+1; i<isetRow;i++){
            usumstep = 0.0;
            for(int m =0;m<k;m++){
                usumstep = usumstep + L[i][m]*U[m][k];
            }
            L[i][k] = (A[i][k] - usumstep) / U[k][k];
        }
    }

    y[0] = b[0];
    double ysumstep = 0.0;
    for(int k = 1; k<isetRow;k++){
        ysumstep = 0.0;
        for(int i =0; i<k; i++){
            ysumstep = ysumstep + L[k][i]* y[i];
        }
        y[k] = b[k] - ysumstep;
    }
    std::cout << "First step result of LU method is: " << std::endl;
    std::cout << "U = " << std::endl;
    for(int i = 0 ; i < isetRow;i++)
    {
        for(int j= 0;j <jsetColumn;j++)
        {
            
            std::cout <<std::setw(13)<<  U[i][j]  << ",     ";
        }
        std::cout << std::endl;
    }

        std::cout << std::endl;

        std::cout << std::endl;
    std::cout << "y = " << std::endl;
    for(int i = 0 ; i < isetRow;i++)
    {

        std::cout <<std::setw(13)<<  y[i] << std::endl;
    }

     std::cout << std::endl;

    x[isetRow-1] = y[isetRow-1]/U[isetRow-1][isetRow-1];

    double xsumstep = 0.0;
    for(int i = isetRow - 2; i>-1;i--){
        xsumstep = 0.0;
        for(int j =i+1; j<isetRow; j++){
            xsumstep = xsumstep + U[i][j]* x[j];
        }
        x[i] = (y[i] - xsumstep)/U[i][i];
    }
    std::cout << "The numerical result of sequential Gaussian method is: " << std::endl;
        std::cout << std::endl;
    std::cout << "x = " << std::endl;
    for(int i = 0 ; i < isetRow;i++)
    {

        std::cout <<std::setw(13)<<  x[i] << std::endl;
    }

     std::cout << std::endl;








    double stepAkk=0.0 ;
    for(int k = 0; k < isetRow-1; k++)
    {
    for(int i = k + 1; i < isetRow;i++)
    {
        stepAkk = A[i][k]/A[k][k];
        for(int j=k + 0;j <jsetColumn;j++)
        {
            
            A[i][j] = A[i][j] - A[k][j]* stepAkk;
        }
        b[i] = b[i]-b[k]*stepAkk;
    }
    }
     std::cout << std::endl;
     std::cout << std::endl;
    std::cout << "First step result of sequential Gaussian method is: " << std::endl;
    std::cout << "A = " << std::endl;
    for(int i = 0 ; i < isetRow;i++)
    {
        for(int j= 0;j <jsetColumn;j++)
        {
            
            std::cout <<std::setw(13)<<  A[i][j]  << ",     ";
        }
        std::cout << std::endl;
    }

        std::cout << std::endl;
    std::cout << "b = " << std::endl;
    for(int i = 0 ; i < isetRow;i++)
    {

        std::cout <<std::setw(13)<<  b[i] << std::endl;
    }

     std::cout << std::endl;

    stepAkk=0.0 ;
    for(int k = isetRow-1; k > 0; k--)
    {
        b[k] = b[k]/A[k][k];
        A[k][k] = 1.0;
    for(int i = k - 1; i > -1;i--)
    {
        stepAkk = A[i][k]/A[k][k];
        for(int j=k - 0;j >i;j--)
        {
            
            A[i][j] = A[i][j] - A[k][j]* stepAkk;
        }
        b[i] = b[i]-b[k]*stepAkk;
    }
    }
        b[0] = b[0]/A[0][0];
        A[0][0] = 1.0;

    std::cout << "The numerical result of sequential Gaussian method is: " << std::endl;
    std::cout << "A = " << std::endl;
    for(int i = 0 ; i < isetRow;i++)
    {
        for(int j= 0;j <jsetColumn;j++)
        {
            
            std::cout <<std::setw(13)<<  A[i][j]  << ",     ";
        }
        std::cout << std::endl;
    }

        std::cout << std::endl;
    std::cout << "b = " << std::endl;
    for(int i = 0 ; i < isetRow;i++)
    {

        std::cout <<std::setw(13)<<  b[i] << std::endl;
    }

     std::cout << std::endl;




    return 0;
}