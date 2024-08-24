
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

    double stepAkk=0.0 ;

    arry1f astep;
    double bstep;
    for(int k = 0; k < isetRow-1; k++)
    {
        
        int imax = k;
        for(int i = k+1; i<isetRow;i++){
            if(A[i][k] > A[k][k]) imax = i;
        }
        if(imax != k){
        for(int j =0;j<isetRow;j++){
            astep[j] = A[k][j];
        }
        bstep = b[k];
        for(int j =0;j<isetRow;j++){
            A[k][j]= A[imax][j];
        }
        for(int j =0;j<isetRow;j++){
            A[imax][j]= astep[j];
        }
        b[k] = b[imax];
        b[imax] = bstep;
        }
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
    std::cout << "First step result of Gauss-Jordan method is: " << std::endl;
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

    std::cout << "The numerical result of Gauss-Jordan method is: " << std::endl;
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