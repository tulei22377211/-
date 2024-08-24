#include <stdio.h>
#include <math.h>
#include <stdlib.h>


double f(double x1,double x2) {
    return (x1-x2)/(x1*x1 + x2*x2 + 2);
}

double dfx1(double x1,double x2) {
    return (-1*x1*x1 + x2*x2 + 2*x1*x2 + 2) / ((x1*x1 + x2*x2 +2)*(x1*x1 + x2*x2 + 2));
}

double dfx2(double x1,double x2) {
    return (x2*x2 - x1*x1 - 2*x1*x2 - 2) / ((x1*x1 + x2*x2 +2)*(x1*x1 + x2*x2 + 2));
    
}

double dfx1x1(double x1,double x2) {
    return (2*(x1*x1*x1 - 3*x1*x1*x2 - 6*x1 - 3*x1*x2*x2 + x2*x2*x2 + 2*x2)) / ((x1*x1 + x2*x2 +2)*(x1*x1 + x2*x2 +2)*(x1*x1 + x2*x2 +2));
}

double dfx1x2(double x1,double x2) {
    return (2*(x1*x1*x1 + 3*x1*x1*x2 + 2*x1 - 3*x1*x2*x2 - x2*x2*x2 - 2*x2)) / ((x1*x1 + x2*x2 +2)*(x1*x1 + x2*x2 +2)*(x1*x1 + x2*x2 +2));
}

double dfx2x1(double x1,double x2) {
    return (-2*(x2*x2*x2 + 3*x2*x2*x1 + 2*x2 - 3*x2*x1*x1 - x1*x1*x1 - 2*x1)) / ((x1*x1 + x2*x2 +2)*(x1*x1 + x2*x2 +2)*(x1*x1 + x2*x2 +2));
}

double dfx2x2(double x1,double x2) {
    return (2*(-1*x2*x2*x2 + 3*x2*x2*x1 + 6*x2 + 3*x2*x1*x1 - x1*x1*x1 - 2*x1)) / ((x1*x1 + x2*x2 +2)*(x1*x1 + x2*x2 +2)*(x1*x1 + x2*x2 +2));
}

double newton_x1(double x1, double x2){
    double x1x1 = dfx1x1(x1,x2);
    double x1x2 = dfx1x2(x1,x2);
    double x2x1 = dfx2x1(x1,x2);
    double x2x2 = dfx2x2(x1,x2);

    double dx1 = dfx1(x1,x2);
    double dx2 = dfx2(x1,x2);

    double det = x1x1*x2x2 - x1x2*x2x1;

    double x1_new = x1 - (1/det)*(dx1*x2x2 - dx2*x2x1);

    return x1_new;
}

double newton_x2(double x1, double x2){
    double x1x1 = dfx1x1(x1,x2);
    double x1x2 = dfx1x2(x1,x2);
    double x2x1 = dfx2x1(x1,x2);
    double x2x2 = dfx2x2(x1,x2);

    double dx1 = dfx1(x1,x2);
    double dx2 = dfx2(x1,x2);

    double det = x1x1*x2x2 - x1x2*x2x1;

    double x2_new = x1 - (1/det)*(dx2*x1x1 - dx1*x1x2  );
    return x2_new;
}


void newton_solve(double *x1, double *x2) {

    double x1x1 = dfx1x1(*x1,*x2);
    double x1x2 = dfx1x2(*x1,*x2);
    double x2x1 = dfx2x1(*x1,*x2);
    double x2x2 = dfx2x2(*x1,*x2);

    double dx1 = dfx1(*x1,*x2);
    double dx2 = dfx2(*x1,*x2);

    double **A = NULL;
    int n = 2;
    

    A = (double **)malloc(n * sizeof(double *));
    for (int i = 0; i < n; i++) {
        A[i] = (double *)malloc(n * sizeof(double));
    }
    if(A == NULL) {
        printf("Memory allocation failed.\n");
        exit(1);
    }
    
    A[0][0] = x1x1;A[0][1] = x1x2;
    A[1][0] = x2x1;A[1][1] = x2x2;

    
    
    
    double* b = (double*)malloc(n * sizeof(double));
    b[0] = -dx1 + (*x1)*x1x1 + (*x2)*x1x2;
    b[1] = -dx2 + (*x1)*x2x1 + (*x2)*x2x2;
    //b[1]出现了误差，发现其误差由dx2引起
    //第一步b[1]应该是  0.56011233393660；
    //但计算结果为      0.5380709


    //列主元Gauss消去法
    for (int i = 0; i < n; i++) {
        int max_index = i;
        for (int j = i + 1; j < n; j++) {//找列最大值
            if (fabs(A[j][i]) > fabs(A[max_index][i])) {
                max_index = j;
            }
        }
        if (A[max_index][i] == 0) {//如果最大值为0，则矩阵为奇异矩阵
            printf("The matrix is singular.\n");
            exit(1);
        }
        if (max_index != i) {//交换两行
            for(int j = i; j < n; j++) {
                double temp = A[i][j];
                A[i][j] = A[max_index][j];
                A[max_index][j] = temp;
            }
            double temp_b = b[i];
            b[i] = b[max_index];
            b[max_index] = temp_b;
        }
        for (int j = i + 1; j < n; j++) {//消去法
            double factor = A[j][i] / A[i][i];
            A[j][i] = 0;
            for (int k = i + 1; k < n; k++) {
                A[j][k] -= factor * A[i][k];
            }
            b[j] -= factor * b[i];
        }
    }
                

    double* x = (double*)malloc(n * sizeof(double));
    // 解线性方程组
    for (int i = n - 1; i >= 0; i--) {
        double ad = b[i];
        for (int j = n-1; j >= i+1; j--) {
            ad -= A[i][j] * x[j];
        }
        x[i] = ad/A[i][i];
    }

    *x1 = x[0];
    *x2 = x[1];

    for (int i = 0; i < n; i++) {
        free(A[i]);
    }
    free(A);
    free(b);
    free(x);
    
}






int main() {
   
    double eps_x = 1e-5;
    double eps_dy = 1e-4;
    double eps_y = 1e-4;

    int count = 1;
    int MAX_COUNT = 10;

    // 牛顿迭代法
    double x1 = -0.3;
    double x2 = 0.2;
    
    double y = f(x1,x2);
    double dy_x1 = dfx1(x1,x2);
    double dy_x2 = dfx2(x1,x2);
    double det = dfx1x1(x1,x2)*dfx2x2(x1,x2) - dfx1x2(x1,x2)*dfx2x1(x1,x2);
    double x1x1 = dfx1x1(x1,x2);
    double x1x2 = dfx1x2(x1,x2);
    double x2x1 = dfx2x1(x1,x2);
    double x2x2 = dfx2x2(x1,x2);
    
    printf("\n");
    printf("Newton's iteration method:\n");
    printf("--------------------------\n");
    printf("iter\tx\t\ty\t\tf(x,y)\t\tdf(x)/dx\tdf(x)/dy\n");
    printf("%d\t%lf\t%lf\t%lf\t%lf\t%lf\n", count, x1,x2, f(x1,x2), dfx1(x1,x2), dfx2(x1,x2));
    


    while ( ((fabs(dfx1(x1,x2)) >= eps_y) || (fabs(dfx2(x1,x2)) >= eps_y)) &&   (count <= MAX_COUNT  )){
        /* //求逆解方程
        x1 =newton_x1(x1, x2);
        x2 = newton_x2(x1, x2); */

        //列主元高斯消去法
        newton_solve(&x1, &x2);

        y = f(x1,x2);
        dy_x1 = dfx1(x1,x2);
        dy_x2 = dfx2(x1,x2);

        double det = dfx1x1(x1,x2)*dfx2x2(x1,x2) - dfx1x2(x1,x2)*dfx2x1(x1,x2);
        double x1x1 = dfx1x1(x1,x2);
        double x1x2 = dfx1x2(x1,x2);
        double x2x1 = dfx2x1(x1,x2);
        double x2x2 = dfx2x2(x1,x2);

        printf("%d\t%lf\t%lf\t%lf\t%lf\t%lf\n", count+1, x1,x2, f(x1,x2), dfx1(x1,x2), dfx2(x1,x2));
        
        count++;
    }

    printf("--------------------------\n");
    if((count > MAX_COUNT) && ((fabs(dfx1(x1,x2)) >= eps_y) || (fabs(dfx2(x1,x2)) >= eps_y))  ) {
        printf("The result was not found within %d iterations.\n", MAX_COUNT);
    } else if(isnan(x1) || isnan(x2) || isnan(f(x1,x2))){
        printf("The result was not found within %d iterations.\n", MAX_COUNT);
    } else {
        printf("root is x=%lf, y=%lf, f(x,y)=%lf\n", x1,x2, f(x1,x2));
    }
    printf("\n");




    return 0;
}