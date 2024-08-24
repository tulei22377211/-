#include <stdio.h>
#include <math.h>
#include <stdlib.h>


double f(double x1,double x2) {
    return 0.5*x1*x1 + 2.5*x2*x2;
}

double dfx1(double x1,double x2) {
    return x1;
}

double dfx2(double x1,double x2) {
    return 5*x2;
    
}

double dfx1x1(double x1,double x2) {
    return 1;
}

double dfx1x2(double x1,double x2) {
    return 0;
}

double dfx2x1(double x1,double x2) {
    return 0;
}

double dfx2x2(double x1,double x2) {
    return 5;
}



    
    
    
    
    







int main() {
   
    double eps_x = 1e-5;
    double eps_dy = 1e-4;
    double eps_y = 1e-4;

    int count = 1;
    int MAX_COUNT = 100;

    //最速下降法
    double x1 = 5;
    double x2 = 1;
    
    double y = f(x1,x2);
    double dy_x1 = dfx1(x1,x2);
    double dy_x2 = dfx2(x1,x2);
    
    double x1x1 = dfx1x1(x1,x2);
    double x1x2 = dfx1x2(x1,x2);
    double x2x1 = dfx2x1(x1,x2);
    double x2x2 = dfx2x2(x1,x2);

    double dx1 = dfx1(x1,x2);
    double dx2 = dfx2(x1,x2);
    
    printf("\n");
    printf("\n");
    printf("--------------------------\n");
    printf("iter\tx\t\ty\t\tf(x,y)\t\tdf(x)/dx\tdf(x)/dy\n");
    printf("%d\t%lf\t%lf\t%lf\t%lf\t%lf\n", count, x1,x2, f(x1,x2), dfx1(x1,x2), dfx2(x1,x2));
    

    double **Q = NULL;
    int n = 2;
    

    Q = (double **)malloc(n * sizeof(double *));
    for (int i = 0; i < n; i++) {
        Q[i] = (double *)malloc(n * sizeof(double));
    }
    if(Q == NULL) {
        printf("Memory allocation failed.\n");
        exit(1);
    }
    
    Q[0][0] = x1x1;Q[0][1] = x1x2;
    Q[1][0] = x2x1;Q[1][1] = x2x2;


    double *d = (double *)malloc(n * sizeof(double));
    if(d == NULL) {
        printf("Memory allocation failed.\n");
        exit(1);    
    }

    


    double * g = (double *)malloc(n * sizeof(double));
    if(g == NULL) {
        printf("Memory allocation failed.\n");
        exit(1);    
    }   
    g[0] = dx1;
    g[1] = dx2;

    d[0] = -g[0];
    d[1] = -g[1];

    int result = 1;


    while (((fabs(g[0])>eps_y) || (fabs(g[1])>eps_y)) &&   (count <= MAX_COUNT  )){
        
        result = 0;
        double x1x1 = dfx1x1(x1,x2);
        double x1x2 = dfx1x2(x1,x2);
        double x2x1 = dfx2x1(x1,x2);
        double x2x2 = dfx2x2(x1,x2);

        double dx1 = dfx1(x1,x2);
        double dx2 = dfx2(x1,x2);

        

        double ak = -1 * (g[0]*d[0] + g[1]*d[1])/((Q[0][0]*d[0] + Q[0][1]*d[1])*d[0] + (Q[1][0]*d[0] + Q[1][1]*d[1]) * d[1]);
        
        x1 = x1 + ak * d[0];
        x2 = x2 + ak * d[1];

        g[0] = dfx1(x1, x2);
        g[1] = dfx2(x1, x2);

        if ((fabs(g[0])<eps_y)&&(fabs(g[1])<eps_y)) {
            result = 1;
            break;
        } else {
            //double bk =  ((Q[0][0]*d[0] + Q[0][1]*d[1])*g[0] + (Q[1][0]*d[0] + Q[1][1]*d[1]) * g[1])/((Q[0][0]*d[0] + Q[0][1]*d[1])*d[0] + (Q[1][0]*d[0] + Q[1][1]*d[1]) * d[1]);
            double bk =  (g[0]*g[0] + g[1]*g[1])/(d[0]*d[0] + d[1]*d[1]);
            d[0] = -g[0] + bk*d[0];
            d[1] = -g[1] + bk*d[1];
        }


        y = f(x1,x2);
        dy_x1 = dfx1(x1,x2);
        dy_x2 = dfx2(x1,x2);

        x1x1 = dfx1x1(x1,x2);
        x1x2 = dfx1x2(x1,x2);
        x2x1 = dfx2x1(x1,x2);
        x2x2 = dfx2x2(x1,x2);

        dx1 = dfx1(x1,x2);
        dx2 = dfx2(x1,x2);


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



    //释放内存
    for (int i = 0; i < n; i++) {
        free(Q[i]);
    }
    free(Q);
    free(d);
    free(g);


    return 0;
}