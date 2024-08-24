#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cstdlib>


using namespace std;
#include "func.hpp"

double f(double x){
    return 0;
}

/* double xx_1(double t){
    ifstream inputFile("cp1.plt"); // 打开文件

    if (!inputFile) {
        cerr << "无法打开文件。\n";
        return 1;
    }

    string line;
    vector<double> X, Y, Cp;

    // 跳过头部信息（TITLE和VARIABLES行）
    getline(inputFile, line); // TITLE行
    getline(inputFile, line); // VARIABLES行

    // 读取数据
    while (getline(inputFile, line)) {
        istringstream iss(line);
        double x, y, cp;
        if (!(iss >> x >> y >> cp)) {
            cerr << "解析数据行出错：" << line << '\n';
            return 1;
        }
        X.push_back(x);
        Y.push_back(y);
        Cp.push_back(cp);
    }

    // 关闭文件
    inputFile.close();

    // 输出数据来验证读取是否正确
     std::cout << "读取到的数据点数量： " << X.size() << '\n';
    for (size_t i = 0; i < X.size(); ++i) {
        std::cout << X[i] << " " << Y[i] << " " << Cp[i] << '\n';
    } 

     // 使用gnuplot绘制图形
 
    // 首先，将数据写入一个临时文件
    std::ofstream tempFile("temp_data.txt");
    for (size_t i = 0; i < X.size(); ++i) {
        tempFile << X[i] << " " << Y[i] << '\n';
    }
    tempFile.close();

    // 调用gnuplot绘图
    std::string command = "gnuplot -p -e \"set terminal wxt size 800,600; plot 'temp_data.txt' with points\" ";
    system(command.c_str());

    // 记得删除临时文件
    std::remove("temp_data.txt");

    // 绘制Cp-X图
    std::ofstream tempFile2("temp_data2.txt");
    for (size_t i = 0; i < X.size(); ++i) {
        tempFile2 << X[i] << " " << Cp[i] << '\n';
    }
    tempFile2.close();

    // 调用gnuplot绘图
    std::string command2 = "gnuplot -p -e \"set terminal wxt size 800,600; plot 'temp_data2.txt' with points\" ";
    system(command2.c_str());

    // 记得删除临时文件
    std::remove("temp_data2.txt"); 

    // 计算拟合曲线
    int n =  X.size()+1;
    double *x_now = (double*)malloc(n * sizeof(double));
    double *y_now = (double*)malloc(n * sizeof(double));
    for (int i = 0; i < n; i++) {
        x_now[i] = i;
    }
    for (int i = 0; i < (n-1)/2; i++) {
        y_now[i] = X[i];
    }
    y_now[n/2] = 0;
    for (int i = (n+1)/2; i < n; i++) {
        y_now[i] = X[i-1];
    }

    double ** A = (double**)malloc(n * sizeof(double*));
    if (A == NULL) {
        printf("Memory allocation failed.\n");
        exit(1);
    }
    for(int i=0;i<n;i++) {
        A[i] = (double*)malloc((n+1) * sizeof(double));
        if (A[i] == NULL) {
            printf("Memory allocation failed.\n");
            exit(1);
        }
    }


    for(int i=0;i<n;i++) {
        for(int j=0;j<n+1;j++) {
            A[i][j] = 0;
        }
    }

    

    for (int i = 0; i < n; i++) {
        A[i][0] = x_now[i];
        A[i][1] = y_now[i];
    }

    for (int j = 2; j < n+1; j++) {
        for (int i = j-1; i < n; i++) {
            
            A[i][j] = (A[i][j-1] - A[i-1][j-1]) / (A[i][0] - A[i-(j-1)][0]);
            
        }
    }

    
    
    double result = A[0][1];
    for (int i = 1; i < n; i++) {
        double sum = 1;
        for (int k = 0; k < i;k++){
            sum = sum * (t - A[k][0]);
        }
        result = result + A[i][i+1]*sum; 
    }


    //释放内存
    for(int i=0;i<n;i++) {  
        free(A[i]);
    }
    free(A);
    free(x_now);
    free(y_now);

    return result;
    
} */

void calcDerivative(const double* x, const double* y, int n, double* dy_dx, double* dx_dy) {
    //通过finite difference method计算导数
    for (int i = 0; i < n; i++) {
        if (i == 0) {
            //边界点处理
            dy_dx[i] = (y[i + 1] - y[i]) / (x[i + 1] - x[i]);
            dx_dy[i] = 1.0 / dy_dx[i];
        } else if (i == n - 1) {
            //边界点处理
            dy_dx[i] = (y[i] - y[i - 1]) / (x[i] - x[i - 1]);
            dx_dy[i] = 1.0 / dy_dx[i];
        } else {
            //内部点处理
            dy_dx[i] = (y[i + 1] - y[i - 1]) / (x[i + 1] - x[i - 1]);
            dx_dy[i] = (x[i + 1] - x[i - 1]) / (y[i + 1] - y[i - 1]);
        }
    }
}

void normal(const double* dy_dx, const double* dx_dy, int n, double* nx, double* ny) {
    //计算单位法向量
    for (int i = 0; i < n/2; i++) {
        if (dy_dx[i] > 0) {
            nx[i] = -1;
            ny[i] = dx_dy[i];
        } else {
            nx[i] = 1;
            ny[i] = -dx_dy[i];
        } 
        double len = sqrt(nx[i] * nx[i] + ny[i] * ny[i]);
        nx[i] /= len;
        ny[i] /= len;
    }
    for (int i = n/2; i < n; i++) {
        if (dy_dx[i] > 0) {
            nx[i] = 1;
            ny[i] = -dx_dy[i];
        }else {
            nx[i] = -1;
            ny[i] = dx_dy[i];
        }
        double len = sqrt(nx[i] * nx[i] + ny[i] * ny[i]);
        nx[i] /= len;
        ny[i] /= len;
    }
}

double f_y(double x){
    ifstream inputFile("cp1.plt"); // 打开文件

    if (!inputFile) {
        cerr << "无法打开文件。\n";
        return 1;
    }

    string line;
    vector<double> X, Y, Cp;

    // 跳过头部信息（TITLE和VARIABLES行）
    getline(inputFile, line); // TITLE行
    getline(inputFile, line); // VARIABLES行

    // 读取数据
    while (getline(inputFile, line)) {
        istringstream iss(line);
        double x, y, cp;
        if (!(iss >> x >> y >> cp)) {
            cerr << "解析数据行出错：" << line << '\n';
            return 1;
        }
        X.push_back(x);
        Y.push_back(y);
        Cp.push_back(cp);
    }

    // 关闭文件
    inputFile.close();

    int n = X.size();
    /* // 输出数据来验证读取是否正确
    std::cout << "读取到的数据点数量： " << X.size() << '\n';
    for (size_t i = 0; i < X.size(); ++i) {
        std::cout << X[i] << " " << Y[i] << " " << Cp[i] << '\n';
    } */

  
    
    double* dy_dx = (double*)malloc(n * sizeof(double));
    double* dx_dy = (double*)malloc(n * sizeof(double));
    
    calcDerivative(X.data(), Y.data(), n, dy_dx, dx_dy);

    /* for (int i = 0; i < n; i++) {
        std::cout << X[i] << " " << Y[i] << " " <<  "dy_dx[" << i << "] = " << dy_dx[i] << ' ';
        std::cout << "dx_dy[" << i << "] = " << dx_dy[i] << '\n';
    } */

    /*   // 使用gnuplot绘制图形

  
 
    // 首先，将数据写入一个临时文件
    std::ofstream tempFile("temp_data.txt");
    for (size_t i = 0; i < X.size(); ++i) {
        tempFile << X[i] << " " << Y[i] << '\n';
    }
    tempFile.close();

    // 调用gnuplot绘图
    std::string command = "gnuplot -p -e \"set terminal wxt size 800,600; plot 'temp_data.txt' with points\" ";
    system(command.c_str());

    // 记得删除临时文件
    std::remove("temp_data.txt");



    // 首先，将数据写入一个临时文件
    std::ofstream tempFile2("temp_data2.txt");
    for (size_t i = 0; i < X.size(); ++i) {
        tempFile2 << X[i] << " " << Cp[i] << '\n';
    }
    tempFile2.close();

    // 调用gnuplot绘图
    std::string command2 = "gnuplot -p -e \"set terminal wxt size 800,600; plot 'temp_data2.txt' with points\" ";
    system(command2.c_str());

    // 记得删除临时文件
    std::remove("temp_data2.txt"); */

    double * nx = (double*)malloc(n * sizeof(double));
    double * ny = (double*)malloc(n * sizeof(double));

    normal(dy_dx, dx_dy, n, nx, ny);

   /*  for (int i = 0; i < n; i++) {   
        std::cout << "nx[" << i << "] = " << nx[i] << " ny[" << i << "] = " << ny[i] << '\n';
    } */
    int i = (int)x;
    double result = Cp[i] * ny[i] * sqrt(1 + dy_dx[i] * dy_dx[i]);

    //释放内存
    free(dy_dx);
    free(dx_dy);
    free(nx);
    free(ny);
    
    return result;
}

double f_x(double x){
    ifstream inputFile("cp1.plt"); // 打开文件

    if (!inputFile) {
        cerr << "无法打开文件。\n";
        return 1;
    }

    string line;
    vector<double> X, Y, Cp;

    // 跳过头部信息（TITLE和VARIABLES行）
    getline(inputFile, line); // TITLE行
    getline(inputFile, line); // VARIABLES行

    // 读取数据
    while (getline(inputFile, line)) {
        istringstream iss(line);
        double x, y, cp;
        if (!(iss >> x >> y >> cp)) {
            cerr << "解析数据行出错：" << line << '\n';
            return 1;
        }
        X.push_back(x);
        Y.push_back(y);
        Cp.push_back(cp);
    }

    // 关闭文件
    inputFile.close();

    int n = X.size();
    /* // 输出数据来验证读取是否正确
    std::cout << "读取到的数据点数量： " << X.size() << '\n';
    for (size_t i = 0; i < X.size(); ++i) {
        std::cout << X[i] << " " << Y[i] << " " << Cp[i] << '\n';
    } */

  
    
    double* dy_dx = (double*)malloc(n * sizeof(double));
    double* dx_dy = (double*)malloc(n * sizeof(double));
    
    calcDerivative(X.data(), Y.data(), n, dy_dx, dx_dy);

    /* for (int i = 0; i < n; i++) {
        std::cout << X[i] << " " << Y[i] << " " <<  "dy_dx[" << i << "] = " << dy_dx[i] << ' ';
        std::cout << "dx_dy[" << i << "] = " << dx_dy[i] << '\n';
    } */

    /*   // 使用gnuplot绘制图形

  
 
    // 首先，将数据写入一个临时文件
    std::ofstream tempFile("temp_data.txt");
    for (size_t i = 0; i < X.size(); ++i) {
        tempFile << X[i] << " " << Y[i] << '\n';
    }
    tempFile.close();

    // 调用gnuplot绘图
    std::string command = "gnuplot -p -e \"set terminal wxt size 800,600; plot 'temp_data.txt' with points\" ";
    system(command.c_str());

    // 记得删除临时文件
    std::remove("temp_data.txt");



    // 首先，将数据写入一个临时文件
    std::ofstream tempFile2("temp_data2.txt");
    for (size_t i = 0; i < X.size(); ++i) {
        tempFile2 << X[i] << " " << Cp[i] << '\n';
    }
    tempFile2.close();

    // 调用gnuplot绘图
    std::string command2 = "gnuplot -p -e \"set terminal wxt size 800,600; plot 'temp_data2.txt' with points\" ";
    system(command2.c_str());

    // 记得删除临时文件
    std::remove("temp_data2.txt"); */

    double * nx = (double*)malloc(n * sizeof(double));
    double * ny = (double*)malloc(n * sizeof(double));

    normal(dy_dx, dx_dy, n, nx, ny);

    /* for (int i = 0; i < n; i++) {   
        std::cout << "nx[" << i << "] = " << nx[i] << " ny[" << i << "] = " << ny[i] << '\n';
    } */
    int i = (int)x;
    double result = Cp[i] * nx[i] * sqrt(1 + dy_dx[i] * dy_dx[i]);

    //释放内存
    free(dy_dx);
    free(dx_dy);
    free(nx);
    free(ny);
    
    return result;
}

int main() {

    ifstream inputFile("cp1.plt"); // 打开文件

    if (!inputFile) {
        cerr << "无法打开文件。\n";
        return 1;
    }

    string line;
    vector<double> X, Y, Cp;

    // 跳过头部信息（TITLE和VARIABLES行）
    getline(inputFile, line); // TITLE行
    getline(inputFile, line); // VARIABLES行

    // 读取数据
    while (getline(inputFile, line)) {
        istringstream iss(line);
        double x, y, cp;
        if (!(iss >> x >> y >> cp)) {
            cerr << "解析数据行出错：" << line << '\n';
            return 1;
        }
        X.push_back(x);
        Y.push_back(y);
        Cp.push_back(cp);
    }

    // 关闭文件
    inputFile.close();

    int n = X.size();
    // 输出数据来验证读取是否正确
    std::cout << "读取到的数据点数量： " << X.size() << '\n';
    /* for (size_t i = 0; i < X.size(); ++i) {
        std::cout << X[i] << " " << Y[i] << " " << Cp[i] << '\n';
    } */

      // 使用gnuplot绘制图形

  
 
    // 首先，将数据写入一个临时文件
    std::ofstream tempFile("temp_data.txt");
    for (size_t i = 0; i < X.size(); ++i) {
        tempFile << X[i] << " " << Y[i] << '\n';
    }
    tempFile.close();

    // 调用gnuplot绘图
    std::string command = "gnuplot -p -e \"set terminal wxt size 800,600; plot 'temp_data.txt' with points\" ";
    system(command.c_str());

    // 记得删除临时文件
    std::remove("temp_data.txt");



    // 首先，将数据写入一个临时文件
    std::ofstream tempFile2("temp_data2.txt");
    for (size_t i = 0; i < X.size(); ++i) {
        tempFile2 << X[i] << " " << Cp[i] << '\n';
    }
    tempFile2.close();

    // 调用gnuplot绘图
    std::string command2 = "gnuplot -p -e \"set terminal wxt size 800,600; plot 'temp_data2.txt' with points\" ";
    system(command2.c_str());

    // 记得删除临时文件
    std::remove("temp_data2.txt");

    double a = 0.0;
    double b = 63.0;
    
    double F_x = 0.0;
    double F_y = 0.0;

    for (double x = a; x < b; x += 3.0) {
        F_x += S_3_8_x(x,x+3.0,4);
        F_y += S_3_8_y(x,x+3.0,4);
    }

    std::cout << "F_x = " << F_x << '\n';
    std::cout << "F_y = " << F_y << '\n';

    return 0;
}
