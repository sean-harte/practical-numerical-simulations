#include <iostream>
#include <iomanip>
#include<cmath>
#include <fstream>
#include<cstdlib>

using namespace std;
double f(double x,double y,double y1,double y2,double y3){ // this is the function for the fourth derivative of y wrt x
    double a = -10*y2 - (1.0/7.0)*x*pow(y,3); // y2 is the second derivative of y wrt x
    return a;
}

// runge kutta below is for the third derivative using the fourth derivative as its slope
double k1_y3(double x,double y,double y1,double y2,double y3) {
    double k1_y3 = f(x,y,y1,y2,y3);
    return k1_y3;
}

double k2_y3(double x,double y,double y1,double y2,double y3,double h){
    double k2_y3 = f(x + (h/2), y + (h/2)*y1, y1 + (h/2)*y2, y2 + (h/2)*y3, y3 + (h/2)*k1_y3(x,y,y1,y2,y3)); /// changed formula
    return k2_y3;
}

double k3_y3(double x,double y,double y1,double y2,double y3,double h){
    double k3_y3 = f(x + (h/2), y + h*(k2_y1(x,y,y1,y2,y3,h)/2), y1 + h*(k2_y2(x,y,y1,y2,y3,h)/2), y2 + h*(k2_y3(x,y,y1,y2,y3,h))/2);
    return k3_y3;
}

double k4_y3(double x,double y,double y1,double y2,double y3,double h){
    double k4_y3 = f(x + h, y + h*(k3_y1(x,y,y1,y2,y3,h)), y1 + h*(k3_y2(x,y,y1,y2,y3,h)), y2 + h*(k3_y3(x,y,y1,y2,y3,h)));
    return k4_y3;
}

double y3_next(double x,double y,double y1,double y2,double y3,double h){
    double y3_next = (h/6)*(k1_y3(x,y,y1,y2,y3) + 2*k2_y3(x,y,y1,y2,y3,h) + 2*k3_y3(x,y,y1,y2,y3,h) + k4_y3(x,y,y1,y2,y3,h));
    return y3_next;
}

// now doing rk method using y3 or the third derivative of y wrt x to get expression for y2 or y''
// it makes sense that y3_next is now our expression for the slope of y2
double k1_y2(double x,double y,double y1,double y2,double y3,double h) {
    double k1_y2 = y3_next(x,y,y1,y2,y3,h);
    return k1_y2;
}

double k2_y2(double x,double y,double y1,double y2,double y3,double h){
    double k2_y2 = y3_next(x + (h/2), y + h*(k1_y1(x,y,y1,y2,y3)/2), y1 + h*(k1_y2(x,y,y1,y2,y3)/2), y2 + h*(k1_y3(x,y,y1,y2,y3))/2);
    return k2_y2;
}

double k3_y2(double x,double y,double y1,double y2,double y3,double h){
    double k3_y2 = y3_next(x + (h/2), y + h*(k2_y1(x,y,y1,y2,y3,h)/2), y1 + h*(k2_y2(x,y,y1,y2,y3,h)/2), y2 + h*(k2_y3(x,y,y1,y2,y3,h))/2);
    return k3_y2;
}

double k4_y2(double x,double y,double y1,double y2,double y3,double h){
    double k4_y2 = y3_next(x + h, y + h*(k3_y1(x,y,y1,y2,y3,h)), y1 + h*(k3_y2(x,y,y1,y2,y3,h)), y2 + h*(k3_y3(x,y,y1,y2,y3,h)));
    return k4_y2;
}

double y2_next(double x,double y,double y1,double y2,double y3,double h){
    double y2_next = (h/6)*(k1_y2(x,y,y1,y2,y3) + 2*k2_y2(x,y,y1,y2,y3,h) + 2*k3_y2(x,y,y1,y2,y3,h) + k4_y2(x,y,y1,y2,y3,h));
    return y2_next;
}

// now we are using y2_next as the slope of y1 or y'
double k1_y1(double x,double y,double y1,double y2,double y3,double h) {
    double k1_y1 = y2_next(x,y,y1,y2,y3,h);
    return k1_y1;
}

double k2_y1(double x,double y,double y1,double y2,double y3,double h){
    double k2_y1 = y2_next(x + (h/2), y + h*(k1_y1(x,y,y1,y2,y3)/2), y1 + h*(k1_y2(x,y,y1,y2,y3)/2), y2 + h*(k1_y3(x,y,y1,y2,y3))/2);
    return k2_y1;
}

double k3_y1(double x,double y,double y1,double y2,double y3,double h){
    double k3_y1 = y2_next(x + (h/2), y + h*(k2_y1(x,y,y1,y2,y3,h)/2), y1 + h*(k2_y2(x,y,y1,y2,y3,h)/2), y2 + h*(k2_y3(x,y,y1,y2,y3,h))/2);
    return k3_y1;
}

double k4_y1(double x,double y,double y1,double y2,double y3,double h){
    double k4_y1 = y2_next(x + h, y + h*(k3_y1(x,y,y1,y2,y3,h)), y1 + h*(k3_y2(x,y,y1,y2,y3,h)), y2 + h*(k3_y3(x,y,y1,y2,y3,h)));
    return k4_y1;
}

double y1_next(double x,double y,double y1,double y2,double y3,double h){
    double y1_next = (h/6)*(k1_y1(x,y,y1,y2,y3) + 2*k2_y1(x,y,y1,y2,y3,h) + 2*k3_y1(x,y,y1,y2,y3,h) + k4_y1(x,y,y1,y2,y3,h));
    return y1_next;
}

// now using y1_next as the slope for y. finally
double k1_y0(double x,double y,double y1,double y2,double y3,double h) {
    double k1_y0 = y1_next(x,y,y1,y2,y3,h);
    return k1_y0;
}

double k2_y0(double x,double y,double y1,double y2,double y3,double h){
    double k2_y0 = y1_next(x + (h/2), y + h*(k1_y1(x,y,y1,y2,y3)/2), y1 + h*(k1_y2(x,y,y1,y2,y3)/2), y2 + h*(k1_y3(x,y,y1,y2,y3))/2);
    return k2_y0;
}

double k3_y0(double x,double y,double y1,double y2,double y3,double h){
    double k3_y0 = y1_next(x + (h/2), y + h*(k2_y1(x,y,y1,y2,y3,h)/2), y1 + h*(k2_y2(x,y,y1,y2,y3,h)/2), y2 + h*(k2_y3(x,y,y1,y2,y3,h))/2);
    return k3_y0;
}

double k4_y0(double x,double y,double y1,double y2,double y3,double h){
    double k4_y0 = y1_next(x + h, y + h*(k3_y1(x,y,y1,y2,y3,h)), y1 + h*(k3_y2(x,y,y1,y2,y3,h)), y2 + h*(k3_y3(x,y,y1,y2,y3,h)));
    return k4_y0;
}

double y0_next(double x,double y,double y1,double y2,double y3,double h){
    double y0_next = (h/6)*(k1_y0(x,y,y1,y2,y3) + 2*k2_y0(x,y,y1,y2,y3,h) + 2*k3_y0(x,y,y1,y2,y3,h) + k4_y0(x,y,y1,y2,y3,h));
    return y0_next;
}

// double y0i = 1.0;
// double y1i = 0.0;
// double y2i = 0.0;
// double y3i = 0.0;
// double xi = 0.0;
// double xf = 30.0;

int main(){
    const int len = 10;
    //double discrepency[len] = {};
    double H[len] = {1.0,0.5,0.25, 0.1,0.05,0.025, 0.01,0.005,0.0025, 0.001};
    // for the error of the numerical value at t=1
    
    double h = H[6]; // H[6] = 0.01
    double y0i = 1.0;
    double y1i = 0.0;
    double y2i = 0.0;
    double y3i = 0.0;
    double xi = 0.0;
    double xf = 30.0;

    const int D = static_cast<int>((30.0-0.0)/0.01);
    double y_array[D] = {};
    double x_array[D] = {};
    int steps = (xf-xi)/h;

    double y0 = y0i;
    double y1 = y1i;
    double y2 = y2i;
    double y3 = y3i;
    double x0 = xi;
    for(int i=0; i<steps; i++){
        y_array[i] = y0;
        x_array[i] = x0;
        y0 = y0_next(x0,y0,y1,y2,y3,h);
        y1 = y1_next(x0,y0,y1,y2,y3,h);
        y2 = y2_next(x0,y0,y1,y2,y3,h);
        y3 = y3_next(x0,y0,y1,y2,y3,h);
        x0 = x0+h;
    }
    //discrepency[j] = x_true - x4_0;
    cout << setprecision(9) << "x = " << x0 << ", y = " << y0 << "\n";
    // RK4 is accurate to 8 significant figures for values of h below but not equal to 0.005
    
    cout << "\n";
}



