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

// runge kutta below is for the third derivative using the fourth derivative as its slope. have this written as N in my notes. and function f is phi
// N1
double k1_y3(double x,double y,double y1,double y2,double y3, double h) {
    double k1_y3 = h*f(x,y,y1,y2,y3);
    return k1_y3;
}
// for the second derivative using the third derivative
// in notes as M
// M1
double k1_y2(double x,double y,double y1,double y2,double y3,double h) {
    double k1_y2 = h*y3; // already updated by runge-kutta
    return k1_y2;
}
// for the first derivative using the second derivative
// in notes as L
// L1
double k1_y1(double x,double y,double y1,double y2,double y3,double h) {
    double k1_y1 = h*y2; // already updated by runge-kutta
    return k1_y1;
}
// for y using the first derivative
// in notes as K
// K1
double k1_y0(double x,double y,double y1,double y2,double y3,double h) {
    double k1_y0 = h*y1; // already updated by runge-kutta
    return k1_y0;
}


// N2
double k2_y3(double x,double y,double y1,double y2,double y3,double h){
    double k2_y3 = h*f(x + h/2, y + (1.0/2.0)*k1_y0(x,y,y1,y2,y3,h), y1 + (1.0/2.0)*k1_y1(x,y,y1,y2,y3,h), y2 + (1.0/2.0)*k1_y2(x,y,y1,y2,y3,h), y3 + (1.0/2.0)*k1_y3(x,y,y1,y2,y3,h));
    return k2_y3;
}
//M2
double k2_y2(double x,double y,double y1,double y2,double y3,double h) {
    double k2_y2 = h*(y3 + (1.0/2.0)*k1_y3(x,y,y1,y2,y3,h)); 
    return k2_y2;
}
//L2
double k2_y1(double x,double y,double y1,double y2,double y3,double h) {
    double k2_y1 = h*(y2 + (1.0/2.0)*k1_y2(x,y,y1,y2,y3,h)); // already updated by runge-kutta
    return k2_y1;
}
//K2
double k2_y0(double x,double y,double y1,double y2,double y3,double h) {
    double k2_y0 = h*(y1 + (1.0/2.0)*k1_y1(x,y,y1,y2,y3,h)); // already updated by runge-kutta
    return k2_y0;
}


//N3
double k3_y3(double x,double y,double y1,double y2,double y3,double h){
    double k3_y3 = h*f(x + h/2, y + (1.0/2.0)*k2_y0(x,y,y1,y2,y3,h), y1 + (1.0/2.0)*k2_y1(x,y,y1,y2,y3,h), y2 + (1.0/2.0)*k2_y2(x,y,y1,y2,y3,h), y3 + (1.0/2.0)*k2_y3(x,y,y1,y2,y3,h));
    return k3_y3;
}
//M3
double k3_y2(double x,double y,double y1,double y2,double y3,double h) {
    double k3_y2 = h*(y3 + (1.0/2.0)*k2_y3(x,y,y1,y2,y3,h)); 
    return k3_y2;
}
//L3
double k3_y1(double x,double y,double y1,double y2,double y3,double h) {
    double k3_y1 = h*(y2 + (1.0/2.0)*k2_y2(x,y,y1,y2,y3,h)); // already updated by runge-kutta
    return k3_y1;
}
//K3
double k3_y0(double x,double y,double y1,double y2,double y3,double h) {
    double k3_y0 = h*(y1 + (1.0/2.0)*k2_y1(x,y,y1,y2,y3,h)); // already updated by runge-kutta
    return k3_y0;
}


// N4
double k4_y3(double x,double y,double y1,double y2,double y3,double h){
    double k4_y3 = h*f(x + h, y + k3_y0(x,y,y1,y2,y3,h), y1 + k3_y1(x,y,y1,y2,y3,h), y2 + k3_y2(x,y,y1,y2,y3,h), y3 + k3_y3(x,y,y1,y2,y3,h));
    return k4_y3;
}
//M4
double k4_y2(double x,double y,double y1,double y2,double y3,double h) {
    double k4_y2 = h*(y3 + k3_y3(x,y,y1,y2,y3,h)); 
    return k4_y2;
}
//L4
double k4_y1(double x,double y,double y1,double y2,double y3,double h) {
    double k4_y1 = h*(y2 + k3_y2(x,y,y1,y2,y3,h)); // already updated by runge-kutta
    return k4_y1;
}
//K4
double k4_y0(double x,double y,double y1,double y2,double y3,double h) {
    double k4_y0 = h*(y1 + k3_y1(x,y,y1,y2,y3,h)); // already updated by runge-kutta
    return k4_y0;
}



double y3_next(double x,double y,double y1,double y2,double y3,double h){
    double y3_next = y3 + (1.0/6.0)*(k1_y3(x,y,y1,y2,y3,h) + 2*k2_y3(x,y,y1,y2,y3,h) + 2*k3_y3(x,y,y1,y2,y3,h) + k4_y3(x,y,y1,y2,y3,h));
    return y3_next;
}

double y2_next(double x,double y,double y1,double y2,double y3,double h){
    double y2_next = y2 + (1.0/6.0)*(k1_y2(x,y,y1,y2,y3,h) + 2*k2_y2(x,y,y1,y2,y3,h) + 2*k3_y2(x,y,y1,y2,y3,h) + k4_y2(x,y,y1,y2,y3,h));
    return y2_next;
}

double y1_next(double x,double y,double y1,double y2,double y3,double h){
    double y1_next = y1 + (1.0/6.0)*(k1_y1(x,y,y1,y2,y3,h) + 2*k2_y1(x,y,y1,y2,y3,h) + 2*k3_y1(x,y,y1,y2,y3,h) + k4_y1(x,y,y1,y2,y3,h));
    return y1_next;
}

double y0_next(double x,double y,double y1,double y2,double y3,double h){
    double y0_next = y + (1.0/6.0)*(k1_y0(x,y,y1,y2,y3,h) + 2*k2_y0(x,y,y1,y2,y3,h) + 2*k3_y0(x,y,y1,y2,y3,h) + k4_y0(x,y,y1,y2,y3,h));
    return y0_next;
}


int main(){
    const double h_fix = 0.000006;
    double h =  h_fix; 
    double y0i = 1.0;
    double y1i = 0.0;
    double y2i = 0.0;
    double y3i = 0.0;
    const double xi = 0.0;
    const double xf = 30.0;

    double y0 = y0i;
    double y1 = y1i;
    double y2 = y2i;
    double y3 = y3i;
    double x0 = xi;
    ofstream dataFile("PNS_RK_method_q2_1.txt");
    while(x0 <= xf){
        dataFile << x0 << " " << y0 << "\n";
        y0 = y0_next(x0,y0,y1,y2,y3,h);
        if (x0 == xf){
            cout << "x0 = " << x0 << ", y(x=30) = " << y0 << "\n";
        }
        y1 = y1_next(x0,y0,y1,y2,y3,h);
        y2 = y2_next(x0,y0,y1,y2,y3,h);
        y3 = y3_next(x0,y0,y1,y2,y3,h);
        x0 = x0+h;
    }
    dataFile.close();
    cout << setprecision(9) << "x = " << x0 << ", y = " << y0 << ", h = " << h << "\n";
    cout << "\n";

     // plotting function y using rk4 for 4th order ODE 
    system("gnuplot -p -e \"set xlabel 'x value'; "
    "set ylabel 'y'; "
    //"set title 'RK4 for 4th Order ODE'; "
    "plot 'PNS_RK_method_q2_1.txt' using 1:2 with lines title 'y using RK4'\""); 
}









