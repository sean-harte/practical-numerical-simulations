#include <iostream>
#include <iomanip>
#include<cmath>
#include <fstream>
#include<cstdlib>

using namespace std;
const double ti = 0.0;
const double tf = 1.0;
const double accuracy = 0.0000001;

double analytic_fn(double t){
    double x = -1 + exp(pow(t,3)/3 - 2*pow(t,2) + 4*t); 
    return x;
}

// all doubles below are the functions that I will need in the for loop that i will iterate using t0 and y0
double f(double t,double x){ // this is the function for the first derivative or slope of x (dx/dt)
    double a = pow((t-2),2)*(x+1);
    return a;
}

// all ki terms are the different orders of the runge-kutta method. doing these seperately so that the 
// code is much easier to read and implement
// all these ki are for the 4th order runge-kutta method
double k1(double t, double x) {
    double k1 = f(t,x);
    return k1;
}

double k2(double t,double x,double h){
    double k2 = f(t + (h/2), x + h*(k1(t,x)/2));
    return k2;
}

double k3(double t,double x,double h){
    double k3 = f(t + (h/2), x + h*(k2(t,x,h)/2));
    return k3;
}

double k4(double t,double x,double h){
    double k4 = f(t + h, x + h*k3(t,x,h));
    return k4;
}

double x4_next(double t,double x, double h){
    double x4_next = x + (h/6)*(k1(t,x) + 2*k2(t,x,h) + 2*k3(t,x,h) + k4(t,x,h));
    return x4_next;
}

//
// using the shooting method to get the exact value of h that yields a result accurate to 8 significant figures
//
double h_g = 0.5; // guess value. this is arbitrary and only used when calling the function
double sf = accuracy; // we know the answer converges to 1<=x<10 from RK4. if answer converged to 10<=x<100 then sf = 1e-6
double a = 1; // starting value for error. this is arbitrary. 
double tol = 0.0001; // this is the tolerance for the ratio of error to desired accuracy. upper bound is 1, lower bound is 1-tol
double h_finder(double sf, double h_guess){ // h_finder finds the minimum value of h that gives an accuracy of between sf = 1e-7 and (1-tol)e-7. for an initial guess of h_g
    double h = h_guess;
    while((a/sf > 1) || (a/sf < 1 - tol)){
        // for the error of the numerical value at t=1
        double x_true = 0;
        double x4_0 = 0.0; // initial y term for the 4th order rk method
        double t0 = ti; // time does not need to be seperate for the different methods
        int steps = (tf-ti)/h;
        for(int i=0; i<steps; i++){
            x4_0 = x4_next(t0,x4_0,h);
            t0 = t0+h;
            x_true = analytic_fn(t0);
        }
    a = abs(x_true - x4_0);
    h = h/2.0;
    //cout << "Too High" << "\n";

    if(a/sf < 1 - tol){
        h = 4.5*h; // I am doing 3.9 as the factor because I want to get it back to almost its previous h but not quite. a factor of 4 would get it back to its previous h value
        // cout << "Too Low" << "\n";
    }
    //cout << "a/sf = " << a/sf << ", h = " << h << "\n";
    }
    cout << "Maximum h for accuracy of 8s.f. = " << h << ", Error/(1e-7) = " << a/sf << "\n";
    return 0;
}


int main(){
    const int len = 11;
    double discrepency[len] = {};
    double H[len] = {1.0,0.5,0.25, 0.1,0.05,0.025, 0.01,0.005,0.004,0.0025, 0.001}; //0.001
    // for the error of the numerical value at t=1
    for(int j=0; j<len; j++){
        double h = H[j];
        double x_true = 0;
        double x4_0 = 0.0; // initial y term for the 4th order rk method
        double t0 = ti; // time does not need to be seperate for the different methods
        int steps = (tf-ti)/h;
        for(int i=0; i<steps; i++){
            x4_0 = x4_next(t0,x4_0,h);
            t0 = t0+h;
            x_true = analytic_fn(t0);
        }
    discrepency[j] = x_true - x4_0;
    cout << setprecision(9) << "h = " << h << ", t = " << t0 << ", x4(t=1) = " << x4_0 << ", x_true(t=1) = " << x_true << ", Error = "  << abs(x4_0 - x_true) <<", steps = " << steps << "\n";
    // RK4 is accurate to 8 significant figures for values of h below but not equal to 0.005
    }
    cout << "\n";
    

    // for the plot of rk4 with h = 0.1
    const int N = 1+(tf-ti)/0.1;
    double x4_array[N] = {};
    double t_array[N] = {};
    double xanalytic_array[N] = {};
    const double h1 = H[3];
    double x4_01 = 0.0;
    double t01 = 0.0;
    for(int i=0;i<N;i++){
        t_array[i] = t01;
        x4_array[i] = x4_01;
        xanalytic_array[i] = analytic_fn(t01);
        x4_01 = x4_next(t01,x4_01,h1);
        t01 = t01 + h1;

    }

    ofstream dataFile("PNS_RK_method.txt");
    ofstream dataFile2("please.txt");
    ofstream functionFile("PNS_RK_func_file.txt");
    for(int i=0;i<len;i++){
        dataFile << H[i] << " " << discrepency[i] << " " << accuracy*0.1 << "\n";
    }
     dataFile.close();

    for(int i=0;i<N;i++){
        dataFile2 << t_array[i] << " " << x4_array[i] << " " << xanalytic_array[i] << " " << xanalytic_array[i] - x4_array[i] << "\n";
    }
    dataFile2.close();

    cout << h_finder(sf,h_g) << "\n";



    // plotting error of RK4 at t=1 for multiple step sizes
    system("gnuplot -p -e \"set xlabel 'h size'; "
    "set ylabel 'Error at t=1'; "
    "set logscale x;"
    "set logscale y;"
    //"set title 'Error between RK4 and Analytic Solution at t=1 for Different h Size'; "
    "plot 'PNS_RK_method.txt' using 1:2 with linespoints title 'Error', ""'PNS_RK_method.txt' using 1:3 with lines title 'Accuracy to 8 s.f.'\""); 

    // plotting RK4 for step size h = 0.1
    system("gnuplot -p -e \"set xlabel 't'; "
    "set ylabel 'x'; "
    //"set title 'Plot of RK4 for Step Size h = 0.1'; "
    "plot 'please.txt' using 1:3 with linespoints title 'Analytic Solution',""'please.txt' using 1:2 with linespoints title 'RK4'\""); 

    // plotting the error of the RK4 for step size h=0.1
    system("gnuplot -p -e \"set xlabel 't'; "
    "set ylabel 'Error'; "
    //"set title 'Error of RK4 for Step Size h = 0.1'; "
    "plot 'please.txt' using 1:4 with linespoints title 'Error'\""); 

    return 0;
}
