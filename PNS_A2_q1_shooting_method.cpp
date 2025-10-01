#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <iomanip>

using namespace std;

// need to create shooting method that obeys the second order ODE corresponding to the function f
double f(double x,double x1,double t,double h){
    return -(t*x*(t+2))/(2+pow(t,2)*pow(x,2));
}

// for the first derivative, x1, using the second derivative given by function f
// L1
double k1_x1(double x,double x1,double t,double h) {
    return h*f(x,x1,t,h); 
}
// for x using the first derivative
// K1
double k1_x0(double x,double x1,double t,double h) {
    return h*x1; 
}

//L2
double k2_x1(double x,double x1,double t,double h) {
    return h*f(x+(1.0/2.0)*k1_x0(x,x1,t,h), x1+(1.0/2.0)*k1_x1(x,x1,t,h), t+h/2, h); 
}
//K2
double k2_x0(double x,double x1,double t,double h) {
    return h*(x1 + (1.0/2.0)*k1_x0(x,x1,t,h)); 
}

//L3
double k3_x1(double x,double x1,double t,double h) {
    return h*f(x+(1.0/2.0)*k2_x0(x,x1,t,h), x1+(1.0/2.0)*k2_x1(x,x1,t,h),t+h/2,h); 
}
//K3
double k3_x0(double x,double x1,double t,double h) {
    return h*(x1 + (1.0/2.0)*k2_x0(x,x1,t,h)); 
}

//L4
double k4_x1(double x,double x1,double t,double h) {
    return h*f(x+k3_x0(x,x1,t,h), x1+k3_x1(x,x1,t,h), t+h, h); 
}
//K4
double k4_x0(double x,double x1,double t,double h) {
    return h*(x1 + k3_x1(x,x1,t,h)); 
}



double x1_next(double x,double x1,double t,double h){
    return x1 + (1.0/6.0)*(k1_x1(x,x1,t,h) + 2*k2_x1(x,x1,t,h) + 2*k3_x1(x,x1,t,h) + k4_x1(x,x1,t,h));
}

double x0_next(double x,double x1,double t,double h){
    return x + (1.0/6.0)*(k1_x0(x,x1,t,h) + 2*k2_x0(x,x1,t,h) + 2*k3_x0(x,x1,t,h) + k4_x0(x,x1,t,h));
}


double error(double a,double b){
    return (a-b);
}

// define hyperparameters
const double ti = 0;
const double tf = 10.0;

const double x0i = 3.0/4.0;
const double x0f = -1.0;

int main(){
    // cout << "\n========================================================================================================================" << 
    // "\nPart 1\n"; // this is just to see what will be plotted for a random initiall guess of x1
    // initial variables for first part
    int n_steps = 10000;
    double h = (tf-ti)/n_steps;

    double x1_guess = 1.0;
    double x1 = x1_guess;
    double x = x0i;
    double t = ti;

    double x0 = x0i;
    double x01 = x1_guess; // these are placeholders so that x0_next and x1_next can be evaluated at the same value for x and x1

    // // this for loop is for finding what x(t=10) is computed for x1(t=0)=x_guess
    // // not doing shooting method yet
    // // want to plot the graph for error vs initial guesses first
    // ofstream dataFile1("PNS_A2_q1_1.txt");
    // dataFile1 << t << " " << x << "\n";
    // for(int i=0;i<n_steps;i++){

    //     x = x0_next(x0,x01,t,h);
    //     x1 = x1_next(x0,x01,t,h); 
    //     t += h;
    //     dataFile1 << t << " " << x << "\n";
    //     x0 = x;
    //     x01 = x1;

    // }
    // dataFile1.close();
    // cout << "\nx(t=10) = " << x << "\nt = " << t << "\nh = " << h << "\nn_steps = " << n_steps <<"\n";

    // system("gnuplot -p -e \"set xlabel 't'; "
    // "set ylabel 'x'; "
    // //"set title 'RK4 for 2nd Order ODE'; "
    // "plot 'PNS_A2_q1_1.txt' using 1:2 with lines title 'x'\""); 

    //
    // Part 1
    //========================================================================================================================
    // Part 2
    //
    cout << "\n========================================================================================================================" << 
    "\nPart 2\n";

    // now i want to plot the error for a range of x1_guess values 
    // initial variables for second part
    double x1_guess_i = -3.5;
    double x1_guess_f = 3.5;
    double x1_g_step = 0.005;

    int n_steps1 = 10000;
    double h1 = (tf-ti)/n_steps1;

    double x1_g = x1_guess_i;

    ofstream dataFile2("PNS_A2_q1_2.txt");
    while(x1_g<=x1_guess_f){
        // reset all my variables to initial values
        double x = x0i;
        double t = ti;

        double x0 = x0i;
        double x01 = x1_g; // these are placeholders so that x0_next and x1_next can be evaluated at the same value for x and x1

        for(int i=0;i<n_steps1;i++){

            x = x0_next(x0,x01,t,h1);
            x1 = x1_next(x0,x01,t,h1);
            t += h1;
            
            x0 = x;
            x01 = x1;

        }
    // cout << "\nxf = " << x << "\nx1_g = " << x1_g << "\nError = " << error(x,x0f) << "\n";
    dataFile2 << x1_g << " " << error(x,x0f) << " " << 0 << "\n";
    x1_g += x1_g_step;
    }
    dataFile2.close();
    // cout << "\nx(t=10) = " << x << "\nt = " << t << "\nh = " << h << "\nn_steps = " << n_steps <<"\n";

    
    system("gnuplot -p -e \"set xlabel 'x1 guess'; "
    "set ylabel 'Error'; "
    //"set title 'RK4 for 2nd Order ODE'; "
    "plot 'PNS_A2_q1_2.txt' using 1:2 with lines title 'Error', ""'PNS_A2_q1_2.txt' using 1:3 with lines title 'Error=0'\""); 

    //
    // Part 2
    //========================================================================================================================
    // Part 3
    //
    cout << "\n========================================================================================================================" << 
    "\nPart 3\n\n";

    // actually doing shooting method now
    // need to find x(tf) for x1_g = A and x1_g = A+d where d is a small difference 
    // initial variables for third part
    double deltai = -3.5;
    double deltaf = 3.5; // deltai/f are the boundaries of where we are looking for solutions of x1
    double d = 0.1; // little delta of how much of a difference there is between the 2 x1 guesses we are comparing

    int n_steps2 = 10000;
    double h2 = (tf-ti)/n_steps2;

    double tolerance = 1e-5;

    double x1_g1 = deltai;
    double x1_g2 = x1_g1+d;

    vector<double> soln_array; // this array stores the values of the solutions when they are found
    
    while(x1_g2 < deltaf){
        // reset all my variables to initial values
        double x0_1 = x0i; // holds the x value for the 1st guess of x1
        double x0_2 = x0i; // holds the x value for the 2nd guess of x1

        double x1_1 = x1_g1; // holds x1 value for 1st guess of x1
        double x1_2 = x1_g2; // holds x1 value for 2nd guess of x1
        
        t = ti;

        // these are placeholders so that x0_next and x1_next can be evaluated at the same value for x and x1
        double x00_1 = x0i; // placeholder for x0_1 i.e. x for the 1st x1 guess
        double x00_2 = x0i; // placeholder for x0_2 i.e. x for the 2nd x1 guess

        double x11_1 = x1_1; // placeholder for for 1st guess of x1
        double x11_2 = x1_2; // placeholder for 2nd guess of x1 

        for(int i=0;i<n_steps2;i++){

            x0_1 = x0_next(x00_1,x11_1,t,h2); // for 1st x1 guess
            x1_1 = x1_next(x00_1,x11_1,t,h2); // for 1st x1 guess

            x0_2 = x0_next(x00_2,x11_2,t,h2); // for 2nd x1 guess
            x1_2 = x1_next(x00_2,x11_2,t,h2); // for 2nd x1 guess

            t += h2; // for both 
            
            x00_1 = x0_1; 
            x00_2 = x0_2;

            x11_1 = x1_1;
            x11_2 = x1_2;

        }
        // cout << "\nPart 3" << 
        // "\nx(t=10) = " << x0_1 << ", " << x0_2 << 
        // "\nError = " << error(x0_1,x0f) << ", " << error(x0_2,x0f) << 
        // "\nx1 guess = " << x1_g1 << ", " << x1_g2 <<
        // "\nn_steps = " << n_steps2 << 
        // "\nt = " << t << "\n";

        // now i need a function that takes in the errors and tells the loop what do do with the initial giesses for x1
        if(error(x0_1,x0f)*error(x0_2,x0f) > 0){ // this is for the case when both solutions are on the same side of the error v x1_guess graph. need to be shifted on 
            x1_g1 = x1_g2;
            x1_g2 = x1_g1 + d;
            //cout << "\nhere1" << "\nx1_g1 = " << x1_g1 << "\nx1_g2 = " << x1_g2 << "\n"; 
        }

        if((error(x0_1,x0f)*error(x0_2,x0f) < 0) && (abs(error(x0_1,x0f)) > tolerance)){ // this is for if they are on opposite sides of the axis and have not satisfied the tolerance yet
            x1_g2 = (x1_g1+x1_g2)/2.0;
            //cout << "\nhere2" << "\n";
        }
        
        if(abs(error(x0_1,x0f)) < tolerance){ 
            cout << "x1 solution = " << setprecision(7) << x1_g1 << "\n";
            soln_array.push_back(x1_g1); 

            x1_g1 = x1_g2; // this moves the method on to find the next solution
            x1_g2 += d;
        }
        
    }
    cout << "\nTolerance = " << tolerance << "\n";

    //
    // Part 3
    //========================================================================================================================
    // Part 4
    //
    cout << "\n========================================================================================================================" << 
    "\nPart 4\n";

    // plotting the solutions found in last part of code

    //cout << "\n" << soln_array[0] << ", " << soln_array[1] << "\n";

    // initial variables for first part
    int ind = 0; // index for the soln_array
    int n_steps3 = 10000;
    double h3 = (tf-ti)/n_steps3;

    x1 = soln_array[ind];
    x = x0i;
    t = ti;

    x0 = x0i;
    x01 = soln_array[ind]; // these are placeholders so that x0_next and x1_next can be evaluated at the same value for x and x1

    ofstream dataFile3("PNS_A2_q1_3.txt");
    dataFile3 << t << " " << x << "\n";
    for(int i=0;i<n_steps3;i++){

        x = x0_next(x0,x01,t,h3);
        x1 = x1_next(x0,x01,t,h3); 
        t += h3;
        dataFile3 << t << " " << x << "\n";
        x0 = x;
        x01 = x1;

    }
    dataFile3.close();
    cout << "\nx(t=10) = " << x << "\nt = " << t << "\nh = " << h3 << "\nn_steps = " << n_steps3 <<"\n";



    stringstream ss;
    ss << "gnuplot -p -e \"set xlabel 't'; "
       << "set ylabel 'x'; "
       << "plot 'PNS_A2_q1_3.txt' using 1:2 with lines title 'Plot for dx/dt = " 
       << soln_array[ind] << "'\"";

    string command = ss.str();

    system(command.c_str());

    



}






