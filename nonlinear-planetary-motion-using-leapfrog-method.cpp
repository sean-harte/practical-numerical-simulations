#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>

using namespace std;

// big G is set to 1
const int G = 1;
// 2d system
const int dim = 2; 

// defining initial conditions for the planets and defining the arrays that i will use like vectors for all the x and v terms
/// planet 0 ///
double m0=2.2, x0xi=-0.5, x0yi=0.1, v0xi=-0.84, v0yi=0.65;
double x0[dim] = {x0xi,x0yi};
double v0[dim] = {v0xi,v0yi};

/// planet 1 ///
double m1=0.8, x1xi=-0.6, x1yi=-0.2, v1xi=1.86, v1yi=0.7;
double x1[dim] = {x1xi,x1yi};
double v1[dim] = {v1xi,v1yi};

/// planet 2 ///
double m2=0.9, x2xi=0.5, x2yi=0.1, v2xi=-0.44, v2yi=-1.5;
double x2[dim] = {x2xi,x2yi};
double v2[dim] = {v2xi,v2yi};

/// planet 3 ///
double m3=0.4, x3xi=0.5, x3yi=0.4, v3xi=1.15, v3yi=-1.6;
double x3[dim] = {x3xi,x3yi};
double v3[dim] = {v3xi,v3yi};

// creating a function that calculates the distance / radius between 2 points expressed as 2 2d arrays
double radius(double arr1[dim], double arr2[dim]){
    double squared_sum = pow((arr1[0]-arr2[0]),2) + pow((arr1[1]-arr2[1]),2);
    double rad = sqrt(squared_sum);
    return rad;
}

// creating function for an r_x for a force acted between 2 planets so that I know the x-direction that the force is acting in
double r_x(double arr1[dim], double arr2[dim]){ // arr1 is the jth planet arr2 is the ith planet (always assuming jth is the planet of interest)
    return arr2[0]-arr1[0];
}
// creating function for an r_y for a force acted between 2 planets so that I know the y-direction that the force is acting in
double r_y(double arr1[dim], double arr2[dim]){ // arr1 is the jth planet arr2 is the ith planet (always assuming jth is the planet of interest)
    return arr2[1]-arr1[1];
}

// now creating the normalised r_hatx and r_haty for the force between two planets
double r_hatx(double arr1[dim], double arr2[dim]){
    return r_x(arr1,arr2)/radius(arr1,arr2);
}
double r_haty(double arr1[dim], double arr2[dim]){
    return r_y(arr1,arr2)/radius(arr1,arr2);
}

// now i am creating a function that calculates the force experienced by the jth planet by the other 3 (i) in the x and y directions
/// first array passed will be the jth planet, other 3 are used to calculate force
double force_x(double arr_xj[dim], double massj, double arr_xi1[dim], double massi1, double arr_xi2[dim], double massi2, double arr_xi3[dim], double massi3){
    double fx01 = ((G*massj*massi1)/(pow(radius(arr_xj,arr_xi1),2)))*r_hatx(arr_xj,arr_xi1);
    double fx02 = ((G*massj*massi2)/(pow(radius(arr_xj,arr_xi2),2)))*r_hatx(arr_xj,arr_xi2);
    double fx03 = ((G*massj*massi3)/(pow(radius(arr_xj,arr_xi3),2)))*r_hatx(arr_xj,arr_xi3);
    return fx01+fx02+fx03;
}
double force_y(double arr_xj[dim], double massj, double arr_xi1[dim], double massi1, double arr_xi2[dim], double massi2, double arr_xi3[dim], double massi3){
    double fy01 = ((G*massj*massi1)/(pow(radius(arr_xj,arr_xi1),2)))*r_haty(arr_xj,arr_xi1);
    double fy02 = ((G*massj*massi2)/(pow(radius(arr_xj,arr_xi2),2)))*r_haty(arr_xj,arr_xi2);
    double fy03 = ((G*massj*massi3)/(pow(radius(arr_xj,arr_xi3),2)))*r_haty(arr_xj,arr_xi3);
    return fy01+fy02+fy03;
}



int main(){
    /// now I will implememt the leap frog algorithm by just changing the existing arrays / vectors corresponding to the position and velocity
    double ti = 0.0;
    double tf = 5.0;
    double t = ti;
    int n_steps = 100000;
    double h = (tf-ti)/n_steps;

    // to get initial positions 
    ofstream dataFile2("PNS_A2_q2_2.txt");
    dataFile2 << t << " " << x0[0] << " " << x0[1] << " " << x1[0] << " " << x1[1] << " " << x2[0] << " " << x2[1] << " " << x3[0] << " " << x3[1] << "\n";
    dataFile2.close();

    // starting data collection for the graph
    ofstream dataFile1("PNS_A2_q2_1.txt");
    dataFile1 << t << " " << x0[0] << " " << x0[1] << " " << x1[0] << " " << x1[1] << " " << x2[0] << " " << x2[1] << " " << x3[0] << " " << x3[1] << "\n";
    
    for(int i=0;i<n_steps;i++){
        // x_{k+1/2} :: first do step for all planets position moving forward by half a step so they are all at the same time 
        //========================================================
        //// 0th planet ////
        x0[0] += (h/2)*v0[0]; 
        x0[1] += (h/2)*v0[1]; 

        //// 1st planet ////
        x1[0] += (h/2)*v1[0]; 
        x1[1] += (h/2)*v1[1]; 

        //// 2nd planet ////
        x2[0] += (h/2)*v2[0]; 
        x2[1] += (h/2)*v2[1]; 

        //// 3rd planet ////
        x3[0] += (h/2)*v3[0]; 
        x3[1] += (h/2)*v3[1]; 
        //========================================================

        // v_{k+1} :: now get next step of v for all the planets which takes a full step and uses the force equations i have made
        //========================================================
        //// 0th planet ///
        v0[0] += (h/m0)*force_x(x0,m0,x1,m1,x2,m2,x3,m3);
        v0[1] += (h/m0)*force_y(x0,m0,x1,m1,x2,m2,x3,m3);

        //// 1st planet ////
        v1[0] += (h/m1)*force_x(x1,m1,x0,m0,x2,m2,x3,m3);
        v1[1] += (h/m1)*force_y(x1,m1,x0,m0,x2,m2,x3,m3);

        //// 2nd planet ////
        v2[0] += (h/m2)*force_x(x2,m2,x0,m0,x1,m1,x3,m3);
        v2[1] += (h/m2)*force_y(x2,m2,x0,m0,x1,m1,x3,m3);

        //// 3rd planet ////
        v3[0] += (h/m3)*force_x(x3,m3,x0,m0,x1,m1,x2,m2);
        v3[1] += (h/m3)*force_y(x3,m3,x0,m0,x1,m1,x2,m2);
        //========================================================

        // x_{k+1} :: now go the other half a step for the position vector to get the full step
        //========================================================
        //// 0th planet ////
        x0[0] += (h/2)*v0[0]; 
        x0[1] += (h/2)*v0[1]; 

        //// 1st planet ////
        x1[0] += (h/2)*v1[0]; 
        x1[1] += (h/2)*v1[1]; 

        //// 2nd planet ////
        x2[0] += (h/2)*v2[0]; 
        x2[1] += (h/2)*v2[1]; 

        //// 3rd planet ////
        x3[0] += (h/2)*v3[0]; 
        x3[1] += (h/2)*v3[1]; 
        //========================================================

        // set time forward by step h
        t += h;

        // add to datafile
        dataFile1 << t << " " << x0[0] << " " << x0[1] << " " << x1[0] << " " << x1[1] << " " << x2[0] << " " << x2[1] << " " << x3[0] << " " << x3[1] << "\n";
    }
dataFile1.close();

// to get final positions 
    ofstream dataFile3("PNS_A2_q2_3.txt");
    dataFile3 << t << " " << x0[0] << " " << x0[1] << " " << x1[0] << " " << x1[1] << " " << x2[0] << " " << x2[1] << " " << x3[0] << " " << x3[1] << "\n";
    dataFile3.close();

cout << "t,h,n_steps = " << t << ", " << h << ", " << n_steps << "\nx0 = (" << x0[0] << ", " << x0[1] << ")" << "\nx1 = (" << x1[0] << ", " << x1[1] << ")" << "\nx2 = (" << x2[0] << ", " << x2[1] << ")" << "\nx3 = (" << x3[0] << ", " << x3[1] << ")" << "\n";

stringstream ss;
double x_min = -30.0;
double x_max = 30.0;
double y_min = -6.0;
double y_max = 5.0;
    ss << "gnuplot -p -e \"set xlabel 'x'; "
       << "set ylabel 'y'; "
    //    << "set xrange [" << x_min << ":" << x_max << "]; "
    //    << "set yrange [" << y_min << ":" << y_max << "]; "
       << "plot 'PNS_A2_q2_1.txt' using 2:3 with lines title 'Planet 0', "
       << "'PNS_A2_q2_1.txt' using 4:5 with lines title 'Planet 1', "
       << "'PNS_A2_q2_1.txt' using 6:7 with lines title 'Planet 2', "
       << "'PNS_A2_q2_1.txt' using 8:9 with lines title 'Planet 3', "

       << "'PNS_A2_q2_2.txt' using 2:3 with points linecolor rgb 'purple' pointtype 7 title 'Planet 0 at t=0', "
       << "'PNS_A2_q2_2.txt' using 4:5 with points linecolor rgb 'green' pointtype 7 title 'Planet 1 at t=0', "
       << "'PNS_A2_q2_2.txt' using 6:7 with points linecolor rgb 'blue' pointtype 7 title 'Planet 2 at t=0', "
       << "'PNS_A2_q2_2.txt' using 8:9 with points linecolor rgb 'orange' pointtype 7 title 'Planet 3 at t=0', "

       << "'PNS_A2_q2_3.txt' using 2:3 with points linecolor rgb 'purple' pointtype 1 title 'Planet 0 at t=5', "
       << "'PNS_A2_q2_3.txt' using 4:5 with points linecolor rgb 'green' pointtype 1 title 'Planet 1 at t=5', "
       << "'PNS_A2_q2_3.txt' using 6:7 with points linecolor rgb 'blue' pointtype 1 title 'Planet 2 at t=5', "
       << "'PNS_A2_q2_3.txt' using 8:9 with points linecolor rgb 'orange' pointtype 1 title 'Planet 3 at t=5'; "
       //<< soln_array[ind] 
       << "\"";

    string command = ss.str();

    system(command.c_str());

}


















// // now creating functions that take in all planets (x,y) coords and calculates the total r_x and r_y components. I will then be sure to normalise these later to get r_hat
// double r_x_tot(double arr1[dim], double arr2[dim], double arr3[dim], double arr4[dim]){ // arr1 is the jth planet
//     return r_x(arr1,arr2) + r_x(arr1,arr3) + r_x(arr1,arr4);
// }
// double r_y_tot(double arr1[dim], double arr2[dim], double arr3[dim], double arr4[dim]){ // arr1 is the jth planet
//     return r_y(arr1,arr2) + r_y(arr1,arr3) + r_y(arr1,arr4);
// }

// // now finally making a fully functional r_hatx and r_hat y that is normalised
// double r_hatx(double arr1[dim], double arr2[dim], double arr3[dim], double arr4[dim]){
//     // if(abs(r_x_tot(arr1,arr2,arr3,arr4))==0 && abs(r_y_tot(arr1,arr2,arr3,arr4))==0){
//     //     return 0;
//     // } else{
//     // return r_x_tot(arr1,arr2,arr3,arr4)*(1/sqrt(pow(r_x_tot(arr1,arr2,arr3,arr4),2) + pow(r_y_tot(arr1,arr2,arr3,arr4),2)));
//     // }
//     return r_x_tot(arr1,arr2,arr3,arr4)*(1/sqrt(pow(r_x_tot(arr1,arr2,arr3,arr4),2) + pow(r_y_tot(arr1,arr2,arr3,arr4),2)));
// }
// double r_haty(double arr1[dim], double arr2[dim], double arr3[dim], double arr4[dim]){
//     // if(abs(r_x_tot(arr1,arr2,arr3,arr4))==0 && abs(r_y_tot(arr1,arr2,arr3,arr4))==0){
//     //     return 0;
//     // } else{
//     // return r_y_tot(arr1,arr2,arr3,arr4)*(1/sqrt(pow(r_x_tot(arr1,arr2,arr3,arr4),2) + pow(r_y_tot(arr1,arr2,arr3,arr4),2)));
//     // }
//     return r_y_tot(arr1,arr2,arr3,arr4)*(1/sqrt(pow(r_x_tot(arr1,arr2,arr3,arr4),2) + pow(r_y_tot(arr1,arr2,arr3,arr4),2)));
// }

