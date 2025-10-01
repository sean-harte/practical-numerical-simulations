#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <cmath>
#include <iomanip> 
#include <math.h>

using namespace std;

double omega = 1.5; // test omega
double domega = 0.01; // amount i want to increase each mega by in loop
double omegai = 1.0; // initial mega for loop 
double omegaf = 2.0-domega; // maximum omega i want to test up to 


int n_iters = 10000; // no of iterations
double tol = 1e-6; // chosen tolerance
double maximum = 1.0; // initialising the maximum. just has to be larger than tolerance
double phi = 0; // initialising phi so it can be used in the loop

const double x_min = 0.0, x_max = 1.0;
const double y_min = 0.0, y_max = 1.0;

const int di = 70 + 1;
const double h = (x_max-x_min)/(di-1);

double tolVdphi_dphi[]={2.18484,2.119,1.89608,1.85547,1.84299,1.84319,1.84334,1.84325,1.84314,1.84313,1.84312,1.84312,1.84312}; // kept di = 11 and therefore h = 0.1 constant throughut this data collection
double tolVdphi_tol[]={0.1,0.05,0.01,0.005,0.001,0.0005,0.0001,5e-05,1e-05,5e-06,1e-06,5e-07,1e-07};// ^
double tolVdphi_omegas[]={1,1.18,1.38,1.4,1.45,1.44,1.4,1.4,1.41,1.42,1.43,1.43,1.43};// ^

double hVdphi_dphi[]={1.84312,1.75986,1.74099,1.73378,1.73007,1.72786,1.7264,1.72539,1.79599,1.7236,1.72329,1.72296,1.72276,1.72251,1.72235,1.72186}; // tolerance kept at 1e-6 throughout whole data collection
double hVdphi_h[]={0.1,0.05,0.0333333,0.025,0.02,0.0166667,0.0142857,0.0125,0.0111111,0.00909091,0.00833333,0.00769231,0.00714286,0.00666667,0.00625,0.005}; // ^ 1/h +1 = di
double hVdphi_omegas[]={1.43,1.66,1.76,1.81,1.85,1.87,1.89,1.9,1.91,1.93,1.93,1.94,1.94,1.95,1.95,1.96}; // ^

double grid[di][di]; // square grid 

// creating function that draws a grid between bounds. if the x bounds or y bounds are the same it will draw a line of constant value 
double lines(int start_x, int end_x, int start_y, int end_y, double value){ // note that all start and end values are the index. not the actual value on the grid. do that conversion yourself
    for(int i=start_x;i<end_x+1;i++){
        for(int j=start_y;j<end_y+1;j++){
            grid[i][j]=value;
        }
    }
return 0;
}

// creating function that draws a hollow rectangle between bounds
double const_rectangle(int start_x, int end_x, int start_y, int end_y, double value){
    for(int j=start_y;j<end_y+1;j++){
        grid[start_x][j]=value;
        grid[end_x][j]=value;
    }
    
    for(int i=start_x;i<end_x+1;i++){
        grid[i][start_y]=value;
        grid[i][end_y]=value;
    }
return 0;
}

// function concerts grid value to index that can be used as i or j in loops
int index(double value){
    return value*(di-1) ;
}

// calculates phi based on algorythem we are using 
double phi_next(int x_index, int y_index, double omega){ // using the SOR method. this is set up so that we have to traverse the plane left to right, bottom to top in order for this to update peoperly. grid[x_index+1][y_index]+grid[x_index][y_index+1] term should still be in the kth iteration while grid[x_index-1][y_index]+grid[x_index][y_index-1] should have already been updated.
    return (1.0-omega)*grid[x_index][y_index] + (omega/4.0)*(grid[x_index-1][y_index]+grid[x_index][y_index-1]+grid[x_index+1][y_index]+grid[x_index][y_index+1]);
}

// calculates 2nd derivative of x
double d2dx(int x_index, int y_index){ // 2nd derivative of x function
    return (1.0/pow(h,2))*(grid[x_index+1][y_index]+grid[x_index-1][y_index]-2*grid[x_index][y_index]);
}

// calculates 2nd derivative of y
double d2dy(int x_index, int y_index){ // 2nd derivative of y function
    return (1.0/pow(h,2))*(grid[x_index][y_index+1]+grid[x_index][y_index-1]-2*grid[x_index][y_index]);
}

// calculates the 1st derivative of x (central difference)
double ddx(int x_index, int y_index){ // 2nd derivative of x function
    return (0.5/h)*(grid[x_index+1][y_index]-grid[x_index-1][y_index]);
}

// calculates the 1st derivative of y (central difference)
double ddy(int x_index, int y_index){ // 2nd derivative of x function
    return (0.5/h)*(grid[x_index][y_index+1]-grid[x_index][y_index-1]);
}


// A boundary conditions
int Ax1 = index(0.2), Ax2 = index(0.4);
int Ay1 = index(0.7), Ay2 = index(0.9);
double Aval = 1.0;
// B boundary conditions
int Bx1 = index(0.8), Bx2 = index(0.8);
int By1 = index(0.1), By2 = index(0.6);
double Bval = 0.0;

// initialises grid based on boundary conditions. resetes grid to default basically
double initialise_grid(){
    // filling the grid with zeros and then boundary points 
    // filling with zeros
    for(int i=0;i<di;i++){
        for(int j=0;j<di;j++){
            grid[i][j] = 0;
            }
    }

    // make the A boundary
    const_rectangle(Ax1,Ax2,Ay1,Ay2,Aval);
    // make B boundary
    lines(Bx1,Bx2,By1,By2,Bval);

    // grid boundary (only need to do top and right side as left and bottom ar both 0 which is taken care of in the initialisation)
    for(int i=0;i<di;i++){
        grid[i][di-1] = float(i)/(di-1);
        grid[di-1][i] = float(i)/(di-1);
    }

return 0;
}



ofstream dataFile1("PNS_A3_q1_1.txt");
int main(){
// ========================================================================================================================================
// code for finding best omega value 

    initialise_grid();

    // for(int j=0;j<di;j++){ // printing initialised grid  
    //     cout << "\n";
    //     for(int i=0;i<di;i++){
    //         cout << "grid[" << i << "]" << "[" << j << "] = " << grid[i][j] << "\n" ;
    //         }
    // }

    cout << "=============================================================================" << "\n";
    double n_counter = 0;
    double best_omega=omegai;
    int best_iters = 1000;
    
    for(omega=omegai; omega<omegaf; omega+=domega){

        n_counter = 0;
        maximum = 1.0;
        
        initialise_grid();

        while(maximum>tol){ // loop for maximum difference of any point in the grid being less than chosen tolerance 
            double max_d = 0.0;
        // for(int n=0;n<n_iters;n++){ // loop for number of iterations
            for(int j=1;j<(di-1);j++){ // top to bottom (y-dir). starting at 1 and finsihing at (di-1) in both loops as to avoid changing the boundary points of the grid
                for(int i=1;i<(di-1);i++){// left to right (x-dir)
                    phi = grid[i][j];
                    if((Ax1==i || i==Ax2) && (Ay1<=j && j<=Ay2)){ // condition not to change x-line of A (A condition is hollow here i.e cannot change boundary of A but can change inside)
                        grid[i][j] = grid[i][j];
                        // cout << "1(x,y) = " << i << ", " << j << "\n";
                    } else if((Ax1<=i && i<=Ax2) && (Ay1==j || j==Ay2)){ // condition not to change y-line of A (A condition is hollow here i.e cannot change boundary of A but can change inside)
                        grid[i][j] = grid[i][j];
                        // cout << "1(x,y) = " << i << ", " << j << "\n";
                    } else if((Bx1<=i && i<=Bx2) && (By1<=j && j<=By2)){ // condition not to change B
                        grid[i][j] = grid[i][j];
                        // cout << "1(x,y) = " << i << ", " << j << "\n";
                    } else{
                        grid[i][j] = phi_next(i,j,omega);
                    }
                    max_d = max(max_d, abs(phi - grid[i][j])); 
                }
            }
            maximum = max_d;
            n_counter += 1;
        }

        // cout << "omega = " << omega << ", iterations = " << n_counter << " , tolerance = " << tol << "\n";
        dataFile1 << omega << " " << n_counter << " " << tol << "\n";
        if(best_iters>n_counter){
            best_iters = n_counter;
            best_omega = omega;
        }
    }
    dataFile1.close();

    ofstream dataFile2("PNS_A3_q1_2.txt");
    dataFile2 << best_omega << " " << best_iters << "\n";
    dataFile2.close();

    cout << "=============================================================================" << "\n";
    cout << "Best Omega = " << best_omega << "\n";
    cout << "=============================================================================" << "\n";

    stringstream ss;
    ss << "gnuplot -p -e \"set xlabel 'omega'; "
       << "set ylabel 'iterations'; "
       << "plot 'PNS_A3_q1_1.txt' using 1:2 with lines title 'Plot of iteration number for "
    //    << tol
       << di 
       << "x"
       << di 
       << " Grid',"
       << "'PNS_A3_q1_2.txt' using 1:2 with points linecolor rgb 'blue' pointtype 7 title 'Best Omega = "
       << best_omega
       << "'\"";

    string command = ss.str();

    system(command.c_str());


// ========================================================================================================================================
// part for computing dphi/dy(2/5,1/2) with good choice of omega

    double tol1 = tol; //1e-5;
    n_counter = 0;
    maximum = 10.0;

    // running code for best omega
    omega = best_omega;
    initialise_grid();
    while(maximum>tol1){ // loop for maximum difference of any point in the grid being less than chosen tolerance 
        double max_d = 0.0;
    // for(int n=0;n<n_iters;n++){ // loop for number of iterations
        for(int j=1;j<(di-1);j++){ // top to bottom (y-dir). starting at 1 and finsihing at (di-1) in both loops as to avoid changing the boundary points of the grid
            for(int i=1;i<(di-1);i++){// left to right (x-dir)
                phi = grid[i][j];
                if((Ax1==i || i==Ax2) && (Ay1<=j && j<=Ay2)){ // condition not to change x-line of A (A condition is hollow here i.e cannot change boundary of A but can change inside)
                    grid[i][j] = grid[i][j];
                    // cout << "1(x,y) = " << i << ", " << j << "\n";
                } else if((Ax1<=i && i<=Ax2) && (Ay1==j || j==Ay2)){ // condition not to change y-line of A (A condition is hollow here i.e cannot change boundary of A but can change inside)
                    grid[i][j] = grid[i][j];
                    // cout << "1(x,y) = " << i << ", " << j << "\n";
                } else if((Bx1<=i && i<=Bx2) && (By1<=j && j<=By2)){ // condition not to change B
                    grid[i][j] = grid[i][j];
                    // cout << "1(x,y) = " << i << ", " << j << "\n";
                } else{
                    grid[i][j] = phi_next(i,j,omega);
                }
                max_d = max(max_d, abs(phi - grid[i][j])); 
            }
        }
        maximum = max_d;
        n_counter += 1;
    }




    cout << "dphi/dy(2/5,1/2) = " << ddy(index(2.0/5.0),index(1.0/2.0)) << "\nwith " << "\nomega = " << omega << "\nh = " << h << "\ntolerance = " << tol1 << "\niterations = " << n_counter << "\nno. grid points = " << di << "x" << di <<  "\n";

    // for(int i=0;i<(di);i++){ // printing final grid 
    //     cout << "\n";
    //     for(int j=0;j<(di);j++){
    //         cout << "grid[" << i << "]" << "[" << j << "] = " << grid[i][j] << "\n" ;
    //         }
    // }

    cout << "=============================================================================" << "\n";

    // ========================================================================================================================================
    // trying to plot 3d grid 
    // set up dataset
    ofstream dataFile3("PNS_A3_q1_3.txt");
    for(int j=0;j<di;j++){
        for(int i=0;i<di;i++){
            dataFile3 << h*i << " " << h*j << " " << grid[i][j] << "\n";
        }
    }
    dataFile3.close();

    ofstream dataFile4("PNS_A3_q1_4.txt");
    for(int i=0;i<di;i++){
        for(int j=0;j<di;j++){
            dataFile4 << h*i << " " << h*j << " " << grid[i][j] << "\n";
        }
    }
    dataFile4.close();

    stringstream ss1;
    
    // Construct the gnuplot command for 3D surface plot from data file
    ss1 << "gnuplot -p -e \"set xlabel 'X-axis'; "
       << "set ylabel 'Y-axis'; "
       << "set zlabel 'Phi Value'; "
       << "set title 'Phi Grid for "
       << di 
       << "x"
       << di 
       << " Points';"
       << "splot 'PNS_A3_q1_3.txt' using 1:2:3 with points title 'Phi', "
    //    << "'PNS_A3_q1_4.txt' using 1:2:3 with points title 'Phi'; "
       << "\"";

    // Convert the stringstream to a string
    string command1 = ss1.str();

    // Output the command (for debugging, if necessary)
    cout << "Generated Command: " << command1 << endl;

    // Execute the gnuplot command
    system(command1.c_str());

    // ========================================================================================================================================
    // unloading data from arrays. tolerance vs dphi

    int length1 = sizeof(tolVdphi_dphi) / sizeof(tolVdphi_dphi[0]);

    ofstream dataFile5("PNS_A3_q1_5.txt");
    for(int i=0;i<length1;i++){
        dataFile5 << tolVdphi_dphi[i] << " " << log10(tolVdphi_tol[i]) << " " << tolVdphi_omegas[i] << " " << abs(tolVdphi_dphi[i]-tolVdphi_dphi[i+1]) << "\n";
    }
    dataFile5.close();

    //plotting log of tolerance vs omega
    stringstream ss2;
    ss2 << "gnuplot -p -e \"set xlabel 'log of Tolerance'; "
       << "set ylabel 'omega'; "
       << "plot 'PNS_A3_q1_5.txt' using 2:3 with linespoints title 'log of Tolerance v omega for "
    //    << tol
       << "11" 
       << "x"
       << "11" 
       << " Grid';"
    //    << "'PNS_A3_q1_2.txt' using 1:2 with points linecolor rgb 'blue' pointtype 7 title 'Best Omega = "
    //    << best_omega
       << "\"";

    string command2 = ss2.str();

    system(command2.c_str());

    // plotting log of tolerance vs solution of dphi
    stringstream ss3;
    ss3 << "gnuplot -p -e \"set xlabel 'log of Tolerance'; "
       << "set ylabel 'dphi/dy'; "
       << "plot 'PNS_A3_q1_5.txt' using 2:1 with linespoints title 'log of Tolerance vs dphi/dy for "
    //    << tol
       << "11" 
       << "x"
       << "11" 
       << " Grid';"
    //    << "'PNS_A3_q1_2.txt' using 1:2 with points linecolor rgb 'blue' pointtype 7 title 'Best Omega = "
    //    << best_omega
       << "\"";

    string command3 = ss3.str();

    system(command3.c_str());

    // plotting tolerance vs error 
    ofstream dataFile51("PNS_A3_q1_51.txt");
    for(int i=0;i<length1;i++){
        dataFile51 << log10(tolVdphi_tol[i]) << " " << log10(abs(tolVdphi_dphi[i]-tolVdphi_dphi[length1-1])) << "\n";
    }
    dataFile51.close();

    stringstream ss4;
    ss4 << "gnuplot -p -e \"set xlabel 'log of Tolerance'; "
       << "set ylabel 'log of Error'; "
       << "plot 'PNS_A3_q1_51.txt' using 1:2 with linespoints title 'log of Tolerance vs log of Error for "
    //    << tol
       << "11" 
       << "x"
       << "11" 
       << " Grid';"
    //    << "'PNS_A3_q1_2.txt' using 1:2 with points linecolor rgb 'blue' pointtype 7 title 'Best Omega = "
    //    << best_omega
       << "\"";

    string command4 = ss4.str();

    system(command4.c_str());

    // ========================================================================================================================================
    // unloading data from arrays. h vs dphi

    int length2 = sizeof(hVdphi_dphi) / sizeof(hVdphi_dphi[0]);

    ofstream dataFile6("PNS_A3_q1_6.txt");
    for(int i=0;i<length2;i++){
        dataFile6 << hVdphi_dphi[i] << " " << log10(hVdphi_h[i]) << " " << hVdphi_omegas[i] << "\n";
    }
    dataFile6.close();


    //plotting log of h vs omega
    stringstream ss5;
    ss5 << "gnuplot -p -e \"set xlabel 'log of h'; "
       << "set ylabel 'omega'; "
       << "plot 'PNS_A3_q1_6.txt' using 2:3 with linespoints title 'log of h v omega for "
       << "tolerance = 1e-6'"
    //    << "11" 
    //    << "x"
    //    << "11" 
    //    << " Grid';"
    //    << "'PNS_A3_q1_2.txt' using 1:2 with points linecolor rgb 'blue' pointtype 7 title 'Best Omega = "
    //    << best_omega
       << "\"";

    string command5 = ss5.str();

    system(command5.c_str());

    // plotting log of h vs solution of dphi
    stringstream ss6;
    ss6 << "gnuplot -p -e \"set xlabel 'log of h'; "
       << "set ylabel 'dphi/dy'; "
       << "plot 'PNS_A3_q1_6.txt' using 2:1 with linespoints title 'log of h vs dphi/dy for "
       << "tolerance = 1e-6'"
    //    << "11" 
    //    << "x"
    //    << "11" 
    //    << " Grid';"
    //    << "'PNS_A3_q1_2.txt' using 1:2 with points linecolor rgb 'blue' pointtype 7 title 'Best Omega = "
    //    << best_omega
       << "\"";

    string command6 = ss6.str();

    system(command6.c_str());

    // plotting tolerance vs error 
    ofstream dataFile61("PNS_A3_q1_61.txt");
    for(int i=0;i<length2;i++){
        dataFile61 << log10(hVdphi_h[i]) << " " << log10(abs(hVdphi_dphi[i]-hVdphi_dphi[length2-1])) << "\n";
    }
    dataFile61.close();

    stringstream ss7;
    ss7 << "gnuplot -p -e \"set xlabel 'log of h'; "
       << "set ylabel 'log of Error'; "
       << "plot 'PNS_A3_q1_61.txt' using 1:2 with linespoints title 'log of h vs log of Error for "
       << "tolerance = 1e-6'"
    //    << "11" 
    //    << "x"
    //    << "11" 
    //    << " Grid';"
    //    << "'PNS_A3_q1_2.txt' using 1:2 with points linecolor rgb 'blue' pointtype 7 title 'Best Omega = "
    //    << best_omega
       << "\"";

    string command7 = ss7.str();

    system(command7.c_str());

return 0;
}

// good to go
