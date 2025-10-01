#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <cmath>
#include <iomanip> 
#include <math.h>
#include <random>
#include <algorithm> // For min function


using namespace std;

// dimension of grid
int l = 24;
int lx = l;
int ly = l;

// how many q values 
int q = 3;

// random number generator initialisation for q value 
random_device rd;
mt19937 gen(rd());
uniform_int_distribution<> distrib(1,q);

// rng for probability random number i.e. randm number between 0 and 1
random_device rd1;
mt19937 gen1(rd1());
uniform_real_distribution<> distrib1(0.0,1.0);

// function to make the grid
vector<vector<int> > grid;
double make_grid(int rows, int cols, int q){
    // initialise all inputs to 0
    grid.resize(rows);
    for(int i=0;i<rows;i++){
        grid[i].resize(cols);
    }

    // choose random q for each site
    int random_number = distrib(gen);
    for(int j=0;j<rows;j++){
        for(int i=0;i<cols;i++){
            grid[i][j] = distrib(gen);
        }
    }

    return 0;
}


// fn counting how many of a value 
int how_many(int a){
    int no = 0;
    for(int j=0;j<ly;j++){
        for(int i=0;i<ly;i++){
            no += (grid[i][j] == a);
        }
    }
    return no;
}

// fn for magnetisation
double mag(){
    double most = 0;
    double number = 0;
    for(int i=1;i<q+1;i++){
        number = how_many(i);
        if(number>most){
            most = number;
        }
    }
    double magnetisation = (q*(most/(lx*ly))-1)/(q-1);
    return magnetisation;
}

double mag_sqd(){
    return pow(mag(),2);
}

// function getting the change in action
double ds(int x_index, int y_index, int prop) {
    double s_old = 0.0;
    double s_new = 0.0;
    int site = grid[x_index][y_index];

    // int dx[] = {1,-1,0,0};
    // int dy[] = {0,0,1,-1};

    // for(int i=0;i<4;i++){
    //     int nx = (x_index+dx[i]+lx)%lx;
    //     int ny = (y_index+dy[i]+ly)%ly;
    //     s_old += (grid[nx][ny] != site);
    //     s_new += (grid[nx][ny] != prop);
    // }

    s_old += (grid[(x_index-1+lx)%lx][y_index]!=site) + (grid[(x_index+1+lx)%lx][y_index]!=site) + (grid[x_index][(y_index-1+ly)%ly]!=site) + (grid[x_index][(y_index+1+ly)%ly]!=site);
    s_new += (grid[(x_index-1+lx)%lx][y_index]!=prop) + (grid[(x_index+1+lx)%lx][y_index]!=prop) + (grid[x_index][(y_index-1+ly)%ly]!=prop) + (grid[x_index][(y_index+1+ly)%ly]!=prop);
    double delta_s = s_new-s_old;
    return delta_s;
}


// function for the acceptance probability
double p_acc(int x_index,int y_index, int prop, double beta){
    double prob = exp(-beta*ds(x_index,y_index,prop));
    return min(1.0,prob);
}

// function that carries out a sweep
double do_sweep(int no_sweeps, double beta){
    for(int n=0;n<no_sweeps;n++){
        for(int j=0;j<ly;j++){
            for(int i=0;i<lx;i++){
                double r = distrib1(gen1);
                int proposal = distrib(gen);
                if(p_acc(i,j,proposal,beta)>r){
                    grid[i][j] = proposal;
                }
            }
        }
    }    
    return 0;
}


double mean_func(vector<double> data){
    int N = data.size();

    double mean = 0.0;
    for(int i=0; i<N; i++){
        mean += data[i];
    }
    mean /= N;

    return mean;
}

double variance_func(vector<double> data, double mean, int tau_i){
    int N = data.size();
    double norm = 0.0;
    double variance = 0.0;
    for(int i=0; i<N; i+=tau_i){
        variance += pow((data[i] - mean),2);
        norm += 1.0;
    }
    variance /= norm;

    return variance;
}

// function that carries out autocorrelation calculation for some t (lag) value
double autoc(vector<double> vector1, double mean, int t){
    int N = vector1.size();
    double sum = 0.0;
    for(int i=0;i<(N-t);i++){
        sum += (vector1[i] - mean)*(vector1[i+t] - mean); // /variance
    }
    return (1.0/(N-t))*sum; 
    
}

double tau(vector<double> data, double mean, int w){ // this is the integration time. this is the periodicity at which we should take measurements in the metropolis aglorithem
    int N = data.size();
    double gamma_0 = autoc(data,mean,0);
    double sum1 = 0.0;
    for(int j=1;j<w;j++){
        sum1 += autoc(data,mean,j)/gamma_0;
    }
    return (1.0/2.0) + sum1;
}

int tau_int(vector<double> data, double mean, int w_range){ // finds the maximum value for tau_int from window sizes between 1 and w_range. worst case scenario tau_int is grater than it should be. not a bad thing. the adjacent measurements are still not correlated. This rounds up to the nearest integer.
    int highest_tau = 0;
    double gamma_0 = autoc(data,mean,0);
    for(int w=1;w<(w_range+1);w++){
        double sum = 0.5;
        for(int t=1;t<(w+1);t++){
            sum += autoc(data,mean,t)/gamma_0;
        }
        if(sum > highest_tau){
            highest_tau = sum;
        }
    }
    return highest_tau + 1;
}

// function that measures an observable using the integration time so that the measurement is more accurate
double measure(vector<double> data, int tau_i){
    int N = data.size();
    double norm = 0.0;
    double sum = 0.0;
    for(int i=0;i<N;i+=tau_i){
        sum += data[i];
        norm += 1.0;
    }
    return sum/norm;
}


int main(){

    make_grid(lx,ly,q);

    int max_sweeps = 1000; // amount of times the do_sweeps fn is called in the loop
    
    int nbeta = 50;
    double betai = 0.5;
    double betaf = 1.5;
    double dbeta = (betaf-betai)/(nbeta-1);

    // array for plotting magnetism as a function of beta
    vector<double> mag_vs_beta;
    mag_vs_beta.resize(nbeta+1);
    vector<double> beta_array;
    beta_array.resize(nbeta+1);
    vector<int> tau_int_array;
    tau_int_array.resize(nbeta+1);
    vector<double> mag_sqd_vs_beta;
    mag_sqd_vs_beta.resize(nbeta+1);
    vector<double> var_mag_vs_beta;
    var_mag_vs_beta.resize(nbeta+1);
    vector<double> var_mag_sqd_vs_beta;
    var_mag_sqd_vs_beta.resize(nbeta+1);

    // initialise arrays to store values to be used in analysis and autocorrelation functions
    // ============================================================
    vector<double> mag_array;
    mag_array.resize(max_sweeps);

    vector<double> mag_sqd_array;
    mag_sqd_array.resize(max_sweeps);

    vector<double> var_mag_array;
    var_mag_array.resize(max_sweeps);

    vector<double> var_mag_sqd_array;
    var_mag_sqd_array.resize(max_sweeps);
    // ============================================================


    double counter = 0.0;
    for(double beta1=betai;beta1<(betaf+dbeta);beta1+=dbeta){ // betai
        double beta = beta1;
        
        int init_sweeps = 10000;

        int no_sweeps = 1; // no of sweeps that the do_sweeps function does itself



        for(int i=0;i<max_sweeps;i++){
            
            do_sweep(no_sweeps,beta);
            
            double magnetism = mag();
            
            mag_array[i] = magnetism;

            double magnetism_sqd = pow(magnetism,2);
            
            mag_sqd_array[i] = magnetism_sqd;


        }

        int w_range = 100;
        double mag_mean = mean_func(mag_array);
        double mag_sqd_mean = mean_func(mag_sqd_array);
        int tau_i = max(tau_int(mag_array,mag_mean, w_range),tau_int(mag_sqd_array,mag_sqd_mean, w_range)); //taking the larger tau_int that was calculated for magnetism and mag_squared so that the same number of samples from both are used to calculate the average
        // int tau_i = 5;
        // double mag_variance = variance_func(mag_array,mag_mean,tau_i);
        // double mag_sqd_variance = variance_func(mag_sqd_array,mag_sqd_mean,tau_i);
        
        int w = 100;
    
        int t = 5;

        // cout << "===========================================================================" << "\n";
        // cout << "No. of Sweeps = " << max_sweeps << "\n";
        // cout << "Average Magnetisation = " << av_mag << "\n";
        // cout << "Average Magnetisation Squared = " << av_mag_sqd << "\n";
        // cout << "Variance of Magnetism = " << var_mag << "\n";
        // cout << "Variance of Magnetism Squared = " << var_mag_sqd << "\n";
        // cout << "Autocorrelation for t=" << t << " of Magnetism = " << autoc(mag_array,mag_mean,mag_variance,t) << "\n";
        // cout << "Tau (for w = " << w << ") = " << tau(mag_array,mag_mean,mag_variance, w) << "\n";
        // cout << "Tau_int (for w range = " << w_range << ") = " << tau_i << "\n";
        // cout << "Mag Measurement (for tau_int = " << tau_i << ") = " << measure(mag_array,tau_i) << "\n";
        // cout << "===========================================================================" << "\n";


        double true_mag_mean = measure(mag_array,tau_i);
        mag_vs_beta[counter] = true_mag_mean;

        double true_mag_sqd_mean = measure(mag_sqd_array,tau_i);
        mag_sqd_vs_beta[counter] = true_mag_sqd_mean;

        beta_array[counter] = beta;
        tau_int_array[counter] = tau_i;
        
        double mag_variance = variance_func(mag_array,true_mag_mean,tau_i);
        double mag_sqd_variance = variance_func(mag_sqd_array,true_mag_sqd_mean,tau_i);

        var_mag_vs_beta[counter] = mag_variance;
        var_mag_sqd_vs_beta[counter] = mag_sqd_variance;


        //===================================================================================================================================
        //===================================================================================================================================
        // FOR MAGNETISM
        //===================================================================================================================================
        //===================================================================================================================================
        // if statement for plots so that I dont get a plot for every beta
        if((counter == 0.0 || counter == 25.0 || counter == nbeta/1.0)){
            //===================================================================================================================================
            //// #### PLOT OF AUTOCORRELATION VS LAG (t) #### ////

            ofstream dataFile1("PNS_A4_q1_1.txt");
            int t1_max = 100; // max_sweeps/
            for(int t1=0;t1<t1_max;t1++){
                dataFile1 << t1 << " " << autoc(mag_array,true_mag_mean,t1) << "\n";
            }
            dataFile1.close();

            stringstream ss;
            ss << "gnuplot -p -e \"set xlabel 't'; "
            << "set ylabel 'Gamma(t) - Autocorrelation'; "
            << "plot 'PNS_A4_q1_1.txt' using 1:2 with lines title 'Autocorrelation of Magnetism for beta ="
            << beta
            << "';"
            //    << "'PNS_A3_q2_2.txt' using 1:2 with points linecolor rgb 'blue' pointtype 7 title 'Best Omega = "
            //    << best_omega
            << "\"";

            string command = ss.str();
            system(command.c_str());

            //===================================================================================================================================
            //// #### PLOT OF MAGNETISM VS SWEEP #### ////

            ofstream dataFile2("PNS_A4_q1_2.txt");
            for(int i=0;i<mag_array.size();i++){
                double sumdude = 0.0;
                for(int j=0;j<(i+1);j++){
                    sumdude += mag_array[j]/(i+1);
                }

                dataFile2 << i << " " << mag_array[i] << " " << sumdude << "\n";
            }
            dataFile2.close();

            stringstream ss1;
            ss1 << "gnuplot -p -e \"set xlabel 'iteration'; "
               << "set ylabel 'magnetism'; "
               << "plot 'PNS_A4_q1_2.txt' using 1:2 with points title 'Magnetism for beta ="
               << beta
               << "';"
            //    << "'PNS_A3_q2_2.txt' using 1:2 with points linecolor rgb 'blue' pointtype 7 title 'Best Omega = "
            //    << best_omega
               << "\"";

            string command1 = ss1.str();
            system(command1.c_str());

            //===================================================================================================================================
            //// #### PLOT OF AV MAG AS FUNCTION OF SWEEPS AND THE FINAL AVERAGE MAGNETISM AS A LINE #### ////

            ofstream dataFile3("PNS_A4_q1_3.txt");
            for(int i=0;i<mag_array.size();i++){
                dataFile3 << i << " " << mag_mean << "\n";
            }
            dataFile3.close();

            stringstream ss2;
            ss2 << "gnuplot -p -e \"set xlabel 'iteration'; "
            << "set ylabel 'average magnetism'; "
            << "plot 'PNS_A4_q1_2.txt' using 1:3 with lines title 'Av Magnetism for beta ="
            << beta
            << "',"
            << "'PNS_A4_q1_3.txt' using 1:2 with lines title 'Average Mag = "
            << mag_mean
            << "'\"";

            string command2 = ss2.str();
            system(command2.c_str());

            //===================================================================================================================================
            //// #### PLOT OF INTEGRATION TIME AS A FUNCTION OF W (WINDOW) #### ////

            ofstream dataFile4("PNS_A4_q1_4.txt");
            for(int w1=1;w1<(500);w1++){
                dataFile4 << w1 << " " << tau(mag_array,true_mag_mean,w1) << "\n";

            }
            dataFile4.close();

            stringstream ss4;
            ss4 << "gnuplot -p -e \"set xlabel 'w'; "
            << "set ylabel 'Tau int'; "
            << "plot 'PNS_A4_q1_4.txt' using 1:2 with lines title 'Tau int for beta ="
            << beta
            << "';"
            //    << "'PNS_A3_q2_2.txt' using 1:2 with points linecolor rgb 'blue' pointtype 7 title 'Best Omega = "
            //    << best_omega
            << "\"";

            string command4 = ss4.str();
            system(command4.c_str());


        }
        //===================================================================================================================================
        //===================================================================================================================================
        //===================================================================================================================================
        //===================================================================================================================================
        //===================================================================================================================================

        //===================================================================================================================================
        //===================================================================================================================================
        // FOR MAGNETISM SQUARED
        //===================================================================================================================================
        //===================================================================================================================================
        // if statement for plots so that I dont get a plot for every beta
        if((counter == 0.0 || counter == 25.0 || counter == nbeta/1.0)){
            //===================================================================================================================================
            //// #### PLOT OF AUTOCORRELATION VS LAG (t) #### ////

            ofstream dataFile11("PNS_A4_q1_11.txt");
            int t1_max = 100; // max_sweeps/
            for(int t1=0;t1<t1_max;t1++){
                dataFile11 << t1 << " " << autoc(mag_sqd_array,true_mag_sqd_mean,t1) << "\n";
            }
            dataFile11.close();

            stringstream ss0;
            ss0 << "gnuplot -p -e \"set xlabel 't'; "
            << "set ylabel 'Gamma(t) - Autocorrelation'; "
            << "plot 'PNS_A4_q1_11.txt' using 1:2 with lines title 'Autocorrelation of Magnetism Squared for beta ="
            << beta
            << "';"
            //    << "'PNS_A3_q2_2.txt' using 1:2 with points linecolor rgb 'blue' pointtype 7 title 'Best Omega = "
            //    << best_omega
            << "\"";

            string command0 = ss0.str();
            system(command0.c_str());

            //===================================================================================================================================
            //// #### PLOT OF MAGNETISM SQUARED VS SWEEP #### ////

            ofstream dataFile22("PNS_A4_q1_22.txt");
            for(int i=0;i<mag_sqd_array.size();i++){
                double sumdude = 0.0;
                for(int j=0;j<(i+1);j++){
                    sumdude += mag_sqd_array[j]/(i+1);
                }

                dataFile22 << i << " " << mag_sqd_array[i] << " " << sumdude << "\n";
            }
            dataFile22.close();

            stringstream ss11;
            ss11 << "gnuplot -p -e \"set xlabel 'iteration'; "
               << "set ylabel 'magnetism squared'; "
               << "plot 'PNS_A4_q1_22.txt' using 1:2 with points title 'Magnetism Squared for beta ="
               << beta
               << "';"
            //    << "'PNS_A3_q2_2.txt' using 1:2 with points linecolor rgb 'blue' pointtype 7 title 'Best Omega = "
            //    << best_omega
               << "\"";

            string command11 = ss11.str();
            system(command11.c_str());

            //===================================================================================================================================
            //// #### PLOT OF AV MAG SQD AS FUNCTION OF SWEEPS AND THE FINAL AVERAGE MAGNETISM SQD AS A LINE #### ////

            ofstream dataFile33("PNS_A4_q1_33.txt");
            for(int i=0;i<mag_sqd_array.size();i++){
                dataFile33 << i << " " << mag_sqd_mean << "\n";
            }
            dataFile33.close();

            stringstream ss22;
            ss22 << "gnuplot -p -e \"set xlabel 'iteration'; "
            << "set ylabel 'average magnetism squared'; "
            << "plot 'PNS_A4_q1_22.txt' using 1:3 with lines title 'Av Magnetism Squared for beta ="
            << beta
            << "',"
            << "'PNS_A4_q1_33.txt' using 1:2 with lines title 'Average Mag Sqd = "
            << mag_sqd_mean
            << "'\"";

            string command22 = ss22.str();
            system(command22.c_str());

            //===================================================================================================================================
            //// #### PLOT OF INTEGRATION TIME AS A FUNCTION OF W (WINDOW) #### ////

            ofstream dataFile44("PNS_A4_q1_44.txt");
            for(int w1=1;w1<(500);w1++){
                dataFile44 << w1 << " " << tau(mag_sqd_array,true_mag_sqd_mean,w1) << "\n";

            }
            dataFile44.close();

            stringstream ss44;
            ss44 << "gnuplot -p -e \"set xlabel 'w'; "
            << "set ylabel 'Tau int'; "
            << "plot 'PNS_A4_q1_44.txt' using 1:2 with lines title 'Tau int for beta ="
            << beta
            << "';"
            //    << "'PNS_A3_q2_2.txt' using 1:2 with points linecolor rgb 'blue' pointtype 7 title 'Best Omega = "
            //    << best_omega
            << "\"";

            string command44 = ss44.str();
            system(command44.c_str());


        }
        //===================================================================================================================================
        //===================================================================================================================================
        //===================================================================================================================================
        //===================================================================================================================================
        //===================================================================================================================================

        counter += 1.0;
    }
    
    //===================================================================================================================================
    //// #### PLOTS #### ////

    ofstream dataFile5("PNS_A4_q1_5.txt"); // uncertainty is just the square root of variance
            cout << "\n### FOR " << max_sweeps << " SWEEPS ###" << "\n";
            for(int i=0;i<(nbeta);i++){
                dataFile5 << beta_array[i] << " " << mag_vs_beta[i] << " " << tau_int_array[i] << " " << mag_sqd_vs_beta[i] << " " << sqrt(var_mag_vs_beta[i]) << " " << sqrt(var_mag_sqd_vs_beta[i]) << "\n";
                cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << "\n";
                cout << "Beta = " << beta_array[i] << " || " << "<M> = " << mag_vs_beta[i] << " || " << " <M^2> = " << mag_sqd_vs_beta[i] << " || " << " Uncertainty of <M> = " << sqrt(var_mag_vs_beta[i]) << " || " << " Uncertainty of <M^2> = " << sqrt(var_mag_sqd_vs_beta[i]) << "\n";
            } 
            dataFile5.close();


            //// #### PLOT OF MAGNETISM FOUND USING INTEGRATED TIME VS BETA #### ////

            stringstream ss5;
            ss5 << "gnuplot -p -e \"set xlabel 'Beta'; "
            << "set ylabel 'Magnetism'; "
            << "plot 'PNS_A4_q1_5.txt' using 1:2 with linespoints title 'Magnetism vs Beta for "
            << max_sweeps
            << " sweeps';"
            //    << "'PNS_A3_q2_2.txt' using 1:2 with points linecolor rgb 'blue' pointtype 7 title 'Best Omega = "
            //    << best_omega
            << "\"";

            string command5 = ss5.str();
            system(command5.c_str());


    //========================================================================================================================================================================================================

            //// #### PLOT OF MAGNETISM SQD FOUND USING INTEGRATED TIME VS BETA #### ////

            stringstream ss55;
            ss55 << "gnuplot -p -e \"set xlabel 'Beta'; "
            << "set ylabel 'Magnetism Squared'; "
            << "plot 'PNS_A4_q1_5.txt' using 1:4 with linespoints title 'Magnetism Squared vs Beta for "
            << max_sweeps
            << " sweeps';"
            //    << "'PNS_A3_q2_2.txt' using 1:2 with points linecolor rgb 'blue' pointtype 7 title 'Best Omega = "
            //    << best_omega
            << "\"";

            string command55 = ss55.str();
            system(command55.c_str());

            
            //===================================================================================================================================
            //// #### PLOT OF MAG UNCERTAINTY #### ////

            // ofstream dataFile6("PNS_A4_q1_6.txt");
            // for(int i=1;i<var_mag_array.size();i++){
            //     dataFile6 << beta_array[i] << " " << sqrt(var_mag_array[i]) << "\n";

            // }
            // dataFile6.close();

            stringstream ss6;
            ss6 << "gnuplot -p -e \"set xlabel 'beta'; "
            << "set ylabel 'Uncertainty'; "
            // << "set yrange [-0.1:0.1];"
            << "plot 'PNS_A4_q1_5.txt' using 1:5 with lines title 'Unertainties of Magnetism "
            // << max_sweeps
            // << " sweeps"
            << "';"
            //    << "'PNS_A3_q2_2.txt' using 1:2 with points linecolor rgb 'blue' pointtype 7 title 'Best Omega = "
            //    << best_omega
            << "\"";

            string command6 = ss6.str();
            system(command6.c_str());

            //===================================================================================================================================
            //// #### PLOT OF MAG SQD UNCERTAINTY #### ////

            // ofstream dataFile66("PNS_A4_q1_66.txt");
            // for(int i=1;i<var_mag_sqd_array.size();i++){
            //     dataFile66 << beta_array[i] << " " << sqrt(var_mag_sqd_array[i]) << "\n";

            // }
            // dataFile66.close();

            stringstream ss66;
            ss66 << "gnuplot -p -e \"set xlabel 'beta'; "
            << "set ylabel 'Uncertainty'; "
            // << "set yrange [-0.1:0.1];"
            << "plot 'PNS_A4_q1_5.txt' using 1:6 with lines title 'Unertainties of Magnetism Squared "
            // << max_sweeps
            // << " sweeps"
            << "';"
            //    << "'PNS_A3_q2_2.txt' using 1:2 with points linecolor rgb 'blue' pointtype 7 title 'Best Omega = "
            //    << best_omega
            << "\"";

            string command66 = ss66.str();
            system(command66.c_str());
    


}
