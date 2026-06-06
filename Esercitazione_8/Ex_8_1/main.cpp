/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include <iomanip>
#include "random.h"
#include "functions.h"

using namespace std;


int main (int argc, char *argv[]){

    Random rnd;
    int seed[4];
    int p1, p2;
    ifstream Primes("Primes");
    if (Primes.is_open()){
        Primes >> p1 >> p2 ;
    } else cerr << "PROBLEM: Unable to open Primes" << endl;
    Primes.close();

    ifstream input("seed.in");
    string property;
    if (input.is_open()){
        while ( !input.eof() ){
            input >> property;
            if( property == "RANDOMSEED" ){
                input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
                rnd.SetRandom(seed,p1,p2);
            }
        }
        input.close();
    } else cerr << "PROBLEM: Unable to open seed.in" << endl;

    //Initializing parameters
    int M = 1000000;             //Total number of throws
    int N = 100;                 //Number of blocks
    int L = int(M/N);            //Number of throws in each block, please use for M a multiple of N

    vector<double> ave(N, 0.0);            //Vector for average in each block
    vector<double> av2(N, 0.0);            //Vector for square average in each block
    vector<double> sum_prog(N, 0.0);       //Vector for progressive sum
    vector<double> su2_prog(N, 0.0);       //Vector for progressive square sum
    vector<double> err_prog(N, 0.0);       //Vector for error propagation

    vector<double> acc_vec(N, 0.0);                         //Vector for acceptance, I want to check that my set_delta works
    int n_bins = 100;                                       //Number of bins used to fill the histogram
    int index_hist = n_bins*N;                              //n_bins for each block
    vector<double> positions_histogram(index_hist, 0.0);    //Vector for the histograms

    vector<double> sum_prog_hist(index_hist, 0.0);       //Vector for progressive sum for the histogram
    vector<double> su2_prog_hist(index_hist, 0.0);       //Vector for progressive square sum for the histogram
    vector<double> err_prog_hist(index_hist, 0.0);       //Vector for error propagation for the histogram

    //Parameters obtained from 8.2
    vector<double> mu_test = {0.2, 0.8, 0.797172};
    vector<double> sigma_test = {1.5, 0.6, 0.569287};

    //Test for two couples (mu, sigma) + definitive couple
    for (int i = 0; i < 3; i++){
        //Starting parameters
        double mu = mu_test[i];
        double sigma = sigma_test[i];
        double initial_position = 0.0;
        double delta = set_delta(rnd, 1.0, initial_position, mu, sigma); //initial delta=1.0
        double bin_width = 10.0/double(n_bins);                          //I know that the wave function will stay inside [-5, 5]: x_min=-5, x_max=+5
        int acceptance = 0;
        positions_histogram.assign(index_hist, 0.0);

        //Equilibration
        for(int j=0; j<2000; j++){ 
            double random_number = rnd.Rannyu();
            double new_position = rnd.Rannyu(initial_position - delta, initial_position + delta);
            initial_position = Metropolis(initial_position, new_position, sigma, mu, random_number, acceptance);
        }

        for (int j=0; j<N; j++) {
            double sum_energy = 0;
            acceptance = 0;
            for (int k = 0; k < L; k++){
                double random_number = rnd.Rannyu();
                double new_position = rnd.Rannyu(initial_position - delta, initial_position + delta);
                initial_position = Metropolis(initial_position, new_position, sigma, mu, random_number, acceptance);
                sum_energy += H_psi(initial_position, mu, sigma); 

                //Build the PDF histogram
                int index_bin = int((initial_position+5)/bin_width);        //Protection for negative index
                if(index_bin >= 0 && index_bin < n_bins){
                    positions_histogram[j*n_bins + index_bin] += 1.0/double(L*bin_width);
                }
            }

            acc_vec[j] = double(acceptance)/double(L);
            ave[j] = sum_energy/L;
            av2[j] = pow((sum_energy/L), 2);
        }

        //Computation of A and <A>^2 in function of N
        double running_sum = 0;
        double running_sum2= 0;
    
        for (int j = 0; j < N; j++) {
            running_sum += ave[j];
            running_sum2 += av2[j];
        
            sum_prog[j] = running_sum / (j+1);
            su2_prog[j] = running_sum2 / (j+1);
            err_prog[j] = error(sum_prog[j], su2_prog[j], j);
        }

        //Computation of A and <A>^2 in function of N for the histogram
        for (int j = 0; j < n_bins; j++){
        
            double running_sum_hist = 0;
            double running_sum2_hist= 0;
    
            for (int k = 0; k < N; k++) {
                running_sum_hist += positions_histogram[k*n_bins + j];
                running_sum2_hist += pow(positions_histogram[k*n_bins + j], 2);
        
                sum_prog_hist[k*n_bins + j] = running_sum_hist / (k+1);
                su2_prog_hist[k*n_bins + j] = running_sum2_hist / (k+1);
                err_prog_hist[k*n_bins + j] = error(sum_prog_hist[k*n_bins + j], su2_prog_hist[k*n_bins + j], k);
            }
        }

        //Create an output file for Python
        string energy_file_name = "energy_test_" + to_string(i + 1) + ".dat";
        ofstream out_energy(energy_file_name);
        if (out_energy.is_open()){
            out_energy << "#" 
                     << setw(11) << "Block"    
                     << setw(12) << "Prog sum" 
                     << setw(12) << "Prog error" << endl;
            for (int j = 0; j < N; j++) {
                out_energy << setw(12) << j << setw(12) << sum_prog[j] << setw(12) << err_prog[j] << endl;
            }
            out_energy.close();
        } else cerr << "PROBLEM: Unable to open "  << energy_file_name << endl;

        string acc_file_name = "acceptance_test_" + to_string(i + 1) + ".dat";
        ofstream out_acc(acc_file_name);
            out_acc << "#" 
                    << setw(11) << "Block"    
                    << setw(12) << "Acceptance" << endl;
        if (out_acc.is_open()){
            for (int j = 0; j < N; j++) {
                out_acc << setw(12) << j << setw(12) << acc_vec[j] << endl;
            }
            out_acc.close();
        } else cerr << "PROBLEM: Unable to open "  << acc_file_name << endl;

        string hist_file_name = "hist_test_" + to_string(i + 1) + ".dat";
        ofstream out_hist(hist_file_name);
        if (out_hist.is_open()){
            out_hist << "#" 
                     << setw(11) << "Block = 100"    
                     << setw(12) << "x" 
                     << setw(12) << "Prog sum" 
                     << setw(12) << "Prog error" << endl;
            for (int j = N-1; j < N; j++) {
                for (int bin_idx = 0; bin_idx < n_bins; bin_idx++) {
                    int prog_idx = j*n_bins + bin_idx;
                    double x_coordinate = -5.0 + (bin_idx + 0.5)*bin_width;
                    out_hist << setw(12) << j + 1 << setw(12) << x_coordinate << setw(12) << sum_prog_hist[prog_idx] << setw(12) << err_prog_hist[prog_idx] << endl;
                }
            }
            out_hist.close();
        } else cerr << "PROBLEM: Unable to open " << hist_file_name << endl;


    }
    

    return 0;
}




/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
