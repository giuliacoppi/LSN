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
#include "random.h"

using namespace std;

//Exercise 1

//Calculates the statistical uncertainty (standard deviation of the mean)
double error(double AV, double AV2, int n) { 
    if (n == 0) {
        return 0;
    }
    return sqrt((AV2 - pow(AV, 2)) / n);
}

//Returns the variance adjusted by -0.5 (simulated 'I' value). Returns 0.0 if the sample size (n) is zero.
double variance(double AV, double AV2, int n) { 
    if (n == 0) {
        return 0;
    } else {
        return (((AV2 - pow(AV, 2)) / n) - 0.5);
    }
}

 
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
    int M = 100000;              //Total number of throws
    int N = 100;                 //Number of blocks
    int L = int(M/N);            //Number of throws in each block, please use for M a multiple of N

    //Task 1

    vector<double> ave(N, 0.0);            //Vector for average in each block
    vector<double> av2(N, 0.0);            //Vector for square average in each block
    vector<double> sum_prog(N, 0.0);       //Vector for progressive sum
    vector<double> su2_prog(N, 0.0);       //Vector for progressive square sum
    vector<double> err_prog(N, 0.0);       //Vector for error propagation

    //Computation of <Ai> and <Ai>^2 for each block
    for (int i = 0; i < N; i++) {
        double sum = 0;
        for (int j = 0; j < L; j++){
            sum += rnd.Rannyu();
        }
        ave[i] = sum/L;
        av2[i] = pow((sum/L), 2);
    }

    //Computation of A and <A>^2 in function of N
    double running_sum = 0;
    double running_sum2 = 0;
    
    for (int i = 0; i < N; i++) {
        running_sum += ave[i];
        running_sum2 += av2[i];
        
        sum_prog[i] = running_sum / (i+1);
        su2_prog[i] = running_sum2 / (i+1);
        err_prog[i] = error(sum_prog[i], su2_prog[i], i);
    }

    //Create an output file for Python
    ofstream out("output_ave.dat");
    if (out.is_open()){
        for (int i = 0; i < N; i++) {
            out << i << " " << sum_prog[i] - 0.5 << " " << err_prog[i] << endl;
        }
        out.close();
    } else cerr << "PROBLEM: Unable to open output_ave.dat" << endl;

    
    //Task 2

    //Clear all vectors
    ave.assign(N, 0.0);
    av2.assign(N, 0.0);
    sum_prog.assign(N, 0.0);
    su2_prog.assign(N, 0.0);
    err_prog.assign(N, 0.0);

    //Computation of sigma_i^2 and its square mean
    for (int i = 0; i < N; i++) {
        double sum=0;
        for (int j = 0; j < L; j++){
            sum += pow((rnd.Rannyu() - 0.5), 2);
        }
        ave[i] = sum/L;
        av2[i] = pow(ave[i], 2);
    }

    //Computation of sigma^2 in funcion of N
    running_sum = 0;
    running_sum2 = 0;
    
    for (int i = 0; i < N; i++) {
        running_sum += ave[i];
        running_sum2 += av2[i];
        
        sum_prog[i] = running_sum / (i+1);
        su2_prog[i] = running_sum2 / (i+1);
        err_prog[i] = error(sum_prog[i], su2_prog[i], i);
    }

    //Create a file for Python
    ofstream out_sigma("output_sigma.dat");
    if (out_sigma.is_open()){
        for (int i=0; i<N; i++) {
            out_sigma << i << " " << sum_prog[i] - (1.0/12) << " " << err_prog[i] << endl;
        }
        out_sigma.close();
    } else cerr << "PROBLEM: Unable to open output_sigma.dat" << endl;


    //Task 3

    //Define new variables 
    int K=100;                          //Number of trials
    int I=100;                          //Number of equal intervals
    int NUM=10000;                      //Number of throws in each trial
    double chi=0;                       //Auxiliary variable

    vector<double> n_i(I, 0.0);         //Vector for ni
    vector<double> chi2(K, 0.0);        //Vector for chi2

    //Computation of chi2
    for (int k = 0; k < K; k++){
        n_i.assign(I, 0.0);

        for (int j = 0; j < NUM; j++){
            double n = rnd.Rannyu();
            int index = n * I;            //Gives back a number [0, 100)

            n_i[index] += 1;
        }

        for(int i = 0; i < I; i++){
            chi += (pow(n_i[i] - (NUM/I), 2)) / (NUM/I);
        }
        
        chi2[k] = chi;
        chi=0;
    }
    
    //Create a file for Python
    ofstream out_chi2("output_chi2.dat");
    if (out_chi2.is_open()){
        for (int i = 0; i < K; i++) {
            out_chi2 << i << " " << chi2[i] << endl;
        }
        out_chi2.close();
    } else {
        cerr << "PROBLEM: Unable to open output_chi2.dat" << endl;
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
