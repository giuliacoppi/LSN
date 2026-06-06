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

    ////////////Simulated Annealing Simulation////////////
    double SA_steps = 200;

    double sigma = 1.0;
    double mu = 1.0;

    vector<double> H = simulate_H(rnd, mu, sigma);

    //Create an output file for Python
    ofstream out("output_8_2.dat");
    out << "#" 
        << setw(11) << "T"    
        << setw(12) << "mu" 
        << setw(12) << "sigma" 
        << setw(12) << "H"
        << setw(12) << "Error" << endl;

    //Selects the T for SA
    for (double T = 2.5; T > 0.01; T *= 0.98){
        //loops on SA_steps
        for (int i = 0; i < SA_steps; i++){
            double delta_SA = 0.5*T;        //SA delta: lets me pick mu and sigma in an appropriate interval based on the current temperature
            double mu_new = fabs(rnd.Rannyu(mu - delta_SA, mu + delta_SA));
            double sigma_new;
            do{
                sigma_new = fabs(rnd.Rannyu(sigma - delta_SA, sigma + delta_SA));
            }while(sigma_new < 0.1);

            vector<double> H_new = simulate_H(rnd, mu_new, sigma_new);

            double random_number = rnd.Rannyu();
            if(Metropolis_SA(T, H_new[0], H[0], random_number)){
                sigma = sigma_new;
                mu = mu_new;
                H = H_new;
            }
        }
        out << setw(12) << T
            << setw(12) << mu
            << setw(12) << sigma
            << setw(12) << H[0]
            << setw(12) << H[1] << endl;
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
