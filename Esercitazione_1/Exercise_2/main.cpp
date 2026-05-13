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

//Exercise 1.2
 
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

    int M=10000;       //Number of iterations

    int values_N[] = {1, 2, 10, 100};

    //Create files for Python
    ofstream out_uni("out_uni.dat");
    ofstream out_exp("out_exp.dat");
    ofstream out_lor("out_lor.dat");

    for(int k = 0; k < M; k++) {

        for(int i = 0; i < 4; i++) { 
            int N = values_N[i];
            double sum_uni = 0;
            double sum_exp = 0;
            double sum_lor = 0;
            
            for(int j = 0; j < N; j++) {
                sum_uni += rnd.Rannyu();
                sum_exp += rnd.Exponential(1);
                sum_lor += rnd.Lorentz(0, 1);
            }

            //Sn
            out_uni << sum_uni/(double)N << "\t";
            out_exp << sum_exp/(double)N << "\t";
            out_lor << sum_lor/(double)N << "\t";
            
            }

        out_uni << endl;
        out_exp << endl;
        out_lor << endl;

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
