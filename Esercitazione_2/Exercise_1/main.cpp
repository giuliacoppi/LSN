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

double error(double AV, double AV2, int n) { // Function for statistical uncertainty estimation
   if (n == 0) {
         return 0;
   } else {
      return sqrt((AV2 - pow(AV, 2)) / n);
   }
}

double variance(double AV, double AV2, int n) { // Function for variance estimation
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



   int M=100000;              //Total number of throws
   int N=100;                 //Number of blocks
   int L=int(M/N);            //Number of throws in each block, please use for M a multiple of N
 
   //uniform
   vector<double> ave_unif(N, 0.0);            //Vector for average in each block
   vector<double> av2_unif(N, 0.0);            //Vector for square average in each block
   vector<double> sum_prog_unif(N, 0.0);       //Vector for progressive sum
   vector<double> su2_prog_unif(N, 0.0);       //Vector for progressive square sum
   vector<double> err_prog_unif(N, 0.0);       //Vector for error propagation

   //importance sampling
   vector<double> ave_imp(N, 0.0);            //Vector for average in each block
   vector<double> av2_imp(N, 0.0);            //Vector for square average in each block
   vector<double> sum_prog_imp(N, 0.0);       //Vector for progressive sum
   vector<double> su2_prog_imp(N, 0.0);       //Vector for progressive square sum
   vector<double> err_prog_imp(N, 0.0);       //Vector for error propagation


   //Computation of <Ai> and <Ai>^2 for each block
   for (int i=0; i<N; i++) {
      double I_unif=0;
      double I_imp=0;
      for (int j = 0; j < L; j++){
         I_unif += 1 * M_PI / 2 * cos(M_PI * rnd.Rannyu() / 2) ;
         double x = rnd.Importance_Sampling();
         I_imp += 1 * M_PI / 2 * cos(M_PI * x / 2) / (2*(1-x)) ;
      }
      ave_unif[i] = I_unif/L;
      av2_unif[i] = pow((I_unif/L), 2);

      ave_imp[i] = I_imp/L;
      av2_imp[i] = pow((I_imp/L), 2);
   }

   //Computation of A and <A>^2 in function of N
   double running_sum_unif=0;
   double running_sum2_unif=0; 
   double running_sum_imp=0;
   double running_sum2_imp=0; 
   
   for (int i=0; i<N; i++) {
      running_sum_unif += ave_unif[i];
      running_sum2_unif += av2_unif[i];
      running_sum_imp += ave_imp[i];
      running_sum2_imp += av2_imp[i];
    
      sum_prog_unif[i] = running_sum_unif / (i+1);
      su2_prog_unif[i] = running_sum2_unif / (i+1);
      err_prog_unif[i] = error(sum_prog_unif[i], su2_prog_unif[i], i);
      sum_prog_imp[i] = running_sum_imp / (i+1);
      su2_prog_imp[i] = running_sum2_imp / (i+1);
      err_prog_imp[i] = error(sum_prog_imp[i], su2_prog_imp[i], i);
   }

   //Create an output file for Python
   ofstream out_unif("output_ave_unif.dat");
   if (out_unif.is_open()){
      for (int i=0; i<N; i++) {
         out_unif << i << " " << sum_prog_unif[i] - 1 << " " << err_prog_unif[i] << endl;
      }
      out_unif.close();
   } else cerr << "PROBLEM: Unable to open output_ave_unif.dat" << endl;

   ofstream out_imp("output_ave_imp.dat");
   if (out_imp.is_open()){
      for (int i=0; i<N; i++) {
         out_imp << i << " " << sum_prog_imp[i] - 1 << " " << err_prog_imp[i] << endl;
      }
      out_imp.close();
   } else cerr << "PROBLEM: Unable to open output_ave_imp.dat" << endl;

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
