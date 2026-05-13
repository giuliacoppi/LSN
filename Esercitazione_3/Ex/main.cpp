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

   vector<double> ave_C(N, 0.0);            //Vector for average in each block
   vector<double> av2_C(N, 0.0);            //Vector for square average in each block
   vector<double> sum_prog_C(N, 0.0);       //Vector for progressive sum
   vector<double> su2_prog_C(N, 0.0);       //Vector for progressive square sum
   vector<double> err_prog_C(N, 0.0);       //Vector for error propagation

   vector<double> ave_P(N, 0.0);            //Vector for average in each block
   vector<double> av2_P(N, 0.0);            //Vector for square average in each block
   vector<double> sum_prog_P(N, 0.0);       //Vector for progressive sum
   vector<double> su2_prog_P(N, 0.0);       //Vector for progressive square sum
   vector<double> err_prog_P(N, 0.0);       //Vector for error propagation

   double S_0 = 100.0;         //Asset price at t=0
   double T = 1.0;             //Delivery time
   double K = 100.0;           //Strike price
   double r = 0.1;             //risk-free interest rate
   double sigma = 0.25;        //Volatility

   //DIRECT

   //Computation of <Ai> and <Ai>^2 for each block
   for (int i=0; i<N; i++) {
      double S=0;     
      double sum_C=0; 
      double sum_P=0;
      for (int j = 0; j < L; j++){
         S = S_0 * exp( (r - 0.5 * sigma * sigma) * T + sigma * rnd.Gauss(0, T) );
         sum_C += exp(-r*T) * max(0.0, (S - K));
         sum_P += exp(-r*T) * max(0.0, (K - S));
      }
      ave_C[i] = sum_C/L;
      av2_C[i] = pow((sum_C/L), 2);
      ave_P[i] = sum_P/L;
      av2_P[i] = pow((sum_P/L), 2);
   }

   //Computation of A and <A>^2 in function of N
   double running_sum_C=0;
   double running_sum2_C=0;
   double running_sum_P=0;
   double running_sum2_P=0;
   
   for (int i=0; i<N; i++) {
      running_sum_C += ave_C[i];
      running_sum2_C += av2_C[i];
      running_sum_P += ave_P[i];
      running_sum2_P += av2_P[i];
    
      sum_prog_C[i] = running_sum_C / (i+1);
      su2_prog_C[i] = running_sum2_C / (i+1);
      err_prog_C[i] = error(sum_prog_C[i], su2_prog_C[i], i);
      sum_prog_P[i] = running_sum_P / (i+1);
      su2_prog_P[i] = running_sum2_P / (i+1);
      err_prog_P[i] = error(sum_prog_P[i], su2_prog_P[i], i);
   }

   //Create an output file for Python
   ofstream out_direct_C("output_direct_C.dat");
   if (out_direct_C.is_open()){
      for (int i=0; i<N; i++) {
         out_direct_C << i << " " << sum_prog_C[i] << " " << err_prog_C[i] << endl;
      }
      out_direct_C.close();
   } else cerr << "PROBLEM: Unable to open output_direct_C.dat" << endl;

   ofstream out_direct_P("output_direct_P.dat");
   if (out_direct_P.is_open()){
      for (int i=0; i<N; i++) {
         out_direct_P << i << " " << sum_prog_P[i] << " " << err_prog_P[i] << endl;
      }
      out_direct_P.close();
   } else cerr << "PROBLEM: Unable to open output_direct_P.dat" << endl;


   //DISCRETIZED

   ave_C.assign(N, 0.0);            
   av2_C.assign(N, 0.0);            
   sum_prog_C.assign(N, 0.0);       
   su2_prog_C.assign(N, 0.0);       
   err_prog_C.assign(N, 0.0);       

   ave_P.assign(N, 0.0);            
   av2_P.assign(N, 0.0);            
   sum_prog_P.assign(N, 0.0);       
   su2_prog_P.assign(N, 0.0);       
   err_prog_P.assign(N, 0.0);   
   


   //Computation of <Ai> and <Ai>^2 for each block
   for (int i=0; i<N; i++) {
     
      double sum_C=0;
      double sum_P=0;

      for (int j = 0; j < L; j++){

         S_0=100.0;

         for (int k = 0; k < 100; k++){
            
            S_0 = S_0 * exp( (r - 0.5 * sigma * sigma) * 0.01 + sigma * rnd.Gauss(0, 1) * sqrt(0.01) );

         }
                  
         sum_C += exp(-r*T) * max(0.0, (S_0 - K));
         sum_P += exp(-r*T) * max(0.0, (K - S_0));
      }
      ave_C[i] = sum_C/L;
      av2_C[i] = pow((sum_C/L), 2);
      ave_P[i] = sum_P/L;
      av2_P[i] = pow((sum_P/L), 2);
   }

   //Computation of A and <A>^2 in function of N
   running_sum_C=0;
   running_sum2_C=0;
   running_sum_P=0;
   running_sum2_P=0;
   
   for (int i=0; i<N; i++) {
      running_sum_C += ave_C[i];
      running_sum2_C += av2_C[i];
      running_sum_P += ave_P[i];
      running_sum2_P += av2_P[i];
    
      sum_prog_C[i] = running_sum_C / (i+1);
      su2_prog_C[i] = running_sum2_C / (i+1);
      err_prog_C[i] = error(sum_prog_C[i], su2_prog_C[i], i);
      sum_prog_P[i] = running_sum_P / (i+1);
      su2_prog_P[i] = running_sum2_P / (i+1);
      err_prog_P[i] = error(sum_prog_P[i], su2_prog_P[i], i);
   }

   //Create an output file for Python
   ofstream out_discretized_C("out_discretized_C.dat");
   if (out_discretized_C.is_open()){
      for (int i=0; i<N; i++) {
         out_discretized_C << i << " " << sum_prog_C[i] << " " << err_prog_C[i] << endl;
      }
      out_discretized_C.close();
   } else cerr << "PROBLEM: Unable to open out_discretized_C.dat" << endl;

   ofstream out_discretized_P("out_discretized_P.dat");
   if (out_discretized_P.is_open()){
      for (int i=0; i<N; i++) {
         out_discretized_P << i << " " << sum_prog_P[i] << " " << err_prog_P[i] << endl;
      }
      out_discretized_P.close();
   } else cerr << "PROBLEM: Unable to open out_discretized_P.dat" << endl;
   

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
