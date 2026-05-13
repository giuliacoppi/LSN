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


   int M=10000;              //Total number of throws
   int N=100;                 //Number of blocks
   int L=int(M/N);            //Number of throws in each block, please use for M a multiple of N. In this case, each throw corresponds to a RW of 100 steps S
   int S=100;                 //Number of steps per RW

   double Nx = 0;
   double Ny = 0;
   double Nz = 0;

   vector<double> ave(S, 0.0);            //Vector for average in each block, function of steps
   vector<double> av2(S, 0.0);            //Vector for square average in each block, function of steps
   vector<double> su2_prog(S, 0.0);       //Vector for progressive square sum
   vector<double> err_prog(S, 0.0);       //Vector for error propagation
   vector<double> sum_prog(S, 0.0);       //Vector for progressive sum

   vector<double> running_sum(S, 0.0);
   vector<double> running_sum2(S, 0.0);
   
   vector<double> rn2_trial (S, 0.0);     //Vector for computing rn2 for each simulation. Here I am going to accumulate the rn^2 of all the RWs in a block a as function of steps

   //CUBIC LATTICE

   //Computation of <Ai> and <Ai>^2 for each block, in this case I need the values for each step. I am coming out of the bloch with S means and S squared means
   for (int i=0; i<N; i++) { //Select Block

      rn2_trial.assign(S, 0.0);
      ave.assign(S, 0.0);
      av2.assign(S, 0.0);
         
      for (int j = 0; j < L; j++){ //Select the RW

         Nx = 0;
         Ny = 0;
         Nz = 0;

         for (int k = 1; k < S; k++){ //100 steps

            double x = rnd.Rannyu(0, 3);

            if (x < 0.5) { Nx++; }
            if (x >= 0.5 && x < 1) { Nx--; }
            if (x >= 1 && x < 1.5) { Ny++; }
            if (x >= 1.5 && x < 2) { Ny--; }
            if (x >= 2 && x < 2.5) { Nz++; }
            if (x >= 2.5) { Nz--; }

            rn2_trial[k] += ( Nx*Nx + Ny*Ny + Nz*Nz );

         }
         
      }//In rn2_trial I have the quadratic distances, as a function of the step, for all the L RWs in this block

      //The observable I need is the sqrt of <rn^2>
      for (int k = 0; k < S; k++){
         rn2_trial[k] = (rn2_trial[k] / L);
         rn2_trial[k] = sqrt(rn2_trial[k]);

         //mean and variance of each step
         running_sum[k] += rn2_trial[k];           //Accumulating the observable, as a function of steps. After simulating every block I have the sum of all the means as a function of steps
         running_sum2[k] += pow(rn2_trial[k], 2);
         ave[k] = running_sum[k] / (i+1);
         av2[k] = running_sum2[k] / (i+1);

      }
      
   }


   //Create an output file for Python
   ofstream out("output_cubic_lattice.dat");
   if (out.is_open()){
      for (int i=0; i<S; i++) {
         out << i << " " << ave[i] << " " << error(ave[i], av2[i], N) << endl;
      }
      out.close();
   } else cerr << "PROBLEM: Unable to open output_cubic_lattice.dat" << endl;


   //CONTINUUM LATTICE

   running_sum.assign(S, 0.0);
   running_sum2.assign(S, 0.0);
   sum_prog.assign(S, 0.0);
   su2_prog.assign(S, 0.0);
   err_prog.assign(S, 0.0);

   //Computation of <Ai> and <Ai>^2 for each block
   for (int i=0; i<N; i++) { //Entro nel blocco

      rn2_trial.assign(S, 0.0);
      ave.assign(S, 0.0);
      av2.assign(S, 0.0);
         
      for (int j = 0; j < L; j++){ //seleziono una simulazione 

         Nx = 0;
         Ny = 0;
         Nz = 0;

         for (int k = 1; k < S; k++){ //faccio 100 passi

            double phi = rnd.Rannyu(0, 2*M_PI);
            double cos_theta = rnd.Rannyu(-1, 1);
            double sin_theta = sqrt(1.0 - cos_theta * cos_theta);

            Nx += sin_theta * cos(phi);
            Ny += sin_theta * sin(phi);
            Nz += cos_theta;

            rn2_trial[k] += ( Nx*Nx + Ny*Ny + Nz*Nz );

         }
         
      }

      for (int k = 0; k < S; k++){
         rn2_trial[k] = (rn2_trial[k] / L);
         rn2_trial[k] = sqrt(rn2_trial[k]);

         //mean and variance of each step
         running_sum[k] += rn2_trial[k];
         running_sum2[k] += pow(rn2_trial[k], 2);
         ave[k] = running_sum[k] / (i+1);
         av2[k] = running_sum2[k] / (i+1);


      }
      
   } //esco da questo ciclo che possiedo nel running sum le medie tutte sommate


   //Create an output file for Python
   ofstream out_continuum("output_continuum_lattice.dat");
   if (out_continuum.is_open()){
      for (int i=0; i<S; i++) {
         out_continuum << i << " " << ave[i] << " " << error(ave[i], av2[i], N) << endl;
      }
      out_continuum.close();
   } else cerr << "PROBLEM: Unable to open output_continuum_lattice.dat" << endl;


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
