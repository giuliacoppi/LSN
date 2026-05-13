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

//Exercise 1.3

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

	int M=100000;              //Total number of throws
	int N=100;                 //Number of blocks
	int L=int(M/N);            //Number of throws in each block, please use for M a multiple of N

	double d=10;               //Distance between each line
	double length=8;           //Needle length

	vector<double> ave(N, 0.0);            //Vector for average in each block (pi)
	vector<double> av2(N, 0.0);            //Vector for square average in each block
	vector<double> sum_prog(N, 0.0);       //Vector for progressive sum
	vector<double> su2_prog(N, 0.0);       //Vector for progressive square sum
	vector<double> err_prog(N, 0.0);       //Vector for error propagation


	//Computation of <Ai> and <Ai>^2 for each block
	for (int i = 0; i < N; i++){
		double Nhit = 0;
		double pi = 0;

		for (int j = 0; j < L; j++){
			double CM = 0;
			double x = 0;
			double y = 0;

			CM = rnd.Rannyu(0, d/2);
			
			do{
				x = rnd.Rannyu();
				y = rnd.Rannyu();
			} while (pow(x, 2) + pow(y, 2) > 1);

			double theta = acos( x / (sqrt(pow(x, 2) + pow(y, 2))));

			if ( (length / 2) * cos(theta) > (d/2) - CM) { Nhit ++; }
		}

		if (Nhit > 0) {
			pi = (2.0 * length * L) / (Nhit * d);
		} else {
			pi = 0; 
		}
		
		ave[i] = pi;
		av2[i] = pow(pi, 2);
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
	ofstream out("output_pi.dat");
	if (out.is_open()){
		for (int i=0; i<N; i++) {
			out << i << " " << sum_prog[i] << " " << err_prog[i] << endl;
		}
		out.close();
	} else cerr << "PROBLEM: Unable to open output_pi.dat" << endl;

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
