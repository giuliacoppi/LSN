#ifndef __FUNCTIONS_H__
#define __FUNCTIONS_H__

#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <iomanip>
#include "random.h"

using namespace std;

// Functions' definitions
double error(double AV, double AV2, int n);
double psi(double x, double sigma, double mu);
double H_psi(double x, double mu, double sigma);
double Metropolis(double old_position, double new_position, double sigma, double mu, double random_number, int & acceptance);
double set_delta(Random& rnd, double delta, double x, double mu, double sigma);
vector<double> simulate_H(Random& rnd, double mu, double sigma);



//////////////////////////////////////////////FUNCTIONS//////////////////////////////////////////////

//Calculates the statistical uncertainty (standard deviation of the mean)
double error(double AV, double AV2, int n) { 
    if (n == 0) {
        return 0;
    }
    return sqrt((AV2 - pow(AV, 2)) / n);
}


//Wave function
double psi(double x, double mu, double sigma){
    return (exp(-((x-mu)*(x-mu))/(2*sigma*sigma)) + exp(-((x+mu)*(x+mu))/(2*sigma*sigma)));
}

//Computes the total local energy through the formula in the notebook
double H_psi(double x, double mu, double sigma){
    return -0.5 *(  -1/(sigma*sigma) *  exp(-(x-mu)*(x-mu)/(2*sigma*sigma)) -1/(sigma*sigma) *  exp(-(x+mu)*(x+mu)/(2*sigma*sigma)) + 1/(pow(sigma,4)) * (x-mu)*(x-mu) * exp(-(x-mu)*(x-mu)/(2*sigma*sigma)) + 1/(pow(sigma,4)) * (x+mu)*(x+mu) * exp(-(x+mu)*(x+mu)/(2*sigma*sigma)) )/psi(x, mu, sigma) + pow(x,4) - 5./2. * pow(x,2);
}

//Metropolis algorithm
double Metropolis(double old_position, double new_position, double sigma, double mu, double random_number, int & acceptance){

    //Probabilities
    double P_old = pow(psi(old_position, mu, sigma), 2);
    double P_new = pow(psi(new_position, mu, sigma), 2);
    double P = min(1.0, P_new/P_old);
    if ( random_number < P){
        acceptance++;
        return new_position;
    }
    return old_position;
}

//Function for tuning delta to obtain an acceptance rate around 50%
double set_delta(Random& rnd, double delta, double x, double mu, double sigma){

    bool calibrated = false;
    int max_iteration = 500;        //max number of iterations to find the best fitting delta
    int iteration = 0;
    int steps = 1000;               //number of equilibration steps
    double expected = 50.0;         //acceptance 50%
    double offset = 5.0;            //I am allowing acceptance to be in [45%-55%]

    while(!calibrated && iteration < max_iteration){
        int acceptance = 0;

        for (int j = 0; j < steps; j++){
            double random_number = rnd.Rannyu();
        
            double new_x = rnd.Rannyu(x - delta, x + delta);
            x = Metropolis(x, new_x, sigma, mu, random_number, acceptance);
        }
        double rate = acceptance/(double)steps*100.0;

        if (rate >= (expected - offset) && rate <= (expected + offset)){
            calibrated = true;
        }else{
            if(rate > 0.0){
                delta = delta * (rate / 50.0);
            }else{
                delta *= 0.5; //protection if the step was too big and I have acceptance too low
            }
        }
        iteration++;  
    }
    return delta;  
}

//This fuction generates one step of our SA simulation
vector<double> simulate_H (Random& rnd, double mu, double sigma){
    //Simulation parameters
    double x = 0.0;       //starting position
    double delta_start = 1.0;   //initial delta to tune

    double delta = set_delta(rnd, delta_start, x, mu, sigma);     //delta is now tuned correctly for the simulation

    int M = 20000;              //Total number of throws
    int N = 100;                 //Number of blocks
    int L = int(M/N);           //Number of throws in each block, please use for M a multiple of N

    vector<double> H(2, 0.0);   //<H> and its uncertainty

    //Accumulators for data blocking
    double running_sum = 0;     
    double running_sum2 = 0;

    //With the starting parameters, it proposes a new postion to Metropolis. In the end it gives back <H> and its uncertainty for one SA step
    for (int i = 0; i < N; i++){
        double energy = 0;
        int acc = 0;
        for (int j = 0; j < L; j++){
            double random_number = rnd.Rannyu();
            double x_new = rnd.Rannyu(x - delta, x + delta);
            x = Metropolis(x, x_new, sigma, mu, random_number, acc);
            energy += H_psi(x, mu, sigma);
        }
        running_sum +=  energy/double(L);
        running_sum2 += pow(energy/double(L), 2);
    }
    //I just need the final value of <H> from data blocking
    H[0] = running_sum/double(N);
    H[1] = error(running_sum/double(N), running_sum2/double(N), N-1);

    return H;
}

#endif

