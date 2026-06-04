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

//////////////////////////////////////////////FUNCTIONS//////////////////////////////////////////////

//Calculates the statistical uncertainty (standard deviation of the mean)
double error(double AV, double AV2, int n) { 
    if (n == 0) {
        return 0;
    }
    return sqrt((AV2 - pow(AV, 2)) / n);
}

//Bohr radius
double a = 0.0529;

//Definitions of the wave functions
double psi_100(vector<double> position){
    double r = sqrt(position[0]*position[0] + position[1]*position[1] + position[2]*position[2]);
    return (1/(M_PI * a*a*a) * exp(-2*r/a)); //definition of probability 
}
double psi_210(vector<double> position){
    double r = sqrt(position[0]*position[0] + position[1]*position[1] + position[2]*position[2]);
    return ((position[2]*position[2]) / (64*M_PI * pow(a, 5)) * exp (-r/a)); //definition of probability
}
//Metropolis algorithm
vector<double> Metropolis_100(vector<double> old_position, vector<double> new_position, double random_number, int & acceptance){
    //Probabilities
    double P_old = psi_100(old_position);
    double P_new = psi_100(new_position);
    if (P_new > P_old || random_number < (P_new/P_old)){
        acceptance++;
        return new_position;
    }
    return old_position;
}
vector<double> Metropolis_210(vector<double> old_position, vector<double> new_position, double random_number, int & acceptance){
    //Probabilities
    double P_old = psi_210(old_position);
    double P_new = psi_210(new_position);
    if (P_new > P_old || random_number < (P_new/P_old)){
        acceptance++;
        return new_position;
    }
    return old_position;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////

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


    ////////////////////////////////////////////////UNIFORM////////////////////////////////////////////////

    //Initializing position files for 3D plots
    ofstream pos_near_unif_100("pos_near_unif_100.dat");
    if (pos_near_unif_100.is_open()){
        pos_near_unif_100 << "x\ty\tz" << endl;
        pos_near_unif_100.close();
    } else cerr << "PROBLEM: Unable to open pos_near_unif_100.dat" << endl;

    ofstream pos_far_unif_100("pos_far_unif_100.dat");
    if (pos_far_unif_100.is_open()){
        pos_far_unif_100 << "x\ty\tz" << endl;
        pos_far_unif_100.close();
    } else cerr << "PROBLEM: Unable to open pos_far_unif_100.dat" << endl;

    ofstream pos_near_unif_210("pos_near_unif_210.dat");
    if (pos_near_unif_210.is_open()){
        pos_near_unif_210 << "x\ty\tz" << endl;
        pos_near_unif_210.close();
    } else cerr << "PROBLEM: Unable to open pos_near_unif_210.dat" << endl;

    ofstream pos_far_unif_210("pos_far_unif_210.dat");
    if (pos_far_unif_210.is_open()){
        pos_far_unif_210 << "x\ty\tz" << endl;
        pos_far_unif_210.close();
    } else cerr << "PROBLEM: Unable to open pos_far_unif_210.dat" << endl;

    ////////////////////////////////NEAR THE ORIGIN////////////////////////////////

    //100
    vector<double> ave_100(N, 0.0);            //Vector for average in each block
    vector<double> av2_100(N, 0.0);            //Vector for square average in each block
    vector<double> sum_prog_100(N, 0.0);       //Vector for progressive sum
    vector<double> su2_prog_100(N, 0.0);       //Vector for progressive square sum
    vector<double> err_prog_100(N, 0.0);       //Vector for error propagation

    //210
    vector<double> ave_210(N, 0.0);            //Vector for average in each block
    vector<double> av2_210(N, 0.0);            //Vector for square average in each block
    vector<double> sum_prog_210(N, 0.0);       //Vector for progressive sum
    vector<double> su2_prog_210(N, 0.0);       //Vector for progressive square sum
    vector<double> err_prog_210(N, 0.0);       //Vector for error propagation

    vector<double> initial_position_100(3, a); //STARTING NEAR THE ORIGIN, P=0 never for 100
    vector<double> initial_position_210(3, a); //STARTING NEAR THE ORIGIN, P=0 in r=0 and xy plane (z=0) for 210

    //Equilibration
    for(int i=0; i<2000; i++){ 
        int acceptance = 0;
        double random_number = rnd.Rannyu();
        //100
        vector<double> new_position_100 = {rnd.Rannyu(initial_position_100[0] - a, initial_position_100[0] + a), rnd.Rannyu(initial_position_100[1] - a, initial_position_100[1] + a), rnd.Rannyu(initial_position_100[2] - a, initial_position_100[2] + a)};
        initial_position_100 = Metropolis_100(initial_position_100, new_position_100, random_number, acceptance);
        //210
        vector<double> new_position_210 = {rnd.Rannyu(initial_position_210[0] - 4.8*a, initial_position_210[0] + 4.8*a), rnd.Rannyu(initial_position_210[1] - 4.8*a, initial_position_210[1] + 4.8*a), rnd.Rannyu(initial_position_210[2] - 4.8*a, initial_position_210[2] + 4.8*a)};
        initial_position_210 = Metropolis_210(initial_position_210, new_position_210, random_number, acceptance);
    }

    //Open files (append) to write the coordinates
    ofstream write_pos_near_unif_100("pos_near_unif_100.dat", ios::app);
    ofstream write_pos_near_unif_210("pos_near_unif_210.dat", ios::app);

    if (write_pos_near_unif_100.is_open() && write_pos_near_unif_210.is_open()){

        for (int i = 0; i < N; i++){
            double final_position_100 = 0;
            int acceptance_100 = 0;
            double final_position_210 = 0;
            int acceptance_210 = 0;
            for (int j = 0; j  < L; j++){
                double random_number = rnd.Rannyu();

                //100
                vector<double> new_position_100 = {rnd.Rannyu(initial_position_100[0] - a*1.2, initial_position_100[0] + a*1.2), rnd.Rannyu(initial_position_100[1] - a*1.2, initial_position_100[1] + a*1.2), rnd.Rannyu(initial_position_100[2] - a*1.2, initial_position_100[2] + a*1.2)};
                write_pos_near_unif_100 << initial_position_100[0] << "\t" 
                                        << initial_position_100[1] << "\t" 
                                        << initial_position_100[2] << endl;
                initial_position_100 = Metropolis_100(initial_position_100, new_position_100, random_number, acceptance_100);

                //210
                vector<double> new_position_210 = {rnd.Rannyu(initial_position_210[0] - 4.8*a, initial_position_210[0] + 4.8*a), rnd.Rannyu(initial_position_210[1] - 4.8*a, initial_position_210[1] + 4.8*a), rnd.Rannyu(initial_position_210[2] - 4.8*a, initial_position_210[2] + 4.8*a)};
                write_pos_near_unif_210 << initial_position_210[0] << "\t" 
                                        << initial_position_210[1] << "\t" 
                                        << initial_position_210[2] << endl;
                initial_position_210 = Metropolis_210(initial_position_210, new_position_210, random_number, acceptance_210);
                
                final_position_100 += sqrt(initial_position_100[0]*initial_position_100[0] + initial_position_100[1]*initial_position_100[1] + initial_position_100[2]*initial_position_100[2]);
                final_position_210 += sqrt(initial_position_210[0]*initial_position_210[0] + initial_position_210[1]*initial_position_210[1] + initial_position_210[2]*initial_position_210[2]);
            }
            //cout << i << "      " << acceptance_100/double(L)*100 << endl;
            //cout << i << "      " << acceptance_210/double(L)*100 << endl;

            ave_100[i] = final_position_100/L;
            av2_100[i] = pow((final_position_100/L), 2); 
            ave_210[i] = final_position_210/L;
            av2_210[i] = pow((final_position_210/L), 2);    
        }
        write_pos_near_unif_100.close();
        write_pos_near_unif_210.close();
    } else {
    cerr << "PROBLEM: Unable to open file for appending!" << endl;
    }

    //Computation of A and <A>^2 in function of N
    double running_sum_100 = 0;
    double running_sum2_100 = 0;
    double running_sum_210 = 0;
    double running_sum2_210 = 0;
    
    for (int i = 0; i < N; i++) {
        running_sum_100 += ave_100[i];
        running_sum2_100 += av2_100[i];
        running_sum_210 += ave_210[i];
        running_sum2_210 += av2_210[i];
        
        sum_prog_100[i] = running_sum_100 / (i+1);
        su2_prog_100[i] = running_sum2_100 / (i+1);
        err_prog_100[i] = error(sum_prog_100[i], su2_prog_100[i], i);
        sum_prog_210[i] = running_sum_210 / (i+1);
        su2_prog_210[i] = running_sum2_210 / (i+1);
        err_prog_210[i] = error(sum_prog_210[i], su2_prog_210[i], i);
    }

    //Create an output file for Python
    ofstream out_100("output_100.dat");
    if (out_100.is_open()){
        for (int i = 0; i < N; i++) {
            out_100 << i << " " << sum_prog_100[i] << " " << err_prog_100[i] << endl;
        }
        out_100.close();
    } else cerr << "PROBLEM: Unable to open output_100.dat" << endl;
    ofstream out_210("output_210.dat");
    if (out_210.is_open()){
        for (int i = 0; i < N; i++) {
            out_210 << i << " " << sum_prog_210[i] << " " << err_prog_210[i] << endl;
        }
        out_210.close();
    } else cerr << "PROBLEM: Unable to open output_210.dat" << endl;





    //STARTING FAR AWAY

    //100
    ave_100.assign(N, 0.0);            //Vector for average in each block
    av2_100.assign(N, 0.0);            //Vector for square average in each block
    sum_prog_100.assign(N, 0.0);       //Vector for progressive sum
    su2_prog_100.assign(N, 0.0);       //Vector for progressive square sum
    err_prog_100.assign(N, 0.0);       //Vector for error propagation

    //210
    ave_210.assign(N, 0.0);            //Vector for average in each block
    av2_210.assign(N, 0.0);            //Vector for square average in each block
    sum_prog_210.assign(N, 0.0);       //Vector for progressive sum
    su2_prog_210.assign(N, 0.0);       //Vector for progressive square sum
    err_prog_210.assign(N, 0.0);       //Vector for error propagation

    initial_position_100.assign(3, 100*a); //STARTING FAR AWAY FROM THE ORIGIN, P=0 never for 100
    initial_position_210.assign(3, 100*a); //STARTING FAR AWAY FROM THE ORIGIN, P=0 in r=0 and xy plane (z=0) for 210

    //Equilibration
    for(int i=0; i<2000; i++){ 
        int acceptance = 0;
        double random_number = rnd.Rannyu();
        //100
        vector<double> new_position_100 = {rnd.Rannyu(initial_position_100[0] - a, initial_position_100[0] + a), rnd.Rannyu(initial_position_100[1] - a, initial_position_100[1] + a), rnd.Rannyu(initial_position_100[2] - a, initial_position_100[2] + a)};
        initial_position_100 = Metropolis_100(initial_position_100, new_position_100, random_number, acceptance);
        //210
        vector<double> new_position_210 = {rnd.Rannyu(initial_position_210[0] - 4.8*a, initial_position_210[0] + 4.8*a), rnd.Rannyu(initial_position_210[1] - 4.8*a, initial_position_210[1] + 4.8*a), rnd.Rannyu(initial_position_210[2] - 4.8*a, initial_position_210[2] + 4.8*a)};
        initial_position_210 = Metropolis_210(initial_position_210, new_position_210, random_number, acceptance);
    }
    
    //Open files (append) to write the coordinates
    ofstream write_pos_far_unif_100("pos_far_unif_100.dat", ios::app);
    ofstream write_pos_far_unif_210("pos_far_unif_210.dat", ios::app);

    if (write_pos_far_unif_100.is_open() && write_pos_far_unif_210.is_open()){
        for (int i = 0; i < N; i++){
            double final_position_100 = 0;
            int acceptance_100 = 0;
            double final_position_210 = 0;
            int acceptance_210 = 0;
            for (int j = 0; j  < L; j++){
                double random_number = rnd.Rannyu();

                //100
                vector<double> new_position_100 = {rnd.Rannyu(initial_position_100[0] - a*1.2, initial_position_100[0] + a*1.2), rnd.Rannyu(initial_position_100[1] - a*1.2, initial_position_100[1] + a*1.2), rnd.Rannyu(initial_position_100[2] - a*1.2, initial_position_100[2] + a*1.2)};
                write_pos_far_unif_100 << initial_position_100[0] << "\t" 
                                       << initial_position_100[1] << "\t" 
                                       << initial_position_100[2] << endl;
                initial_position_100 = Metropolis_100(initial_position_100, new_position_100, random_number, acceptance_100);
                
                //210
                vector<double> new_position_210 = {rnd.Rannyu(initial_position_210[0] - 4.8*a, initial_position_210[0] + 4.8*a), rnd.Rannyu(initial_position_210[1] - 4.8*a, initial_position_210[1] + 4.8*a), rnd.Rannyu(initial_position_210[2] - 4.8*a, initial_position_210[2] + 4.8*a)};
                write_pos_far_unif_210 << initial_position_210[0] << "\t" 
                                       << initial_position_210[1] << "\t" 
                                       << initial_position_210[2] << endl;
                initial_position_210 = Metropolis_210(initial_position_210, new_position_210, random_number, acceptance_210);
                
                final_position_100 += sqrt(initial_position_100[0]*initial_position_100[0] + initial_position_100[1]*initial_position_100[1] + initial_position_100[2]*initial_position_100[2]);
                final_position_210 += sqrt(initial_position_210[0]*initial_position_210[0] + initial_position_210[1]*initial_position_210[1] + initial_position_210[2]*initial_position_210[2]);
            }
            //cout << i << "      " << acceptance_100/double(L)*100 << endl;            disable the commeting if you want to see the acceptance rate (circa 50% for both)
            //cout << i << "      " << acceptance_210/double(L)*100 << endl;

            ave_100[i] = final_position_100/L;
            av2_100[i] = pow((final_position_100/L), 2); 
            ave_210[i] = final_position_210/L;
            av2_210[i] = pow((final_position_210/L), 2);    
        } 
        write_pos_far_unif_100.close();
        write_pos_far_unif_210.close();
    } else {
    cerr << "PROBLEM: Unable to open file for appending!" << endl;
    }

    //Computation of A and <A>^2 in function of N
    running_sum_100 = 0;
    running_sum2_100 = 0;
    running_sum_210 = 0;
    running_sum2_210 = 0;
    
    for (int i = 0; i < N; i++) {
        running_sum_100 += ave_100[i];
        running_sum2_100 += av2_100[i];
        running_sum_210 += ave_210[i];
        running_sum2_210 += av2_210[i];
        
        sum_prog_100[i] = running_sum_100 / (i+1);
        su2_prog_100[i] = running_sum2_100 / (i+1);
        err_prog_100[i] = error(sum_prog_100[i], su2_prog_100[i], i);
        sum_prog_210[i] = running_sum_210 / (i+1);
        su2_prog_210[i] = running_sum2_210 / (i+1);
        err_prog_210[i] = error(sum_prog_210[i], su2_prog_210[i], i);
    }

    //Create an output file for Python
    ofstream out_100_far("output_100_far.dat");
    if (out_100_far.is_open()){
        for (int i = 0; i < N; i++) {
            out_100_far << i << " " << sum_prog_100[i] << " " << err_prog_100[i] << endl;
        }
        out_100_far.close();
    } else cerr << "PROBLEM: Unable to open output_100_far.dat" << endl;
    ofstream out_210_far("output_210_far.dat");
    if (out_210_far.is_open()){
        for (int i = 0; i < N; i++) {
            out_210_far << i << " " << sum_prog_210[i] << " " << err_prog_210[i] << endl;
        }
        out_210_far.close();
    } else cerr << "PROBLEM: Unable to open output_210_far.dat" << endl;







    ////////////////////////////////////////////////GAUSS////////////////////////////////////////////////

    //Initializing position files for 3D plots
    ofstream pos_near_gauss_100("pos_near_gauss_100.dat");
    if (pos_near_gauss_100.is_open()){
        pos_near_gauss_100 << "x\ty\tz" << endl;
        pos_near_gauss_100.close();
    } else cerr << "PROBLEM: Unable to open pos_near_gauss_100.dat" << endl;

    ofstream pos_far_gauss_100("pos_far_gauss_100.dat");
    if (pos_far_gauss_100.is_open()){
        pos_far_gauss_100 << "x\ty\tz" << endl;
        pos_far_gauss_100.close();
    } else cerr << "PROBLEM: Unable to open pos_far_gauss_100.dat" << endl;

    ofstream pos_near_gauss_210("pos_near_gauss_210.dat");
    if (pos_near_gauss_210.is_open()){
        pos_near_gauss_210 << "x\ty\tz" << endl;
        pos_near_gauss_210.close();
    } else cerr << "PROBLEM: Unable to open pos_near_gauss_210.dat" << endl;

    ofstream pos_far_gauss_210("pos_far_gauss_210.dat");
    if (pos_far_gauss_210.is_open()){
        pos_far_gauss_210 << "x\ty\tz" << endl;
        pos_far_gauss_210.close();
    } else cerr << "PROBLEM: Unable to open pos_far_gauss_210.dat" << endl;

    ////////////////////////////////NEAR THE ORIGIN////////////////////////////////

    //100
    ave_100.assign(N, 0.0);            //Vector for average in each block
    av2_100.assign(N, 0.0);            //Vector for square average in each block
    sum_prog_100.assign(N, 0.0);       //Vector for progressive sum
    su2_prog_100.assign(N, 0.0);       //Vector for progressive square sum
    err_prog_100.assign(N, 0.0);       //Vector for error propagation

    //210
    ave_210.assign(N, 0.0);            //Vector for average in each block
    av2_210.assign(N, 0.0);            //Vector for square average in each block
    sum_prog_210.assign(N, 0.0);       //Vector for progressive sum
    su2_prog_210.assign(N, 0.0);       //Vector for progressive square sum
    err_prog_210.assign(N, 0.0);       //Vector for error propagation

    initial_position_100.assign(3, a); //STARTING FAR NEAR THE ORIGIN, P=0 never for 100
    initial_position_210.assign(3, a); //STARTING FAR NEAR THE ORIGIN, P=0 in r=0 and xy plane (z=0) for 210

    //Equilibration
    for(int i=0; i<2000; i++){ 
        int acceptance = 0;
        double random_number = rnd.Rannyu();
        //100
        vector<double> new_position_100 = {rnd.Gauss(initial_position_100[0], 0.75*a), rnd.Gauss(initial_position_100[1], 0.75*a), rnd.Gauss(initial_position_100[2], 0.75*a)};
        initial_position_100 = Metropolis_100(initial_position_100, new_position_100, random_number, acceptance);
        //210
        vector<double> new_position_210 = {rnd.Gauss(initial_position_210[0], 1.87*a), rnd.Gauss(initial_position_210[1], 1.87*a), rnd.Gauss(initial_position_210[2], 1.87*a)};
        initial_position_210 = Metropolis_210(initial_position_210, new_position_210, random_number, acceptance);
    }

    //Open files (append) to write the coordinates
    ofstream write_pos_near_gauss_100("pos_near_gauss_100.dat", ios::app);
    ofstream write_pos_near_gauss_210("pos_near_gauss_210.dat", ios::app);

    if (write_pos_near_gauss_100.is_open() && write_pos_near_gauss_210.is_open()){

        for (int i = 0; i < N; i++){
            double final_position_100 = 0;
            int acceptance_100 = 0;
            double final_position_210 = 0;
            int acceptance_210 = 0;
            for (int j = 0; j  < L; j++){
                double random_number = rnd.Rannyu();

                //100
                vector<double> new_position_100 = {rnd.Gauss(initial_position_100[0], 0.75*a), rnd.Gauss(initial_position_100[1], 0.75*a), rnd.Gauss(initial_position_100[2], 0.75*a)};
                write_pos_near_gauss_100 << initial_position_100[0] << "\t" 
                                         << initial_position_100[1] << "\t" 
                                         << initial_position_100[2] << endl;
                initial_position_100 = Metropolis_100(initial_position_100, new_position_100, random_number, acceptance_100);

                //210
                vector<double> new_position_210 = {rnd.Gauss(initial_position_210[0], 1.87*a), rnd.Gauss(initial_position_210[1], 1.87*a), rnd.Gauss(initial_position_210[2], 1.87*a)};
                write_pos_near_gauss_210 << initial_position_210[0] << "\t" 
                                         << initial_position_210[1] << "\t" 
                                         << initial_position_210[2] << endl;
                initial_position_210 = Metropolis_210(initial_position_210, new_position_210, random_number, acceptance_210);
                
                final_position_100 += sqrt(initial_position_100[0]*initial_position_100[0] + initial_position_100[1]*initial_position_100[1] + initial_position_100[2]*initial_position_100[2]);
                final_position_210 += sqrt(initial_position_210[0]*initial_position_210[0] + initial_position_210[1]*initial_position_210[1] + initial_position_210[2]*initial_position_210[2]);
            }
            //cout << i << "      " << acceptance_100/double(L)*100 << endl;
            //cout << i << "      " << acceptance_210/double(L)*100 << endl;

            ave_100[i] = final_position_100/L;
            av2_100[i] = pow((final_position_100/L), 2); 
            ave_210[i] = final_position_210/L;
            av2_210[i] = pow((final_position_210/L), 2);    
        }
        write_pos_near_gauss_100.close();
        write_pos_near_gauss_210.close();
    } else {
    cerr << "PROBLEM: Unable to open file for appending!" << endl;
    }

    //Computation of A and <A>^2 in function of N
    running_sum_100 = 0;
    running_sum2_100 = 0;
    running_sum_210 = 0;
    running_sum2_210 = 0;
    
    for (int i = 0; i < N; i++) {
        running_sum_100 += ave_100[i];
        running_sum2_100 += av2_100[i];
        running_sum_210 += ave_210[i];
        running_sum2_210 += av2_210[i];
        
        sum_prog_100[i] = running_sum_100 / (i+1);
        su2_prog_100[i] = running_sum2_100 / (i+1);
        err_prog_100[i] = error(sum_prog_100[i], su2_prog_100[i], i);
        sum_prog_210[i] = running_sum_210 / (i+1);
        su2_prog_210[i] = running_sum2_210 / (i+1);
        err_prog_210[i] = error(sum_prog_210[i], su2_prog_210[i], i);
    }

    //Create an output file for Python
    ofstream out_100_gauss("output_100_gauss.dat");
    if (out_100_gauss.is_open()){
        for (int i = 0; i < N; i++) {
            out_100_gauss << i << " " << sum_prog_100[i] << " " << err_prog_100[i] << endl;
        }
        out_100_gauss.close();
    } else cerr << "PROBLEM: Unable to open output_100_gauss.dat" << endl;
    ofstream out_210_gauss("output_210_gauss.dat");
    if (out_210_gauss.is_open()){
        for (int i = 0; i < N; i++) {
            out_210_gauss << i << " " << sum_prog_210[i] << " " << err_prog_210[i] << endl;
        }
        out_210_gauss.close();
    } else cerr << "PROBLEM: Unable to open output_210_gauss.dat" << endl;





    //STARTING FAR AWAY

    //100
    ave_100.assign(N, 0.0);            //Vector for average in each block
    av2_100.assign(N, 0.0);            //Vector for square average in each block
    sum_prog_100.assign(N, 0.0);       //Vector for progressive sum
    su2_prog_100.assign(N, 0.0);       //Vector for progressive square sum
    err_prog_100.assign(N, 0.0);       //Vector for error propagation

    //210
    ave_210.assign(N, 0.0);            //Vector for average in each block
    av2_210.assign(N, 0.0);            //Vector for square average in each block
    sum_prog_210.assign(N, 0.0);       //Vector for progressive sum
    su2_prog_210.assign(N, 0.0);       //Vector for progressive square sum
    err_prog_210.assign(N, 0.0);       //Vector for error propagation

    initial_position_100.assign(3, 100*a); //STARTING FAR AWAY FROM THE ORIGIN, P=0 never for 100
    initial_position_210.assign(3, 100*a); //STARTING FAR AWAY FROM THE ORIGIN, P=0 in r=0 and xy plane (z=0) for 210

    //Equilibration
    for(int i=0; i<2000; i++){ 
        int acceptance = 0;
        double random_number = rnd.Rannyu();
        //100
        vector<double> new_position_100 = {rnd.Gauss(initial_position_100[0], 0.75*a), rnd.Gauss(initial_position_100[1], 0.75*a), rnd.Gauss(initial_position_100[2], 0.75*a)};
        initial_position_100 = Metropolis_100(initial_position_100, new_position_100, random_number, acceptance);
        //210
        vector<double> new_position_210 = {rnd.Gauss(initial_position_210[0], 1.87*a), rnd.Gauss(initial_position_210[1], 1.87*a), rnd.Gauss(initial_position_210[2], 1.87*a)};
        initial_position_210 = Metropolis_210(initial_position_210, new_position_210, random_number, acceptance);
    }
    
    //Open files (append) to write the coordinates
    ofstream write_pos_far_gauss_100("pos_far_gauss_100.dat", ios::app);
    ofstream write_pos_far_gauss_210("pos_far_gauss_210.dat", ios::app);

    if (write_pos_far_gauss_100.is_open() && write_pos_far_gauss_210.is_open()){
        for (int i = 0; i < N; i++){
            double final_position_100 = 0;
            int acceptance_100 = 0;
            double final_position_210 = 0;
            int acceptance_210 = 0;
            for (int j = 0; j  < L; j++){
                double random_number = rnd.Rannyu();

                //100
                vector<double> new_position_100 = {rnd.Gauss(initial_position_100[0], 0.75*a), rnd.Gauss(initial_position_100[1], 0.75*a), rnd.Gauss(initial_position_100[2], 0.75*a)};
                write_pos_far_gauss_100 << initial_position_100[0] << "\t" 
                                        << initial_position_100[1] << "\t" 
                                        << initial_position_100[2] << endl;
                initial_position_100 = Metropolis_100(initial_position_100, new_position_100, random_number, acceptance_100);
                
                //210
                vector<double> new_position_210 = {rnd.Gauss(initial_position_210[0], 1.87*a), rnd.Gauss(initial_position_210[1], 1.87*a), rnd.Gauss(initial_position_210[2], 1.87*a)};
                write_pos_far_gauss_210 << initial_position_210[0] << "\t" 
                                        << initial_position_210[1] << "\t" 
                                        << initial_position_210[2] << endl;
                initial_position_210 = Metropolis_210(initial_position_210, new_position_210, random_number, acceptance_210);
                
                final_position_100 += sqrt(initial_position_100[0]*initial_position_100[0] + initial_position_100[1]*initial_position_100[1] + initial_position_100[2]*initial_position_100[2]);
                final_position_210 += sqrt(initial_position_210[0]*initial_position_210[0] + initial_position_210[1]*initial_position_210[1] + initial_position_210[2]*initial_position_210[2]);
            }
            //cout << i << "      " << acceptance_100/double(L)*100 << endl;            disable the commeting if you want to see the acceptance rate (circa 50% for both)
            //cout << i << "      " << acceptance_210/double(L)*100 << endl;

            ave_100[i] = final_position_100/L;
            av2_100[i] = pow((final_position_100/L), 2); 
            ave_210[i] = final_position_210/L;
            av2_210[i] = pow((final_position_210/L), 2);    
        } 
        write_pos_far_gauss_100.close();
        write_pos_far_gauss_210.close();
    } else {
    cerr << "PROBLEM: Unable to open file for appending!" << endl;
    }

    //Computation of A and <A>^2 in function of N
    running_sum_100 = 0;
    running_sum2_100 = 0;
    running_sum_210 = 0;
    running_sum2_210 = 0;
    
    for (int i = 0; i < N; i++) {
        running_sum_100 += ave_100[i];
        running_sum2_100 += av2_100[i];
        running_sum_210 += ave_210[i];
        running_sum2_210 += av2_210[i];
        
        sum_prog_100[i] = running_sum_100 / (i+1);
        su2_prog_100[i] = running_sum2_100 / (i+1);
        err_prog_100[i] = error(sum_prog_100[i], su2_prog_100[i], i);
        sum_prog_210[i] = running_sum_210 / (i+1);
        su2_prog_210[i] = running_sum2_210 / (i+1);
        err_prog_210[i] = error(sum_prog_210[i], su2_prog_210[i], i);
    }

    //Create an output file for Python
    ofstream out_100_far_gauss("output_100_far_gauss.dat");
    if (out_100_far_gauss.is_open()){
        for (int i = 0; i < N; i++) {
            out_100_far_gauss << i << " " << sum_prog_100[i] << " " << err_prog_100[i] << endl;
        }
        out_100_far_gauss.close();
    } else cerr << "PROBLEM: Unable to open output_100_far_gauss.dat" << endl;
    ofstream out_210_far_gauss("output_210_far_gauss.dat");
    if (out_210_far_gauss.is_open()){
        for (int i = 0; i < N; i++) {
            out_210_far_gauss << i << " " << sum_prog_210[i] << " " << err_prog_210[i] << endl;
        }
        out_210_far_gauss.close();
    } else cerr << "PROBLEM: Unable to open output_210_far_gauss.dat" << endl;

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
