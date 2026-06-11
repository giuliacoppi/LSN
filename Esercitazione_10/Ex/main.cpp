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
#include "mpi.h"
#include "random.h"
#include "functions.h"

using namespace std;


int main (int argc, char *argv[]){

    //Activates the parallel environment
    MPI_Init(&argc, &argv);

    int size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    Random rnd;
    int seed[4];
    int p1, p2;
    ifstream Primes("Primes");
    if (Primes.is_open()){
        int pippo, pluto;
        for (int i = 0; i < rank; i++){
            Primes >> pippo >> pluto ;
        }
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

    //Parameters
    int N_city = 110;            //Total number of cities to visit
    int N_individuals = 500;     //Total number of individuals in my population
    int N_generations = 1000;    //Number of generations iterated from the first one
    int N_migrations = 50;       //Number of generations created between each migration
    int N_migrants = 30;         //Number of migrants in each migration
    double p_m = 0.08;           //Mutation probability
    double p_c = 0.75;           //Crossover probability

    double x_city[110];
    double y_city[110];

    if(rank == 0){
        ifstream file("cap_prov_ita.dat");
        if (file.is_open()) {
            for (int i = 0; i < N_city; i++){
                file >> x_city[i] >> y_city[i];
            }  
            file.close();
        }else{ cerr << "Problem: unable to open cap_prov_ita.dat" << endl; }
    }

    MPI_Bcast(x_city, N_city, MPI_DOUBLE, 0, MPI_COMM_WORLD); //Send to other continents the x coords of the cities
    MPI_Bcast(y_city, N_city, MPI_DOUBLE, 0, MPI_COMM_WORLD); //Send to other continents the y coords of the cities

    //Saves the cities in a vector<vec> city, this way my functions will still work
    vector<vec> city(N_city);
    for (int i = 0; i < N_city; i++){
        city[i].set_size(2);        //each armadillo vec has two components (x, y)
        city[i](0) = x_city[i];
        city[i](1) = y_city[i];
    }

    ofstream out_best;
    ofstream out_mean_best;

    //opened only by rank = 0
    if (rank == 0) {
        out_best.open("L_global_best.dat");
        if(out_best.is_open()){
            out_best << "#" 
                     << setw(11) << "Generation"   
                     << setw(12) << "L(1)" << endl;
        }
        
        out_mean_best.open("L_mean_global_best.dat");
        if(out_mean_best.is_open()){
            out_mean_best << "#" 
                          << setw(11) << "Generation" 
                          << setw(12) << "<L(1)>" << endl;
        }
    }
    

    //Create the starting population
    vector<Individual> population = creation(N_city, N_individuals, rnd, city);
    //order the population
    sort(population.begin(), population.end(), [](const Individual& x, const Individual& y) { return x.cost < y.cost; } );

    //New generations
    for (int i = 0; i < N_generations; i++){
        vector<Individual> new_population;
        while ((int) new_population.size() < N_individuals){
            evolve(population, new_population, city, rnd, N_individuals, N_city, p_c, p_m);
        }
        if((int)new_population.size() > N_individuals) { new_population.resize(N_individuals); }        //Protection in case the new_population has size N_individuals + 1 (in the evolve function I create two individuals to check, I will never know how may (0, 1, 2) have actualy benn pushed back inside the new population)
        sort(new_population.begin(), new_population.end(), [](const Individual& x, const Individual& y) { return x.cost < y.cost; } );
        population = new_population;

        //Migration
        if (i != 0 && i % N_migrations == 0){
         //   migration(N_city, city, size, rank, rnd, population, N_individuals, N_migrants);
        }

        double local_best_cost = population[0].cost;
        double global_best_cost = 0.0;

        double local_mean = mean_L(population, N_individuals);
        double global_mean = 0.0;

        //All continents call it simultaneously
        MPI_Reduce(&local_best_cost, &global_best_cost, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
        MPI_Reduce(&local_mean, &global_mean, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

        if(rank == 0){
            out_best << setw(8) << i << setw(12) << global_best_cost << endl;
            out_mean_best << setw(11) << i
                          << setw(12) << global_mean/size << endl;
        }

    }

    //Print the best route
    double final_local_cost = population[0].cost;
    vector<double> all_final_costs;
    vector<int> all_final_routes;
    if (rank == 0) {
        all_final_costs.resize(size);
        all_final_routes.resize(size * N_city);
    }
    //gather all costs in rank 0
    MPI_Gather(&final_local_cost, 1, MPI_DOUBLE, 
               all_final_costs.data(), 1, MPI_DOUBLE, 
               0, MPI_COMM_WORLD);
    //gather all routes in rank 0
    MPI_Gather(population[0].chromosome.data(), N_city, MPI_INT, 
               all_final_routes.data(), N_city, MPI_INT, 
               0, MPI_COMM_WORLD);

    if (rank == 0) {
        int best_rank_index = 0;
        double min_cost = all_final_costs[0];
    
        for (int j = 1; j < size; j++) {
            if (all_final_costs[j] < min_cost) {
                min_cost = all_final_costs[j];
                best_rank_index = j; 
            }
        }
        ofstream out_route("best_route.dat");
        if (out_route.is_open()) {
            int start_index = best_rank_index * N_city;
            for (int i = 0; i < N_city; i++) {
                out_route << all_final_routes[start_index + i] << endl;
            }
            out_route.close();
        }
    }


    out_best.close();
    out_mean_best.close();

    MPI_Finalize();
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
