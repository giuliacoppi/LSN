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

    //Parameters
    int N_city = 34;            //Total number of cities to visit
    int N_individuals = 200;    //Total number of individuals in my population
    int N_generations = 500;    //Number of generations iterated from the first one
    double p_m = 0.08;          //Mutation probability
    double p_c = 0.75;          //Crossover probability

    //Create my cities
    vector<vec> cities_on_circle = city_circle(rnd, N_city);
    vector<vec> cities_on_square = city_square(rnd, N_city);

    vector<vec> city = cities_on_square;    //remember to change the output filename (line 67) when changing city!

    //Create the starting population
    vector<Individual> population = creation(N_city, N_individuals, rnd, city);
    //order the population
    sort(population.begin(), population.end(), [](const Individual& x, const Individual& y) { return x.cost < y.cost; } );


    //file for L(1) to plot
    ofstream out("output_L_1_square.dat");
    out << "#" 
        << setw(11) << "Generation"   
        << setw(53) << "Route" 
        << setw(53) << "L(1)" << endl;

    //file for <L(1)> to plot
    ofstream out_mean("output_L_1_mean_square.dat");
    out_mean << "#" 
             << setw(11) << "Generation" 
             << setw(12) << "<L(1)>" << endl;

    //New generations
    //sort(new_population.begin(), new_population.end(), [](const Individual& x, const Individual& y) { return x.cost < y.cost; });
    for (int i = 0; i < N_generations; i++){
        vector<Individual> new_population;
        while ((int) new_population.size() < N_individuals){
            Individual mother = population[selection(rnd, N_individuals)];
            Individual father = population[selection(rnd, N_individuals)];
            //Crossover, generates 2 children
            pair<vector<int>, vector<int>> children;
            if(rnd.Rannyu() < p_c){
                children = crossover(mother.chromosome, father.chromosome, N_city, rnd);
            }else {
                children = make_pair(mother.chromosome, father.chromosome); //No crossover, they stay the same
                }
            Individual son[2];
            son[0].chromosome = children.first;
            son[1].chromosome = children.second;
            //Mutate the chosen individual of the old population
            for (int k = 0; k < 2; k++){
                if(rnd.Rannyu() < p_m){
                    double p_m_select = rnd.Rannyu();   //Each son should have no more than one mutation
                    if(p_m_select < 0.25) son[k].chromosome = pair_permutation(son[k].chromosome, N_city, rnd);
                    else if(p_m_select < 0.5) son[k].chromosome = shift(son[k].chromosome, N_city, rnd);
                    else if(p_m_select < 0.75) son[k].chromosome = permutation(son[k].chromosome, N_city, rnd);
                    else if(p_m_select < 1) son[k].chromosome = inversion(son[k].chromosome, N_city, rnd);
                }
                son[k].cost = loss(son[k].chromosome, city);
                bool checked = check(son[k].chromosome);
                if(checked) new_population.push_back(son[k]);   
            }
        }
        sort(new_population.begin(), new_population.end(), [](const Individual& x, const Individual& y) { return x.cost < y.cost; } );
        population = new_population;

        //data for L(1) and 2D plot
        out << setw(8) << i << setw(8);
            for(int city : population[0].chromosome) {
                out << city << " "; 
            }
        out << setw(12) << population[0].cost << endl;

        //data for <L(1)>
        out_mean << setw(11) << i
                 << setw(12) << mean_L(population, N_individuals) << endl;

    }
    out.close();
    out_mean.close();
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
