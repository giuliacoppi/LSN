#ifndef __FUNCTIONS_H__
#define __FUNCTIONS_H__

#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <iomanip>
#include <algorithm>
#include <armadillo>
#include "random.h"

using namespace std;
using namespace arma;

//little struct for every individual, this way I won't compute the loss function for each individual too many times when sorting the starting population
struct Individual {
    vector<int> chromosome;
    double cost;
};

// Functions' definitions
double error(double AV, double AV2, int n);
double loss(vector<int> chromosome, vector<vec>& city);
bool check(vector<int> chromosome);
double mean_L(const vector<Individual>& population, int N_individuals);

//City coordinates
vector<vec> city_circle(Random& rnd, int N_city);        //Contains (x, y) for each city (matrix N_cityx2)
vector<vec> city_square(Random& rnd, int N_city);        //Contains (x, y) for each city (matrix N_cityx2)


//Initial population creation and individual selection
vector<Individual> creation(int N_city, int N_individuals, Random& rnd, vector<vec>& city);
int selection(Random& rnd, int N_individuals);

//Genetic-mutation operator
vector<int> pair_permutation(vector<int>& chromosome, int N_city, Random& rnd);
vector<int> shift(vector<int>& chromosome, int N_city, Random& rnd);
vector<int> permutation(vector<int>& chromosome, int m, int start_1, int start_2, Random& rnd);
vector<int> inversion(vector<int>& chromosome, int N_city, Random& rnd);

//Crossover operator
pair<vector<int>, vector<int>> crossover(const vector<int>& mother, const vector<int>& father, int N_city, Random& rnd);




//////////////////////////////////////////////FUNCTIONS//////////////////////////////////////////////

//Calculates the statistical uncertainty (standard deviation of the mean)
double error(double AV, double AV2, int n) { 
    if (n == 0) {
        return 0;
    }
    return sqrt((AV2 - pow(AV, 2)) / n);
}

//Calculates the cost funciton
double loss(vector<int> chromosome, vector<vec>& city){
    double cost = 0.0;
    for (int i = 0; i < (int)chromosome.size() - 1; i++){
        int c1 = chromosome[i] - 1;
        int c2 = chromosome[i+1] - 1;
        cost += arma::norm((city[c1] - city[c2]));
    }
    int last = chromosome[chromosome.size() - 1] - 1;
    int first = chromosome[0] - 1;
    cost += arma::norm(city[last] - city[first]);
    
    return cost;
}

//Checs if the individual fulfills the bonds
bool check(vector<int> chromosome){
    bool approved = false;
    sort(chromosome.begin() + 1, chromosome.end());
    if(chromosome[0] == 1){
        for (int i = 1; i < (int)chromosome.size(); i++){
            if(chromosome[i] != i+1){
                return approved;
            }
            approved = true;
        }  
    }
    return approved;
}

double mean_L(const vector<Individual>& population, int N_individuals){
    double sum = 0.0;
    for (int i = 0; i < N_individuals/2 ; i++){
        sum += population[i].cost;
    }
    return (sum/((double)N_individuals/2.0));
}



//City coordinates
vector<vec> city_circle(Random& rnd, int N_city){
    ofstream out("output_circle_city.dat");
    out << "#" 
        << setw(11) << "city"   
        << setw(12) << "x" 
        << setw(12) << "y" << endl;
    vector<vec> city_positions;
    for(int i = 0; i < N_city; i++){
        double theta = rnd.Rannyu(0, 2*M_PI);
        vec city = { cos(theta), sin(theta) };
        city_positions.push_back(city);
        out << setw(11) << i   
            << setw(12) << city[0]
            << setw(12) << city[1] << endl;
    }
    return city_positions;
}


vector<vec> city_square(Random& rnd, int N_city){
    ofstream out("output_square_city.dat");
    out << "#" 
        << setw(11) << "city"   
        << setw(12) << "x" 
        << setw(12) << "y" << endl;
    vector<vec> city_positions;
    for(int i = 0; i < N_city; i++){
        vec city = { rnd.Rannyu(-1.0, 1.0), rnd.Rannyu(-1.0, 1.0)};
        city_positions.push_back(city);
        out << setw(11) << i   
            << setw(12) << city[0]
            << setw(12) << city[1] << endl;
    }
    return city_positions;
}






//Initial population creation
vector<Individual> creation(int N_city, int N_individuals, Random& rnd, vector<vec>& city){
    vector<Individual> population;
    for (int i = 0; i < N_individuals; i++){ 
        //Creation of the first chromosome, ordered from 1 to N_city
        vector<int> parent(N_city, 0);
        for (int i = 0; i < N_city; i++){
            parent[i] = i + 1;
        }  
        do{
            //Creation of a new chromosome
            for(int j = 0; j < 30; j++){
                parent = pair_permutation(parent, N_city, rnd);
            }
        }while(check(parent) == false);     //check if the new chromosome generated fulfills the bonds, otherwise create another one
        Individual ind;
        ind.chromosome = parent;
        ind.cost = loss(parent, city);

        population.push_back(ind);
    }
    return population;
}

//Individual selection
int selection(Random& rnd, int N_individuals){
    double random_number = rnd.Rannyu();
    int j = int(N_individuals * pow(random_number, 3));         //I have chosen p=3
    return j;
}




//Genetic mutations

//pair permutation
vector<int> pair_permutation(vector<int>& chromosome, int N_city, Random& rnd){
    bool ok = false;
    do{
        int city_1 = rnd.Rannyu(1, N_city - 1);
        int city_2 = rnd.Rannyu(1, N_city - 1);
        swap(chromosome[city_1], chromosome[city_2]);
        ok = check(chromosome);
    }while(ok == false);
    return chromosome;
}

//shifting m contiguous cities of n positions in the chromosome, m < N-1
vector<int> shift(vector<int>& chromosome, int N_city, Random& rnd){
    vector<int> chromosome_copy = chromosome;
    int L = chromosome.size() - 1;                  //cities I can move, excludes the fisrt one
    int n = int(rnd.Rannyu(1, N_city - 1));
    for (int i = 1; i < (int)chromosome.size(); i++){              
        int new_pos = (i - 1 + n) % L + 1;          //must use n > 0
        chromosome_copy[new_pos] = (chromosome[i]);
    }
    return chromosome_copy;    
}

//permutation among m contiguous cities in the chromosome, m < N/2
vector<int> permutation(vector<int>& chromosome, int N_city, Random& rnd){
    int m = rnd.Rannyu(2, (N_city - 1)/2 + 1);      //Starts from 2 because otherwise it's the same as pair_permutation
    int start_1 = int(rnd.Rannyu(1, N_city + 1 - 2*m));
    int start_2 = int(rnd.Rannyu(start_1 + m, N_city + 1 - m));
    for (int i = 0; i < m; i++){
        swap(chromosome[start_1 + i], chromosome[start_2 + i]);
    }
    return chromosome;
}

//Inversion of the order in which they appear in the path of m cities, m <= N
vector<int> inversion(vector<int>& chromosome, int N_city, Random& rnd){
    int index_start = int(rnd.Rannyu(1, N_city - 1));
    int m = int(rnd.Rannyu(2, N_city - index_start + 1));
    for (int i = 0; i < m/2; i++){
        swap(chromosome[index_start + i], chromosome[index_start + m - 1 - i]);
    } 
    return chromosome;
}


//Crossover function
pair<vector<int>, vector<int>> crossover(const vector<int>& mother, const vector<int>& father, int N_city, Random& rnd){
    vector<int> son_1;
    vector<int> son_2;
    int m = rnd.Rannyu(2, N_city - 1);      //2 and N_city - 1 because outerwise the cut wouldn'make sense and there would be no crossover
    for (int i = 0; i < m; i++){
        son_1.push_back(mother[i]);
        son_2.push_back(father[i]);
    }
    for (int i = 1; i < N_city; i++){
        bool is_missing_1 = true;
        bool is_missing_2 = true;

        for (int j = 0; j < m; j++){
            if(son_1[j] == father[i]){
                is_missing_1 = false;
                break;
            }
        }
        if(is_missing_1 == true) son_1.push_back(father[i]);

        for (int j = 0; j < m; j++){
            if(son_2[j] == mother[i]){
                is_missing_2 = false;
                break;
            }
        }
        if(is_missing_2 == true) son_2.push_back(mother[i]);
    }
    return make_pair(son_1, son_2);
}



/////////////////////////////////////////////////////////////////////////////////////////////////////

#endif

