#include <iostream>
#include <vector>
#include <algorithm>
#include <random>

using namespace std;

const int POPULATION_SIZE = 50; // population size
const int NUM_GENERATIONS = 1000; // number of generations
const double MUTATION_RATE = 0.1; // mutation rate
const double CROSSOVER_RATE = 0.5; // crossover rate
const double ELITE_RATE = 0.1; // elite rate

double rosenbrock(vector<double> individual) {
    double x = individual[0];
    double y = individual[1];
    return pow(1-x,2) + 100*pow(y-x*x,2);
}

vector<vector<double>> initializePopulation(int size, int dimension) {
    vector<vector<double>> population;

    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<> dis(-1, 1);

    for (int i = 0; i < size; i++) {
        vector<double> individual;
        for (int j = 0; j < dimension; j++) {
            individual.emplace_back(dis(gen));
        }
        population.emplace_back(individual);
    }

    return population;
}

vector<vector<double>> evolvePopulation(vector<vector<double>> population) {
    
    vector<vector<double>> newPopulation;
    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<> dis(0, 1);
    uniform_int_distribution<> idx(0, population.size() - 1);

    sort(population.begin(), population.end(), [](vector<double> a, vector<double> b) { return rosenbrock(a) < rosenbrock(b); });
    int eliteSize = POPULATION_SIZE * ELITE_RATE;
    vector<vector<double>> elite(population.begin(), population.begin() + eliteSize);

    while (newPopulation.size() < POPULATION_SIZE - eliteSize) {
        vector<double> parent1 = population[idx(gen)];
        vector<double> parent2 = population[idx(gen)];

        vector<double> offspring;
        if (dis(gen) < CROSSOVER_RATE) {
            int crossoverPoint = idx(gen);
            for (int i = 0; i < parent1.size(); i++) {
                if (i < crossoverPoint) {
                    offspring.push_back(parent1[i]);
                } else {
                    offspring.push_back(parent2[i]);
                }
            }
        } else {
            offspring = parent1;
        }

        for (int i = 0; i < offspring.size(); i++) {
            if (dis(gen) < MUTATION_RATE) {
                offspring[i] += dis(gen);
            }
        }

        newPopulation.push_back(offspring);
    }

    newPopulation.insert(newPopulation.end(), elite.begin(), elite.end());

    return newPopulation;
}


int main() {
    vector<vector<double>> population = initializePopulation(POPULATION_SIZE, 2);
    
    for (int i = 0; i < NUM_GENERATIONS; i++)  {
        population = evolvePopulation(population);
        vector<double> best = population[0];
        double bestFitness = rosenbrock(best);
        std::cout << "Generation: " << i << ", Best fitness: " << bestFitness << endl;
    }

return 0;

}