#include <iostream>
#include <vector>
#include <algorithm>
#include <random>
#include <map>
#include <fstream>
#include <mutex>

const int MAX_FUNCTION_CALLS = 1000000;
using namespace std;


class GeneticAlgorithm {
public:
    GeneticAlgorithm(int populationSize, int numGenerations, double mutationRate, double crossoverRate, double eliteRate, double accuracy) :
        populationSize(populationSize), numGenerations(numGenerations), mutationRate(mutationRate), crossoverRate(crossoverRate), eliteRate(eliteRate), accuracy(accuracy) {}

    virtual ~GeneticAlgorithm() {}
    virtual void run();

    int functionCalls;

    vector<double> overallBestIndividual;
    double overallBestFitness = 1000.0;

protected:
    virtual vector<double> initializeIndividual() = 0;
    virtual double fitness(vector<double>& individual) = 0;
    virtual vector<vector<double>> evolvePopulation() = 0;
    vector<vector<double>> initializePopulation();
    void printBest(vector<double> best) {
        for (auto& x: best) {
            cout<< x << " ";
        }
    }


    vector<vector<double>> population;
    int populationSize;
    int numGenerations;
    double mutationRate;
    double crossoverRate;
    double eliteRate;
    double accuracy;
    int badLocalMinCounter;

};

void GeneticAlgorithm:: run() {
        population = initializePopulation();
        badLocalMinCounter = 0;
        for (int i = 0; i < numGenerations; i++) {
            population = evolvePopulation();
            vector<double> best = population[0];
            double bestFitness = fitness(best);
            if (bestFitness < overallBestFitness){
                overallBestFitness = bestFitness;
                overallBestIndividual = best;
            }
            cout << "Generation: " << i << ", Best fitness: " << bestFitness <<", with coefficients: ";
            printBest(best); 
            cout <<endl;
            if (bestFitness < accuracy) {
                cout << "Reached accuracy of: " << accuracy  << " in " << i << " generations, with " << functionCalls <<" function calls! " << endl;
                return;
            }

            if (functionCalls > MAX_FUNCTION_CALLS) {
                cout << "Reached MAX function calls of: " << MAX_FUNCTION_CALLS  << "! Overall Best found: " << overallBestFitness << ", with coefficients: ";
                printBest(overallBestIndividual);
                cout << "!"  << endl;
                return;
            }

            if (bestFitness > 5)
                populationSize = 800;

            else if (bestFitness < 5 && bestFitness > 1)
                populationSize = 400;

            else if (bestFitness < 1 && bestFitness > 0.1)
                populationSize = 200;

            else if (bestFitness < 0.1)
                populationSize = 100;


            if ( populationSize > 100)
                badLocalMinCounter++;
        }

        cout << "Overall Best: " << overallBestFitness << ", with coefficients: ";
        printBest(overallBestIndividual);
        cout <<" and "<< functionCalls <<" function calls!" <<endl;

    }


vector<vector<double>> GeneticAlgorithm:: initializePopulation() {
        
        vector<vector<double>> population;
        for (int i = 0; i < populationSize; i++)
        {
            population.push_back(initializeIndividual());
        }
        
        return population;
    }

// Rosenbrock
class RosenbrockGA : public GeneticAlgorithm {
public:
    RosenbrockGA(int populationSize, int numGenerations, double mutationRate, double crossoverRate, double eliteRate, double accuracy) :
        GeneticAlgorithm(populationSize, numGenerations, mutationRate, crossoverRate, eliteRate, accuracy) {}

protected:
    vector<double>initializeIndividual();
    double fitness(vector<double>individual);
    vector<vector<double>> evolvePopulation();
};


vector<double>RosenbrockGA:: initializeIndividual() {
        vector<double>individual(2);
        random_device rd;
        mt19937 gen(rd());
        uniform_real_distribution<> dis(-3, 3);
        individual[0] = dis(gen);
        individual[1] = dis(gen);
        return individual;
    }


double RosenbrockGA:: fitness(vector<double>individual) {
        double x = individual[0];
        double y = individual[1];
        return pow(1-x,2) + 100*pow(y-x*x,2);
    }

vector<vector<double>> RosenbrockGA:: evolvePopulation() {
        
        vector<vector<double>> newPopulation;
        random_device rd;
        mt19937 gen(rd());
        uniform_real_distribution<> dis(0, 1);
        uniform_real_distribution<> crRate(-0.25, 1.25);
        uniform_int_distribution<> idx(0, population.size() - 1);

        sort(population.begin(), population.end(), [this] (vector<double> a, vector<double> b) { return fitness(a) < fitness(b); });
        int eliteSize = populationSize * eliteRate;
        vector<vector<double>> elite(population.begin(), population.begin() + eliteSize);

        while (newPopulation.size() < populationSize - eliteSize) {
            
            vector<double>parent1 = population[idx(gen)];
            vector<double>parent2 = population[idx(gen)];

            vector<double>offspring = parent1;

            if (dis(gen) < crossoverRate) {
                for (int i = 0; i < 2; i++) {
                    double crossoverPoint = crRate(gen);
                    offspring[i] = crossoverPoint * parent1[i] + (1 - crossoverPoint) * parent2[i];
                }
            }

            for (int i = 0; i < 2; i++) {
                if (dis(gen) < mutationRate) {
                    double a = 1;
                    while (offspring[i] + a > 6 || offspring[i] - a < -6)
                        a /= 2;
                    uniform_real_distribution<> mutationRate(-a, a);
                    offspring[i] += mutationRate(gen);
                }
            }
            newPopulation.push_back(offspring);
        }
        newPopulation.insert(newPopulation.end(), elite.begin(), elite.end());
        return newPopulation;
}

class PredatorPreyGA : public GeneticAlgorithm {
public:
    PredatorPreyGA(int populationSize, int numGenerations, double mutationRate, double crossoverRate, double eliteRate, double accuracy) :
        GeneticAlgorithm(populationSize, numGenerations, mutationRate, crossoverRate, eliteRate, accuracy) {
            functionCalls = 0;
            loadData();
        }
    
    
    void printData(){
        for (int i=0; i<100; i++)
        {
          cout << realPValues[i] << " " << realQValues[i] << std::endl;
        }

    }

   
protected:
    vector<double>initializeIndividual();
    double fitness(vector<double>& individual);
    vector<vector<double>> evolvePopulation();

private: 

    double calculateAbsoluteRelativeError(double, double);
    void loadData();
    vector<double> realPValues;
    vector<double> realQValues;

    double r, K, s; 
    double u, v; 


};

void PredatorPreyGA:: loadData() {

    std::ifstream file("data.txt");  // Open the input file
    float a, b;                       // Variables to store the data

    // Read data from the file
    while (file >> a >> b) {
        // Do something with the data
        realPValues.emplace_back(a);
        realQValues.emplace_back(b);
    }

    file.close();  // Close the file
}

double PredatorPreyGA:: calculateAbsoluteRelativeError(double prediction, double real_value){
    return fabs(prediction - real_value ) / real_value;
}


vector<double>PredatorPreyGA:: initializeIndividual() {
        vector<double>individual(5);
        random_device rd;
        mt19937 gen(rd());
        uniform_real_distribution<> dis(0.1, 3);

        for (int i = 0; i < 5; i++){
            individual[i] = dis(gen);
        }

        return individual;
    }

double PredatorPreyGA:: fitness(vector<double>& individual) {

        // if (badLocalMinCounter > 250)
        // {
        //     individual[0] -= 0.08;
        //     individual[4] += 0.08;
        // }
        double r = individual[0];
        double K = individual[1];
        double s = individual[2];
        double u = individual[3];
        double v = individual[4];


        double currentP = realPValues[0];
        double currentQ = realQValues[0];
        double deltaP = (r * (1 - realPValues[0] / K) - s*realQValues[0]) * realPValues[0];
        double deltaQ = (-u + v*realPValues[0])*realQValues[0];

        vector<double> errorsOfP;
        vector<double> errorsOfQ;
        
        for (int i = 1; i < 101; i++){
            currentP = currentP + deltaP;
            currentQ = currentQ + deltaQ;
            double errorP = calculateAbsoluteRelativeError(currentP, realPValues[i]);
            double errorQ = calculateAbsoluteRelativeError(currentQ, realQValues[i]);
            
            // std::cout<<"Error of P: " << errorP << ", " << "Error of Q: " << errorQ << std::endl;

            errorsOfP.emplace_back(errorP);
            errorsOfQ.emplace_back(errorQ);

            deltaP = (r * (1 - currentP / K) - s*currentQ) * currentP;
            deltaQ = (-u + v*currentP)*currentQ;
        }

        double maxErrorP = *max_element(errorsOfP.begin(), errorsOfP.end());
        double maxErrorQ = *max_element(errorsOfQ.begin(), errorsOfQ.end());

        errorsOfP.clear();
        errorsOfQ.clear();
        // cout<< maxErrorP << maxErrorQ << endl;
        functionCalls++;
        return max(maxErrorP, maxErrorQ);
    }

vector<vector<double>> PredatorPreyGA:: evolvePopulation() {
        
        vector<vector<double>> newPopulation;
        random_device rd;
        mt19937 gen(rd());
        uniform_real_distribution<> dis(0, 1);
        uniform_real_distribution<> crRate(-0.5, 1.5);
        uniform_int_distribution<> idx(0, population.size() - 1);
        map<vector<double>, double> fitnessScore;
        for (auto &x: population){
            fitnessScore.insert({ x, fitness(x)});
        }
        sort(population.begin(), population.end(), [this, fitnessScore](vector<double>a, vector<double>b) { return fitnessScore.at(a) < fitnessScore.at(b); });
        int eliteSize = populationSize * eliteRate;
        vector<vector<double>> elite(population.begin(), population.begin() + eliteSize);

        while (newPopulation.size() < populationSize - eliteSize) {
            
            vector<double>parent1 = population[idx(gen)];
            vector<double>parent2 = population[idx(gen)];

            vector<double>offspring = parent1;

            if (dis(gen) < crossoverRate) {
                for (int i = 0; i < 5; i++) {
                    double crossoverPoint = crRate(gen);
                    offspring[i] = crossoverPoint * parent1[i] + (1 - crossoverPoint) * parent2[i];

                    if (offspring[i] > 3) {
                        offspring[i] = 3;
                    }

                    if (offspring[i] < 0.1) {
                        offspring[i] = 0.1;
                    }
                }
            }

            for (int i = 0; i < 5; i++) {
                if (dis(gen) < mutationRate) {
                    int coef = 10;
                    if (badLocalMinCounter > 250) {
                        coef = 100 * sqrt(numGenerations)* log(numGenerations);
                        badLocalMinCounter = 0;
                    }
                    // double a = min(pow(5, bestOfGeneration) / (sqrt(numGenerations) * log(numGenerations)), 3.0);
                    double a = (1.5 * coef) / (sqrt(numGenerations)* log(numGenerations));
                    // double a = min(pow(1.5, bestOfGeneration), 1.5) / log(numGenerations);
                    int count = 0;
                    while (offspring[i] + a > 3 || offspring[i] - a < 0.1) {
                        if (count > 5){
                            break;
                        }
                        a /= 2;
                        cout << "bisect" << endl;
                        count++;
                    }
                    uniform_real_distribution<> mutationRate(-a, a);

                    offspring[i] += mutationRate(gen);

                    if (offspring[i] > 3) {
                        offspring[i] = 3;
                    }

                    if (offspring[i] < 0.1)
                    {
                        offspring[i] = 0.1;
                    }
                }
            }
            newPopulation.push_back(offspring);
        }
        newPopulation.insert(newPopulation.end(), elite.begin(), elite.end());
        return newPopulation;
}




