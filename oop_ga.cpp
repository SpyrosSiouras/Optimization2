#include "geneticAlgorithm.hpp"
#include "csvWriter.cpp"


int main() {

    std::string csvFile = "resultsGeneticAlgorithmNoNoise.csv";
    // RosenbrockGA ga(1000, 3000, 0.05, 0.7, 0.1, 0.00005);
    //             ga(POPULATION_SIZE, NUMBER_OF_GENERATIONS, MUTATION_RATE, CROSSOVER_RATE, ELITE_RATE, ACCURACY);
    
    for (int i=0; i < 30; i++) {
        
        PredatorPreyGA ga(800, 10000, 0.05, 0.7, 0.1, 0.035);
        // ga.printData();
        ga.run();

        if(!fileExists(csvFile))
            writeCsvFile(csvFile, "fitness", "r", "K", "s", "u", "v", "functionCalls");
        double r = ga.overallBestIndividual[0];
        double K = ga.overallBestIndividual[1];
        double s = ga.overallBestIndividual[2];
        double u = ga.overallBestIndividual[3];
        double v = ga.overallBestIndividual[4];
        if (!writeCsvFile(csvFile, ga.overallBestFitness,r, K, s, u, v, ga.functionCalls)) {
            std::cerr << "Failed to write to file: " << csvFile << "\n";
            }
    }
    return 0;
}

