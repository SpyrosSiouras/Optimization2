#include "geneticAlgorithm.hpp"

bool fileExists(std::string& fileName) {
    return static_cast<bool>(std::ifstream(fileName));
}

template <typename filename, typename T1, typename T2, typename T3,typename T4 ,typename T5 ,typename T6,typename T7>
bool writeCsvFile(filename &fileName, T1 column1, T2 column2, T3 column3, T4 column4, T5 column5, T6 column6, T7 column7) {
    std::lock_guard<std::mutex> csvLock(logMutex);
    std::fstream file;
    file.open (fileName, std::ios::out | std::ios::app);
    if (file) {
        file << column1 << ",";
        file << column2 << ",";
        file << column3 << ",";
        file << column4 << ",";
        file << column5 << ",";
        file << column6 << ",";
        file << column7 << " ";

        file <<  std::endl;
        return true;
    } else {
        return false;
    }
}

int main() {

    std::string csvFile = "resultsGeneticAlgorithm.csv";
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

