#include <iostream>
#include <mutex>
#include <fstream>

std::mutex logMutex;

bool fileExists(std::string& fileName) {
    return static_cast<bool>(std::ifstream(fileName));
}

template <typename filename, typename T1, typename T2, typename T3>
bool writeCsvFile(filename &fileName, T1 column1, T2 column2, T3 column3) {
    std::lock_guard<std::mutex> csvLock(logMutex);
    std::fstream file;
    file.open (fileName, std::ios::out | std::ios::app);
    if (file) {
        file << "\"" << column1 << "\",";
        file << "\"" << column2 << "\",";
        file << "\"" << column3 << "\"";
        file <<  std::endl;
        return true;
    } else {
        return false;
    }
}

int main() {
    std::string csvFile = "results.csv";
    std::string naam = "Hallo";
    if(!fileExists(csvFile))
        writeCsvFile(csvFile, "header1", "header2", "header3");

    for (int i = 1; i < 10; ++i) {
        if (!writeCsvFile(csvFile, i, naam, static_cast<float>(i * 3.5))) {
            std::cerr << "Failed to write to file: " << csvFile << "\n";
        }
    }
    return 0;
}
