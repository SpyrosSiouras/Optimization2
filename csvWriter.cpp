#include <mutex>
#include <fstream>

std::mutex logMutex;

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