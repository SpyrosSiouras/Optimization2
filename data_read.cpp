#include <fstream>
#include <iostream>

int main() {
    std::ifstream file("data.txt");  // Open the input file
    float a, b;                       // Variables to store the data

    // Read data from the file
    while (file >> a >> b) {
        // Do something with the data
        std::cout << a << " " << b << std::endl;
    }

    file.close();  // Close the file

    return 0;
}
