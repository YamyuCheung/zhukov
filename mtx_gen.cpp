#include <iostream>
#include <fstream>
#include <cstdlib>
#include <ctime>
#include <string>

int main(int argc, char* argv[]) {
    if (argc != 3) {
        std::cerr << "Usage: " << argv[0] << " <matrix size> <output file>" << std::endl;
        return 1;
    }

    int n = std::stoi(argv[1]);
    std::string output_file = argv[2];

    std::srand(std::time(0));
    std::ofstream result(output_file);

    if (!result.is_open()) {
        std::cerr << "Error opening file: " << output_file << std::endl;
        return 1;
    }

    result << n << std::endl;

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            int random_value = std::rand() % 100;
            result << random_value;
            if (j != n - 1) {
                result << ' ';
            }
        }
        result << '\n';
    }

    result.close();
    std::cout << "Random " << n << "x" << n << " matrix with size " << n << " saved to " << output_file << std::endl;

    return 0;
}