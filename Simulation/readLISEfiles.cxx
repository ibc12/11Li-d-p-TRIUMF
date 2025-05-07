#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>

void readLISEfiles() {
    std::ifstream file("../LISE files/11Li_silicon.txt");
    if (!file.is_open()) {
        std::cerr << "No se pudo abrir el archivo." << std::endl;
        return;
    }

    std::string line;
    // Saltar las 3 primeras líneas
    for (int i = 0; i < 3; ++i)
        std::getline(file, line);

    std::vector<double> col1, col4, col12;

    while (std::getline(file, line)) {
        std::istringstream ss(line);
        std::string token;
        std::vector<std::string> tokens;

        // Separar por espacios/tabs
        while (ss >> token) {
            tokens.push_back(token);
        }

        if (tokens.size() >= 12) {
            col1.push_back(std::stod(tokens[0]));
            col4.push_back(std::stod(tokens[3]));
            col12.push_back(std::stod(tokens[11]));
        }
    }

    file.close();

    // Ejemplo de impresión de los primeros valores
    for (size_t i = 0; i < col1.size(); ++i) {
        std::cout << "Fila " << i+1 << ": "
                  << "Col1 = " << col1[i] << ", "
                  << "Col4 = " << col4[i] << ", "
                  << "Col12 = " << col12[i] << std::endl;
    }
}