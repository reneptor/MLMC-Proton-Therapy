#include "HU2material.hpp"


HU2MaterialTable::HU2MaterialTable(const std::string file_path){

    std::ifstream file(file_path);  
    assertWithMessage(file.is_open(), "Failed to open file: " + file_path);

    std::string line;
    while (std::getline(file, line)) {
        if (line.empty() || line[0] == '#') {
            continue;
        }
        std::stringstream ss(line);
        std::string item;
        int HU_val {};
        float RSP {}, density {}, X {};

        // Parse the three comma-separated values
        if (std::getline(ss, item, ',')) HU_val = std::stoi(item);      // Houndsfield unit
        if (std::getline(ss, item, ',')) RSP = std::stof(item);         // Relative stopping power to water
        if (std::getline(ss, item, ',')) density = std::stof(item);     // Density [g/cm3]
        if (std::getline(ss, item, ',')) X = std::stof(item);           // Radiation length [mm]

        // Add to map
        table[HU_val] = std::make_tuple(RSP, X, density);
    }
    file.close();
}


float HU2MaterialTable::RSP(int HU_val) const {
    return std::get<0>(table.at(HU_val));
}


float HU2MaterialTable::density(int HU_val) const {
    return std::get<2>(table.at(HU_val));
}


float HU2MaterialTable::X(int HU_val) const {
    return std::get<1>(table.at(HU_val));
}


