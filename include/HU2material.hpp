#pragma once
    
#include <unordered_map>
#include <tuple>
#include <string>
#include <fstream>     
#include <sstream>     
#include <iostream>   
#include <filesystem> 
#include "error_logging.hpp"


struct HU2MaterialTable {
    // Map from HU → (RSP, X, density)
    std::unordered_map<int, std::tuple<float, float, float>> table;

    HU2MaterialTable(const std::string file_path);

    float RSP(int HU_val) const;
    float X(int HU_val) const;
    float density(int HU_val) const;
};
