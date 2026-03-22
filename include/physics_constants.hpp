#pragma once

#include <cmath>


namespace physics_constants {

constexpr float ELECTRON_MASS = 0.511;             // [MeV/c^2]
constexpr float PROTON_MASS = 938.272;             // [MeV/c^2]
constexpr double PI = 3.141592653589793;
constexpr double C = 3e8;                          // [m/s]
constexpr double ELECTRON_CHARGE = 1.602e-19;      // [C]
constexpr double ELECTRON_MASS_KG = 9.109e-31;     // [kg]
constexpr double ELECTRON_RADIUS = 2.817e-15;      // [m]
constexpr double I_POT_WATER = 74.0;               // [eV]
constexpr double PERMITTIVITY = 8.854e-12;         // [F/m]
constexpr double WATER_E_RHO = 3.34e23;            // [electrons/m^3]
constexpr double ETA_CONSTANT = 
2*PI*ELECTRON_MASS_KG*ELECTRON_RADIUS*WATER_E_RHO; // [C^2/m^2]

constexpr double SP_CONSTANT_WATER = 
    (0.001 / ELECTRON_CHARGE) *
    ((4.0 * PI) / (ELECTRON_MASS_KG * C * C)) * WATER_E_RHO *
    std::pow((ELECTRON_CHARGE * ELECTRON_CHARGE) / (4.0 * PI * PERMITTIVITY), 2); // [MeV/mm c^2]


inline float gamma(float E){
    return 1+(E/PROTON_MASS);
}


inline float MSAngleHighland(float r){ // Highland formula (squared), disregarding the energy dependance for now
    return r*(1 + 0.038*std::log(r))*(1 + 0.038*std::log(r));
}


inline float HighlandEnergyDependance(float E){
    //return std::pow(E+PROTON_MASS, 2)/(std::pow(E+PROTON_MASS, 4) - std::pow(PROTON_MASS, 4));
    //return (E+PROTON_MASS)/std::sqrt(std::pow(E+PROTON_MASS, 4) - std::pow(PROTON_MASS, 4));
    return (E+PROTON_MASS)/(E*E + 2*E*PROTON_MASS);
}


}