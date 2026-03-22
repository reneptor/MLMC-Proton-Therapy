#include <string>   // for std::string and std::stoi
#include <cstdlib>  // optional, for std::exit if you add error checks
#include <random>

#include "phantom.hpp"
#include "grids.hpp"
#include "HU2material.hpp"
#include "pencilbeam.hpp"


int main(int argc, char *argv[]){
    // Usage: ./executable (phantom_path) (dim_x) (dim_y) (dim_z)

    const std::string phantom_path = "phantom.dat"; //argv[1];
    const unsigned int phantom_n_x = 100; //std::stoi(argv[2]);
    const unsigned int phantom_n_y = 100; //std::stoi(argv[3]);
    const unsigned int phantom_n_z = 1000; //std::stoi(argv[4]);
    // const std::string phantom_path = argv[1];
    // const unsigned int phantom_n_x = std::stoi(argv[2]);
    // const unsigned int phantom_n_y = std::stoi(argv[3]);
    // const unsigned int phantom_n_z = std::stoi(argv[4]);
    auto phantom = loadPhantom(phantom_path, phantom_n_x, phantom_n_y, phantom_n_z);

    HU2MaterialTable hu2mat = HU2MaterialTable("HU2material.csv");
    

    unsigned int n_primaries = 100000;
    unsigned int level = 0;
    float initial_step_factor = 0.02;

    MediumGrid mediumGrid = MediumGrid(&hu2mat);
    mediumGrid.loadPhantom(phantom);
    ScoringGrid primaryScoringGrid = ScoringGrid(level, n_primaries, true);
    ScoringGrid secondaryScoringGrid = ScoringGrid(level, n_primaries, false);
    
    std::mt19937 gen(123456);  // Magic number (for now)
    PencilBeam<std::mt19937> pb0(n_primaries, initial_step_factor,
    &primaryScoringGrid, &secondaryScoringGrid, &mediumGrid, 1.0, 'z', gen);
    pb0.initParticleHistories(false); // Init shadows set to false 
    pb0.simulatePrimaryBeam();

    return 0;
}