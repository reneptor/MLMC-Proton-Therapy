#pragma once

#include "grids.hpp"
#include "treatmentplan.hpp"
#include <random>
#include <vector>
#include <map>
#include <memory>
#include <string>
#include <cstdint>    
#include <string>      
#include <fstream>     
#include <stdexcept>   
#include <filesystem>
#include "error_logging.hpp"


template <typename RandomGen>
struct MLMCprotons {
 
    MediumGrid* mediumGrid;
    ScoringGrid* combinedScoringGrid;
    HU2MaterialTable hu2mat; 
    RandomGen gen;
    // const unsigned int maxThreads = std::thread::hardware_concurrency();
    static const unsigned int maxThreads;
    std::vector<unsigned int> levels;
    std::map<unsigned int, TreatmentPlan<RandomGen>*> treatmentPlans;
    
    MLMCprotons(unsigned int seed, const std::string tablePath);
    ~MLMCprotons();
    void loadPhantom(const std::string phantomPath, size_t n_x, size_t n_y, size_t n_z);
    void addTreatmentPlan(unsigned int level);
    void addPencilBeam(unsigned int level, float nPrimShare, double initialStep, 
               char entranceDir = dc::ENTRANCE_DIR, float beamWidth = dc::BEAM_WIDTH, 
               double dir_x = dc::DIR_X, double dir_y = dc::DIR_Y, double dir_z = dc::DIR_Z,
               double x_0 = dc::ENTRY_X, double y_0 = dc::ENTRY_Y, double z_0 = dc::ENTRY_Z,
               double E = dc::BEAM_ENERGY, double spread_E = dc::BEAM_ENERGY_SPREAD, float alpha=dc::ALPHA);
    void simulateTreatmentPlan(unsigned int level, unsigned int numPrimaries, unsigned int numThreads = maxThreads);
    void renderCombinedScoringGrid();
    std::unique_ptr<ScoringGrid> yieldDoseAtLevel(unsigned int level);
    std::unique_ptr<ScoringGrid> yieldDoseCombined();

};



