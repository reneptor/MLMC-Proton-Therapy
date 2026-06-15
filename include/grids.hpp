#pragma once

#include <array>
#include <cmath>
#include "HU2material.hpp"
#include "simulation_constants.hpp"
#include "error_logging.hpp"


namespace gc = grid_constants;


struct MediumGrid{

    size_t n_x, n_y, n_z;
    HU2MaterialTable* hu2mat;
    int16_t* grid;

    MediumGrid(HU2MaterialTable* hu2mat, size_t n_x=gc::GRID_N_X, size_t n_y=gc::GRID_N_Y, size_t n_z=gc::GRID_N_Z);
    ~MediumGrid();	

    // void loadPhantom(const Phantom phantom);
    const float RSP(unsigned int i, unsigned int j, unsigned int k) const; 
    const float density(unsigned int i, unsigned int j, unsigned int k) const; 
    const float X(unsigned int i, unsigned int j, unsigned int k) const; 
};


struct ScoringGrid{
    unsigned int level {};
    //std::array<std::array<std::array<float, gc::GRID_N_X>, gc::GRID_N_Y>, gc::GRID_N_Z> scoreGrid {{{{{0.0f}}}}};

    size_t n_x, n_y, n_z;
    unsigned int n_samples {};
    float* grid;
    float* squareGrid;
    const MediumGrid* mediumGrid;

    ScoringGrid(unsigned int level, size_t n_x=gc::GRID_N_X, size_t n_y=gc::GRID_N_Y, size_t n_z=gc::GRID_N_Z, const MediumGrid* mediumGrid=nullptr);
    ~ScoringGrid();

    void initSquareGrid(); 
    void assignMediumGrid(const MediumGrid* mediumGrid);
    void depositEnergyLocal(unsigned int i, unsigned int j, unsigned int k, float dE);
    void addSample();
    std::unique_ptr<ScoringGrid> getDose() const;
    std::unique_ptr<ScoringGrid> getDoseVariance() const;
    float getVoxel(unsigned int i, unsigned int j, unsigned int k) const;
    float getVoxelVariance(unsigned int i, unsigned int j, unsigned int k) const;
    float getTotalVariance() const;
    float getVoxelDose(unsigned int i, unsigned int j, unsigned int k) const;
    float getVoxelDoseVariance(unsigned int i, unsigned int j, unsigned int k) const;
    float getTotalDoseVariance() const; 
    void setVoxel(float val, unsigned int i, unsigned int j, unsigned int k);
    void addVoxel(float val, unsigned int i, unsigned int j, unsigned int k);
};



