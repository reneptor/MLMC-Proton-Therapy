#include "grids.hpp"


namespace cc = control_constants;
namespace gc = grid_constants;


MediumGrid::MediumGrid(HU2MaterialTable* hu2mat, size_t n_x, size_t n_y, size_t n_z)
    : hu2mat(hu2mat), n_x(n_x), n_y(n_y), n_z(n_z) {
    grid = new int16_t[n_x*n_y*n_z]();
}


MediumGrid::~MediumGrid() {
    delete[] grid;
}


const float MediumGrid::RSP(unsigned int i, unsigned int j, unsigned int k) const {
    return hu2mat->RSP(grid[i*n_y*n_z + j*n_z + k]);
}


const float MediumGrid::density(unsigned int i, unsigned int j, unsigned int k) const {
    return hu2mat->density(grid[i*n_y*n_z + j*n_z + k]);
}


const float MediumGrid::X(unsigned int i, unsigned int j, unsigned int k) const {
    return hu2mat->X(grid[i*n_y*n_z + j*n_z + k]);
}


ScoringGrid::ScoringGrid(unsigned int level, size_t n_x, size_t n_y, size_t n_z)
    : level(level), n_x(n_x), n_y(n_y), n_z(n_z) {
    grid = new float[n_x*n_y*n_z]();
}


void ScoringGrid::depositEnergyLocal(unsigned int i, unsigned int j, unsigned int k, float dE) {
    assertWithMessage(std::abs(dE) < cc::ENERGY_DEPOSIT_THRESHOLD, "Energy deposit exceeds threshold");
    assertWithMessage(i < n_x && j < n_y && k < n_z, "Grid index out of bounds");
    // assertWithMessage(dE > 0, "Energy deposit must be positive");
    grid[i*n_y*n_z + j*n_z + k] += dE; 
}


float ScoringGrid::getVoxel(unsigned int i, unsigned int j, unsigned int k) {
    assertWithMessage(i < n_x && j < n_y && k < n_z, "Grid index out of bounds");
    return grid[i*n_y*n_z + j*n_z + k];
}


void ScoringGrid::setVoxel(float val, unsigned int i, unsigned int j, unsigned int k) {
    assertWithMessage(i < n_x && j < n_y && k < n_z, "Grid index out of bounds");
    grid[i*n_y*n_z + j*n_z + k] = val;
}


void ScoringGrid::addVoxel(float val, unsigned int i, unsigned int j, unsigned int k) {
    assertWithMessage(i < n_x && j < n_y && k < n_z, "Grid index out of bounds");
    grid[i*n_y*n_z + j*n_z + k] += val;
}

std::unique_ptr<ScoringGrid> ScoringGrid::convertToDose(const MediumGrid* mediumGrid) {

    std::unique_ptr<ScoringGrid> doseGrid = std::make_unique<ScoringGrid>(level, n_x, n_y, n_z);

    for (size_t i = 0; i < n_x; ++i) {
        for (size_t j = 0; j < n_y; ++j) {
            for (size_t k = 0; k < n_z; ++k) {
                float val = this->getVoxel(i, j, k)/ 
                (mediumGrid->density(i, j, k)*gc::GRID_SPACING_X*gc::GRID_SPACING_Y*gc::GRID_SPACING_Z);
                doseGrid->setVoxel(val, i, j, k);
            }
        }
    }
    return doseGrid;
}


ScoringGrid::~ScoringGrid() {
    delete[] grid;
}

