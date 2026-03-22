#include "simulation.hpp"


template struct MLMCprotons<std::mt19937>;



template <typename RandomGen>
const unsigned int MLMCprotons<RandomGen>::maxThreads = std::thread::hardware_concurrency();


template <typename RandomGen>
MLMCprotons<RandomGen>::MLMCprotons(unsigned int seed, const std::string tablePath) : hu2mat(tablePath), gen(seed) {
    mediumGrid = new MediumGrid(&hu2mat);
    combinedScoringGrid = new ScoringGrid(0);
}


template <typename RandomGen>
MLMCprotons<RandomGen>::~MLMCprotons(){
    // for (const auto& s : scoringGrids){
    //     delete s;
    // }
    for (const auto& [_, tp] : treatmentPlans) {
        delete tp;
    }
    delete mediumGrid;
    delete combinedScoringGrid;
}


template <typename RandomGen>
void MLMCprotons<RandomGen>::renderCombinedScoringGrid() {

    for (unsigned int i = 0; i < gc::GRID_N_X; ++i) {
        for (unsigned int j = 0; j < gc::GRID_N_Y; ++j) {
            for (unsigned int k = 0; k < gc::GRID_N_Z; ++k) {

                combinedScoringGrid->setVoxel(0, i, j, k);
                for (auto& [_, tp] : treatmentPlans) {
                    combinedScoringGrid->addVoxel(tp->scoringGrid->getVoxel(i, j, k), i, j, k);
                }
            }
        }
    }
}


template <typename RandomGen>
void MLMCprotons<RandomGen>::loadPhantom(const std::string phantomPath, size_t n_x, size_t n_y, size_t n_z){

    assertWithMessage((gc::GRID_N_X%n_x == 0 && gc::GRID_N_Y%n_y == 0 && gc::GRID_N_Z%n_z == 0),
    "Missmatch between grid and phantom dimensions");

    std::filesystem::path fullPath = std::filesystem::current_path() / "input" / phantomPath;

    std::ifstream in(fullPath, std::ios::binary);
    assertWithMessage(in.is_open(), "Could not open file: " + phantomPath);

    int16_t* phantom = new int16_t[n_x*n_y*n_z];
    in.read(reinterpret_cast<char*>(phantom), n_x * n_y * n_z * sizeof(int16_t));

    for (unsigned int i = 0; i < gc::GRID_N_X; ++i) {
        for (unsigned int j = 0; j < gc::GRID_N_Y; ++j) {
            for (unsigned int k = 0; k < gc::GRID_N_Z; ++k) {

                unsigned int phantom_i = i / gc::GRID_VOXELS_PER_PHANTOM_X;
                unsigned int phantom_j = j / gc::GRID_VOXELS_PER_PHANTOM_Y;
                unsigned int phantom_k = k / gc::GRID_VOXELS_PER_PHANTOM_Z;

                mediumGrid->grid[i*gc::GRID_N_Y*gc::GRID_N_Z + j*gc::GRID_N_Z + k] = phantom[phantom_i * n_y * n_z + phantom_j * n_z + phantom_k];
            }
        }
    }
    delete[] phantom;
}


template <typename RandomGen>
void MLMCprotons<RandomGen>::addTreatmentPlan(unsigned int level){
    treatmentPlans[level] = new TreatmentPlan<RandomGen>(level, mediumGrid);
}


template <typename RandomGen>
void MLMCprotons<RandomGen>::addPencilBeam(unsigned int level, float nPrimShare, double initialStep, 
               char entranceDir, float beamWidth, 
               double dir_x, double dir_y, double dir_z,
               double x_0, double y_0, double z_0,
               double E, double spread_E, float alpha){
    treatmentPlans[level]->addPencilBeam(nPrimShare, initialStep, 
               entranceDir, beamWidth, 
               dir_x, dir_y, dir_z,
               x_0, y_0, z_0,
               E, spread_E, alpha);
}


template <typename RandomGen>
void MLMCprotons<RandomGen>::simulateTreatmentPlan(unsigned int level, unsigned int numPrimaries,
        unsigned int numThreads) {
    if (level == 0) {
        if (numThreads == 1) {
            treatmentPlans[0]->simulatePlan(numPrimaries);
        }
        else {
            treatmentPlans[0]->simulatePlanParallel(numPrimaries, numThreads);
        }
    }
    else {
        if (numThreads == 1) {
            treatmentPlans[level]->simulatePlanWithShadow(numPrimaries);
        }
        else {
            treatmentPlans[level]->simulatePlanWithShadowParallel(numPrimaries, numThreads);
        }
    }
}


template <typename RandomGen>
std::unique_ptr<ScoringGrid> MLMCprotons<RandomGen>::yieldDoseAtLevel(unsigned int level) {
    // assertWithMessage(index < scoringGrids.size(), "Grid index out of bounds");
    std::unique_ptr<ScoringGrid> doseAtLevel = treatmentPlans[level]->scoringGrid->convertToDose(mediumGrid);
    return doseAtLevel;
}


template <typename RandomGen>
std::unique_ptr<ScoringGrid> MLMCprotons<RandomGen>::yieldDoseCombined() {
    std::unique_ptr<ScoringGrid> doseCombined = combinedScoringGrid->convertToDose(mediumGrid);
    return doseCombined;
}




