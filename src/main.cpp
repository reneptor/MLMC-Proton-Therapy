#include "simulation.hpp"
#include <filesystem>


int main(){

    unsigned int seed = 12345678;

    std::string tablePath = (std::filesystem::current_path().parent_path() / "HU2material.csv").string();
    std::string phantomPath = (std::filesystem::current_path().parent_path() / "input" / "phantom.dat").string();

    MLMCprotons<std::mt19937> protonSim(seed, tablePath);
    protonSim.loadPhantom(phantomPath, 100, 100, 500);

    protonSim.addTreatmentPlan(0);

    unsigned int nPrimaries = 10000;
    unsigned int level = 0;
    float initialStep = 0.250050001;
    char entranceDir = 'z';
    float beamWidth = 3.0;
    double dir_x = 0.0;
    double dir_y = 0.0;
    double dir_z = 1.0;
    double x_0 = 50.0;
    double y_0 = 50.0;
    double z_0 = 0.0;
    double E = 230.0;
    double spreadE = 2.3;
    float alpha = 0.25;

    protonSim.addPencilBeam(0, 1.0, initialStep, entranceDir, beamWidth,
    dir_x, dir_y, dir_z, x_0, y_0, z_0, E, spreadE, alpha);
    protonSim.simulateTreatmentPlan(0, 10000);
    protonSim.renderCombinedScoringGrid();

    std::unique_ptr<ScoringGrid> dose = protonSim.yieldDoseCombined();

    std::cout << "Completed\n";
    
    return 0;
}