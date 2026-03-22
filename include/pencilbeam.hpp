#pragma once


#include "grids.hpp"
#include "particles.hpp"
#include "error_logging.hpp"
#include "simulation_constants.hpp"
#include <random>
#include <mutex>
#include <thread>


namespace dc = default_constants;


template <typename RandomGen>
struct PencilBeam {

    float numPrimShare;
    double initialStep;
    double dir_x, dir_y, dir_z;
    double x_0, y_0, z_0;
    double spread_x, spread_y, spread_z;
    double E_0;
    double spread_E;
    float alpha;
    RandomGen& gen;

    float beamWidth {};
    char entranceDir {};

    PencilBeam(float numPrimShare, double initialStep, RandomGen& gen, 
               char entranceDir, float beamWidth, 
               double dir_x, double dir_y, double dir_z,
               double x_0, double y_0, double z_0,
               double E, double spread_E, float alpha);

    bool validEntrance(double x_0, double y_0, double z_0);
    bool validEntranceDir(char entranceDir);

    void initBeamSpread();
};




