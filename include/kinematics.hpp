#pragma once

#include "physics_constants.hpp"
#include "voxeltraversal.hpp"
#include "grids.hpp"
#include <array>
#include <random>


template<typename Grid, typename Segs>
float effectiveRSPlengthFromMedium(const Grid* mediumGrid, const Segs& segments);


template<typename Grid, typename Segs>
float effectiveRMSangleFromMedium(const Grid* mediumGrid, const Segs& segments);


template <typename RandomGen>
std::array<double, 3> sampleDir(float MSangle, double dir_x, double dir_y, double dir_z, double E, RandomGen& gen);


template <typename RandomGen>
std::array<double, 3> sampleDirAlt(float MSangle, double dir_x, double dir_y, double dir_z, double E, RandomGen& gen);


template <typename RandomGen>
double sampleEnergyLoss(double RSPlength, double E, RandomGen& gen, double E_real=0, float alpha=0);