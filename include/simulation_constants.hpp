#pragma once

#include <array>

// Grid dimension specific constants. Assume these account for all grids

namespace grid_constants {

constexpr unsigned int PHANTOM_N_X = 100;
constexpr unsigned int PHANTOM_N_Y = 100;
constexpr unsigned int PHANTOM_N_Z = 500;
constexpr unsigned int GRID_VOXELS_PER_PHANTOM_X = 1;
constexpr unsigned int GRID_VOXELS_PER_PHANTOM_Y = 1;
constexpr unsigned int GRID_VOXELS_PER_PHANTOM_Z = 1;
constexpr float PHANTOM_SPACING_X = 1.0;                                        // mm
constexpr float PHANTOM_SPACING_Y = 1.0;                                        // mm
constexpr float PHANTOM_SPACING_Z = 1.0;                                        // mm
constexpr unsigned int GRID_N_X = PHANTOM_N_X*GRID_VOXELS_PER_PHANTOM_X;
constexpr unsigned int GRID_N_Y = PHANTOM_N_Y*GRID_VOXELS_PER_PHANTOM_Y;
constexpr unsigned int GRID_N_Z = PHANTOM_N_Z*GRID_VOXELS_PER_PHANTOM_Z;
constexpr float GRID_SPACING_X = PHANTOM_SPACING_X/GRID_VOXELS_PER_PHANTOM_X;   // mm
constexpr float GRID_SPACING_Y = PHANTOM_SPACING_Y/GRID_VOXELS_PER_PHANTOM_Y;   // mm
constexpr float GRID_SPACING_Z = PHANTOM_SPACING_Z/GRID_VOXELS_PER_PHANTOM_Z;   // mm
constexpr float PHANTOM_LENGTH_X = PHANTOM_N_X*PHANTOM_SPACING_X;               // mm
constexpr float PHANTOM_LENGTH_Y = PHANTOM_N_Y*PHANTOM_SPACING_Y;               // mm
constexpr float PHANTOM_LENGTH_Z = PHANTOM_N_Z*PHANTOM_SPACING_Z;               // mm



inline std::array<unsigned int, 3> gridIndices(double x, double y, double z) {
    return {
        static_cast<unsigned int>(x / GRID_SPACING_X),
        static_cast<unsigned int>(y / GRID_SPACING_Y),
        static_cast<unsigned int>(z / GRID_SPACING_Z)
    };
}
}


namespace default_constants {

constexpr float BEAM_ENERGY = 100.0;        // MeV
constexpr float BEAM_ENERGY_SPREAD = 1.0;   // MeV (1% of energy)
constexpr float BEAM_WIDTH = 1.0;          // mm
constexpr char  ENTRANCE_DIR = 'z';        // 'x', 'y', or 'z' - direction of beam entrance
constexpr float ENTRY_X = 50.0;            // mm
constexpr float ENTRY_Y = 50.0;            // mm            
constexpr float ENTRY_Z = 0.0;              // mm
constexpr float DIR_X = 0.0;
constexpr float DIR_Y = 0.0;
constexpr float DIR_Z = 1.0;
constexpr float ALPHA = 1.0;

}


namespace control_constants {

constexpr float MINIMUM_THRESHOLD_ENERGY = 1.0; // Remaining particle energy is deposited in local voxel
// constexpr unsigned int SEED = 12345678;
constexpr float ENERGY_DEPOSIT_THRESHOLD = 10.0; // MeV, for catching numerical errors

}
