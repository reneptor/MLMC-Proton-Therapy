#pragma once

#include "simulation_constants.hpp"
#include "kinematics.hpp"
#include "voxeltraversal.hpp"
#include "error_logging.hpp"
#include "grids.hpp"


struct Particle{
    double E;                          // MeV - total energy
    double x, y, z;                    // mm
    double dir_x, dir_y, dir_z;        // Normal vector
    double step;                       // mm
    float effectiveRSPlength;           
    unsigned int index {};             // Particle index (for debugging)
    Particle* next = nullptr;          // Reference to next particle in history
    unsigned int particleHistoryCount {};   // Total number of particles in given history

    std::vector<std::tuple<float, unsigned int, unsigned int, unsigned int>> dEhistory;  // Energy history

    Particle(double E, double x, double y, double z,
             double dir_x, double dir_y, double dir_z,
             double step, float effectiveRSPlength, 
             unsigned int index, unsigned int particleHistoryCount);
};


struct ParticleShadow{
    double E;                           // MeV - total energy
    unsigned int index {};              // Particle index (for debugging)
    ParticleShadow* next = nullptr;     // Reference to next particle in coarse history
    std::vector<std::tuple<float, unsigned int, unsigned int, unsigned int>> dEhistory;  // Energy history
    
    ParticleShadow(double E, unsigned int index);
};


template <typename P, typename Segs>
void addEnergyHistory(P* p, double dE, double totalLength, const Segs& segments);


template <typename Grid, typename P>
void depositParticleEnergy(Grid* scoringGrid, P* p);


template <typename P>
void deleteParticleHistory(P* p);


template <typename Grid>
void depositParticleEnergiesMean(Grid* meanScoringGrid, Particle* p, ParticleShadow* ps);


template <typename VarTracker, typename Grid>
float computeParticleEnergiesVariance(const VarTracker* varTracker, const Grid* scoringGrid, Particle* p, ParticleShadow* ps);


template <typename VarTracker, typename Grid>
float computeParticleEnergiesVariance(const VarTracker* varTracker, const Grid* scoringGrid, Particle* p);


template <typename Grid, typename RandomGen>
void generateParticleHistory(const Grid* mediumGrid, Particle* initial, RandomGen& gen);


template <typename Grid, typename RandomGen>
void generateParticleShadowHistory(const Grid* mediumGrid, Particle* initial, ParticleShadow* initialShadow, RandomGen& gen, float alpha);


template <typename Grid, typename RandomGen>
Particle* particleStep(const Grid* mediumGrid, const Particle* p, RandomGen& gen);


template <typename Grid, typename RandomGen>
ParticleShadow* particleShadowStep(const Grid* mediumGrid, const Particle* p_prev, const Particle* p, const Particle* p_next, const ParticleShadow* ps_prev, RandomGen& gen, float alpha);


// template <typename Grid, typename RandomGen>
// ParticleShadow* particleShadowStep(const Grid* mediumGrid, const Particle* p, const ParticleShadow* ps, RandomGen& gen, float alpha);