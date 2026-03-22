#pragma once

#include "grids.hpp"
#include "pencilbeam.hpp"
#include "particles.hpp"
#include <random>
#include <vector>
#include <random>
#include "error_logging.hpp"


struct VarianceTracker{
    unsigned int nSamples {};
    float sumSampleSquares {};

    VarianceTracker();
    void addSamples(unsigned int n);
    void incrementSum(float sampleSquare);
    float getVariance();
};


template <typename RandomGen>
struct TreatmentPlan {

    unsigned int level;
    MediumGrid* mediumGrid;
    ScoringGrid* scoringGrid;
    VarianceTracker* varTracker;
    unsigned int nBeams {};
    std::vector<PencilBeam<RandomGen>*> pencilBeams;
    std::vector<float> shareParticlesPerBeam;

    Particle** particleHistories {};
    ParticleShadow** particleShadowHistories {};

    TreatmentPlan(unsigned int level, MediumGrid* mediumGrid);
    ~TreatmentPlan();
    void addPencilBeam(float nPrimShare, double initialStep,
    char entranceDir, float beamWidth, double dir_x, double dir_y, double dir_z, 
    double x_0, double y_0, double z_0,
    double E, double spreadE, float alpha);

    void initParticleHistories(unsigned int numPrimaries);
    void initParticleShadowHistories(unsigned int numPrimaries);

    void simulatePlanWithShadow(unsigned int nPrimaries);
    void simulatePlanWithShadowParallel(unsigned int nPrimaries, unsigned int numThreads);
    void simulatePlan(unsigned int nPrimaries);
    void simulatePlanParallel(unsigned int nPrimaries, unsigned int numThreads);
};


