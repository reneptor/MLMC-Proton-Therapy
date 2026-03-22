#include "treatmentplan.hpp" 



template struct TreatmentPlan<std::mt19937>;


template <typename RandomGen>
TreatmentPlan<RandomGen>::TreatmentPlan(unsigned int level, MediumGrid* mediumGrid):
    level(level), MediumGrid(mediumGrid) 


void TreatmentPlan::initParticleHistories(unsigned int nPrimaries) {

    particleHistories = new Particle*[nPrimaries];

    float sumShares = 0;
    for (auto& pb : pencilBeams) {
        sumShares += pb->nPrimShare;
    }
    // assert sumShares == 1

    unsigned int planPrimaryIndex = 0;
    for (auto& pb : pencilBeams) {
        nPrimariesPencilbeam = nPrimaries * pb->nPrimShare;

        static thread_local std::normal_distribution<double> dist_x(pb->x_0, pb->spread_x);
        static thread_local std::normal_distribution<double> dist_y(pb->y_0, pb->spread_y);
        static thread_local std::normal_distribution<double> dist_z(pb->z_0, pb->spread_z);
        static thread_local std::normal_distribution<double> dist_E(pb->E_0, pb->spread_E);

        double x {}, y {}, z {}, E {}, step {};
        for (int n = 0; n < nPrimariesPencilbeam; n++) {
            x = dist_x(gen);
            y = dist_y(gen);
            z = dist_z(gen);
            E = dist_E(gen);

            particleHistories[planPrimaryIndex] = new Particle(E, x, y, z, pb->dir_x, pb->dir_y, pb->dir_z, pb->initialStep, 0, planPrimaryIndex, 0);
            planPrimaryIndex += 1;
        }
    }
}


template <typename RandomGen>
void TreatmentPlan<RandomGen>::initParticleShadowHistories(unsigned int nPrimaries) {

    particleShadowHistories = new ParticleShadow*[nPrimaries];

    for (int n = 0; n < nPrimaries; n++){
        Particle* p = particleHistories[n];
        particleShadowHistories[n] = new ParticleShadow(p->E, p->index);
    }
}


template <typename RandomGen>
void TreatmentPlan<RandomGen>::addPencilBeam(float nPrimShare, double initialStep,
    char entranceDir, float beamWidth, 
    double dir_x, double dir_y, double dir_z, 
    double x_0, double y_0, double z_0,
    double E, double spreadE, float alpha) {

    PencilBeam<RandomGen>* pb = new PencilBeam<RandomGen>(nPrimShare, initialStep,
        gen, entranceDir, beamWidth, 
        dir_x, dir_y, dir_z,
        x_0, y_0, z_0,
        E, spread_E, alpha);
    pb->initBeamSpread(
    entranceDir, beamWidth, dir_x, dir_y, dir_z,
    x_0, y_0, z_0,
    E, spreadE, 
    alpha, gen);

    pencilBeams.emplace_back(pb);
    shareParticlesPerBeam.emplace_back(share);
    nBeams += 1;
}


template <typename RandomGen>
TreatmentPlan<RandomGen>::~TreatmentPlan() {
    for (auto& pb : pencilBeams) {
        delete pb
    }
    delete scoringGrid;
}


template <typename RandomGen>
void TreatmentPlan<RandomGen>::simulatePlanWithShadow(unsigned int nPrimaries) { 

    initParticleHistories(nPrimaries);
    initParticleShadowHistories(nPrimaries);
    varTracker.addSamples(nPrimaries);
    float sumSampleSquares = 0;

    for (unsigned int n = 0; n < nPrimaries; n++) {

        Particle* p = particleHistories[n];
        generateParticleHistory<MediumGrid, RandomGen>(mediumGrid, p, gen);
        ParticleShadow* ps  = particleShadowHistories[n];
        generateParticleShadowHistory<MediumGrid, RandomGen>(mediumGrid, p, ps, gen, alpha);
        
        depositParticleEnergiesMean<ScoringGrid>(scoringGrid, ps->next, p->next);
    }
    for (unsigned int n = 0; n < nPrimaries; n++) {

        Particle* p = particleHistories[n];
        ParticleShadow* ps = particleShadowHistories[n];
        //sumSampleSquaresThread += computeParticleEnergiesVariance(varTracker, scoringGrid, ps->next, p->next);
        sumSampleSquares += computeParticleEnergiesVariance(varTracker, scoringGrid, p, ps);

        deleteParticleHistory(p); // This particle is not refered to among shadow histories
        deleteParticleHistory(ps);
    }
    varTracker.incrementSum(sumSampleSquares);
    delete[] particleHistories;
    delete[] particleShadowHistories;
}


template <typename RandomGen>
void TreatmentPlan<RandomGen>::simulatePlanWithShadowParallel(unsigned int nPrimaries, unsigned int numThreads) {

    initParticleHistories(nPrimaries);
    initParticleShadowHistories(nPrimaries);
    varTracker.addSamples(nPrimaries);
    const unsigned int chunkSize = (nPrimaries + numThreads - 1) / numThreads;

    std::mutex depositMutex;
    std::vector<std::thread> threads;

    for (unsigned int t = 0; t < numThreads; ++t) {
        threads.emplace_back([&, t]() {
            RandomGen localGen = gen; // Thread-local generator copy
            unsigned int start = t * chunkSize;
            unsigned int end = std::min(start + chunkSize, nPrimaries);
            float sumSampleSquaresThread = 0;

            for (unsigned int n = start; n < end; ++n) { // Must compute mean
                Particle* p = particleHistories[n];
                generateParticleHistory<MediumGrid, RandomGen>(mediumGrid, p, localGen);
                ParticleShadow* ps  =  particleShadowHistories[n];
                generateParticleShadowHistory<MediumGrid, RandomGen>(mediumGrid, p, ps, localGen, alpha);

                std::lock_guard<std::mutex> lock(depositMutex);
                depositParticleEnergiesMean<ScoringGrid>(scoringGrid, ps->next, p->next);
            };
            for (unsigned int n = start; n < end; ++n) { // Before computing variance
                Particle* p = particleHistories[n];
                ParticleShadow* ps  =  particleShadowHistories[n];
                // sumSampleSquaresThread += computeParticleEnergiesVariance(varTracker, scoringGrid, ps->next, p->next);
                sumSampleSquaresThread += computeParticleEnergiesVariance(varTracker, scoringGrid, p, ps);

                deleteParticleHistory(p); // This particle is not refered to among shadow histories
                deleteParticleHistory(ps);
            }
            std::lock_guard<std::mutex> lock(depositMutex);
            varTracker.incrementSum(sumSampleSquaresThread);
        }    
    );
    }
    for (auto& thread : threads) {
        thread.join();
    }
    delete[] particleHistories;
    delete[] particleShadowHistories;
}


template <typename RandomGen>
void TreatmentPlan<RandomGen>::simulatePlan(unsigned int nPrimaries) { 
    initParticleHistories(nPrimaries);
    varTracker.addSamples(nPrimaries);
    float sumSampleSquares = 0;

    for (unsigned int n = 0; n < nPrimaries; n++){

        Particle* p = particleHistories[n];
        generateParticleHistory<MediumGrid, RandomGen>(mediumGrid, p, gen);
        depositParticleEnergy<ScoringGrid, Particle>(scoringGrid, p);
    }
    for (unsigned int n = start; n < end; ++n) { // Before computing variance
            Particle* p = particleHistories[n];
            ParticleShadow* ps  =  particleShadowHistories[n];
            // sumSampleSquaresThread += computeParticleEnergiesVariance(varTracker, scoringGrid, ps->next, p->next);
            sumSampleSquares += computeParticleEnergiesVariance(varTracker, scoringGrid, p);

            deleteParticleHistory(p); 
    }
    varTracker.incrementSum(sumSampleSquares);
    delete[] particleHistories;
}


template <typename RandomGen>
void TreatmentPlan<RandomGen>::simulatePlanParallel(unsigned int nPrimaries, unsigned int numThreads) {
    
    initParticleHistories(nPrimaries);
    varTracker.addSamples(nPrimaries);
    const unsigned int chunkSize = (nPrimaries + numThreads - 1) / numThreads;

    std::mutex depositMutex;
    std::vector<std::thread> threads;

    for (unsigned int t = 0; t < numThreads; ++t) {
        threads.emplace_back([&, t]() {
            RandomGen localGen = gen; // thread-local generator copy

            unsigned int start = t * chunkSize;
            unsigned int end = std::min(start + chunkSize, nPrimaries);
            float sumSampleSquaresThread = 0;

            for (unsigned int n = start; n < end; ++n) {
                Particle* p = particleHistories[n];
                generateParticleHistory<MediumGrid, RandomGen>(mediumGrid, p, localGen);
                {
                    std::lock_guard<std::mutex> lock(depositMutex);
                    depositParticleEnergy<ScoringGrid, Particle>(scoringGrid, p);
                }
            }
            for (unsigned int n = start; n < end; ++n) { // Before computing variance
                Particle* p = particleHistories[n];
                // sumSampleSquaresThread += computeParticleEnergiesVariance(varTracker, scoringGrid, ps->next, p->next);
                sumSampleSquaresThread += computeParticleEnergiesVariance(varTracker, scoringGrid, p);

                deleteParticleHistory(p); 
            }
            std::lock_guard<std::mutex> lock(depositMutex);
            varTracker.incrementSum(sumSampleSquaresThread);
        }
    );
    }
    for (auto& thread : threads) {
        thread.join();
    }
    delete[] particleHistories;
}


VarianceTracker::VarianceTracker() {}


float VarianceTracker::getVariance() {
    return sumSampleSquares / nSamples;
}


void VarianceTracker::incrementSum(float sampleSquare) {
    sumSampleSquares += sampleSquare;
}


void VarianceTracker::addSamples(unsigned int n) {
    nSamples += n;
}
