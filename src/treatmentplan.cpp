#include "treatmentplan.hpp" 




template <typename RandomGen>
TreatmentPlan<RandomGen>::TreatmentPlan(unsigned int level, const MediumGrid* mediumGrid, RandomGen gen):
    level(level), mediumGrid(mediumGrid), gen(gen) {

    scoringGrid = new ScoringGrid(level);
    scoringGrid->initSquareGrid();
    scoringGrid->assignMediumGrid(mediumGrid);
    // // varTracker = new VarianceTracker();
}


template <typename RandomGen>
void TreatmentPlan<RandomGen>::initParticleHistories(unsigned int nPrimaries) {

    particleHistories = new Particle*[nPrimaries];

    float sumShares = 0;
    for (auto& pb : pencilBeams) {
        sumShares += pb->numPrimShare;
    }
    // assert sumShares == 1

    unsigned int planPrimaryIndex = 0;
    for (auto& pb : pencilBeams) {
        unsigned int nPrimariesPencilbeam = nPrimaries * pb->numPrimShare;

        static thread_local std::normal_distribution<double> dist_x(pb->x_0, pb->spread_x);
        static thread_local std::normal_distribution<double> dist_y(pb->y_0, pb->spread_y);
        static thread_local std::normal_distribution<double> dist_z(pb->z_0, pb->spread_z);
        static thread_local std::normal_distribution<double> dist_E(pb->E_0, pb->spread_E);

        double x {}, y {}, z {}, E {}, step {};
        for (unsigned int n = 0; n < nPrimariesPencilbeam; n++) {
            x = dist_x(gen);
            y = dist_y(gen);
            z = dist_z(gen);
            E = dist_E(gen);

            particleHistories[planPrimaryIndex] = new Particle(E, x, y, z, pb->dir_x, pb->dir_y, pb->dir_z, 
                pb->initialStep, 0, planPrimaryIndex, 0, pb->alpha);
            planPrimaryIndex++;
            
            if (planPrimaryIndex >= nPrimaries) {
                break;
            }
        }
    }
}


template <typename RandomGen>
void TreatmentPlan<RandomGen>::initParticleShadowHistories(unsigned int nPrimaries) {

    particleShadowHistories = new ParticleShadow*[nPrimaries];

    for (unsigned int n = 0; n < nPrimaries; n++){
        Particle* p = particleHistories[n];
        particleShadowHistories[n] = new ParticleShadow(p->E, p->index);
    }
}


template <typename RandomGen>
void TreatmentPlan<RandomGen>::addPencilBeam(float numPrimShare, double initialStep,
    char entranceDir, float beamWidth, 
    double dir_x, double dir_y, double dir_z, 
    double x_0, double y_0, double z_0,
    double E, double spreadE, float alpha) {

    PencilBeam<RandomGen>* pb = new PencilBeam<RandomGen>(numPrimShare, initialStep,
        gen, entranceDir, beamWidth, 
        dir_x, dir_y, dir_z,
        x_0, y_0, z_0,
        E, spreadE, alpha);
    pb->initBeamSpread();

    pencilBeams.emplace_back(pb);
    shareParticlesPerBeam.emplace_back(numPrimShare);
    nBeams++;
}


template <typename RandomGen>
TreatmentPlan<RandomGen>::~TreatmentPlan() {
    for (auto& pb : pencilBeams) {
        delete pb;
    }
    delete scoringGrid;
    // delete // varTracker;
}


template <typename RandomGen>
void TreatmentPlan<RandomGen>::simulatePlanWithShadow(unsigned int nPrimaries) { 

    initParticleHistories(nPrimaries);
    initParticleShadowHistories(nPrimaries);
    // varTracker->addSamples(nPrimaries);
    float sumSampleSquares = 0;

    for (unsigned int n = 0; n < nPrimaries; n++) {

        Particle* p = particleHistories[n];
        generateParticleHistory<MediumGrid, RandomGen>(mediumGrid, p, gen);
        ParticleShadow* ps  = particleShadowHistories[n];
        generateParticleShadowHistory<MediumGrid, RandomGen>(mediumGrid, p, ps, gen);

        depositParticleEnergies<ScoringGrid>(scoringGrid, p->next, ps->next);
        deleteParticleHistory(p); // This particle is not refered to among shadow histories
        deleteParticleHistory(ps);
    }
    // for (unsigned int n = 0; n < nPrimaries; n++) {

    //     Particle* p = particleHistories[n];
    //     ParticleShadow* ps = particleShadowHistories[n];
    //     //sumSampleSquaresThread += computeParticleEnergiesVariance(// varTracker, scoringGrid, ps->next, p->next);
    //     sumSampleSquares += computeParticleEnergiesVariance(// varTracker, scoringGrid, p, ps);

    //     deleteParticleHistory(p); // This particle is not refered to among shadow histories
    //     deleteParticleHistory(ps);
    // }
    // varTracker->incrementSum(sumSampleSquares);
    delete[] particleHistories;
    delete[] particleShadowHistories;
}


template <typename RandomGen>
void TreatmentPlan<RandomGen>::simulatePlanWithShadowParallel(unsigned int nPrimaries, unsigned int numThreads) {

    initParticleHistories(nPrimaries);
    initParticleShadowHistories(nPrimaries);
    // varTracker->addSamples(nPrimaries);
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
                generateParticleShadowHistory<MediumGrid, RandomGen>(mediumGrid, p, ps, localGen);

                std::lock_guard<std::mutex> lock(depositMutex);
                depositParticleEnergies<ScoringGrid>(scoringGrid, p->next, ps->next);
                deleteParticleHistory(p); // This particle is not refered to among shadow histories
                deleteParticleHistory(ps);
            };
            // for (unsigned int n = start; n < end; ++n) { // Before computing variance
            //     Particle* p = particleHistories[n];
            //     ParticleShadow* ps  =  particleShadowHistories[n];
            //     // sumSampleSquaresThread += computeParticleEnergiesVariance(// varTracker, scoringGrid, ps->next, p->next);
            //     sumSampleSquaresThread += computeParticleEnergiesVariance(// varTracker, scoringGrid, p, ps);

            //     deleteParticleHistory(p); // This particle is not refered to among shadow histories
            //     deleteParticleHistory(ps);
            // }
            // std::lock_guard<std::mutex> lock(depositMutex);
            // // varTracker->incrementSum(sumSampleSquaresThread);
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
    // varTracker->addSamples(nPrimaries);
    float sumSampleSquares = 0;

    // unsigned int hct = historyChunkSize;
    // unsigned int nHistoryChunks = nPrimaries/historyChunkSize;

    for (unsigned int n = 0; n < nPrimaries; n++) {
        // for (unsigned int n = s*hct; n < (s+1)*hct; n++){
        //     Particle* p = particleHistories[n];
        //     generateParticleHistory<MediumGrid, RandomGen>(mediumGrid, p, gen);
        //     depositParticleEnergy<ScoringGrid, Particle>(scoringGrid, p);
        // }
        Particle* p = particleHistories[n];
        generateParticleHistory<MediumGrid, RandomGen>(mediumGrid, p, gen);

        depositParticleEnergy<ScoringGrid, Particle>(scoringGrid, p);
        deleteParticleHistory(p); 
        // for (unsigned int n = s*hct; n < (s+1)*hct; n++) { // Before computing variance
        //         Particle* p = particleHistories[n];
        //         // sumSampleSquaresThread += computeParticleEnergiesVariance(// varTracker, scoringGrid, ps->next, p->next);
        //         sumSampleSquares += computeParticleEnergiesVariance(// varTracker, scoringGrid, p);

        //         deleteParticleHistory(p); 
        // }
    }

    
    // varTracker->incrementSum(sumSampleSquares);
    delete[] particleHistories;
}


template <typename RandomGen>
void TreatmentPlan<RandomGen>::simulatePlanParallel(unsigned int nPrimaries, unsigned int numThreads) {
    
    initParticleHistories(nPrimaries);
    // varTracker->addSamples(nPrimaries);
    const unsigned int chunkSize = (nPrimaries + numThreads - 1) / numThreads;

    std::mutex depositMutex;
    std::vector<std::thread> threads;

    for (unsigned int t = 0; t < numThreads; ++t) {
        threads.emplace_back([&, t]() {
            RandomGen localGen = gen; // thread-local generator copy

            unsigned int start = t * chunkSize;
            unsigned int end = std::min(start + chunkSize, nPrimaries);
            // float sumSampleSquaresThread = 0;

            for (unsigned int n = start; n < end; n++) {
                Particle* p = particleHistories[n];
                generateParticleHistory<MediumGrid, RandomGen>(mediumGrid, p, localGen);
            
                std::lock_guard<std::mutex> lock(depositMutex);
                depositParticleEnergy<ScoringGrid, Particle>(scoringGrid, p);
                deleteParticleHistory(p); 
    
                // for (unsigned int n = start + t*hct; n < min((t+1)*hct, end); ++n) { // Before computing variance
                //     Particle* p = particleHistories[n];
                //     // sumSampleSquaresThread += computeParticleEnergiesVariance(// varTracker, scoringGrid, ps->next, p->next);
                //     sumSampleSquaresThread += computeParticleEnergiesVariance(// varTracker, scoringGrid, p);

                //     deleteParticleHistory(p); 
                // }
                // std::lock_guard<std::mutex> lock(depositMutex);
                // // varTracker->incrementSum(sumSampleSquaresThread);
            }
        }
    );
    }
    for (auto& thread : threads) {
        thread.join();
    }
    delete[] particleHistories;
}


// VarianceTracker::VarianceTracker() {}


// float VarianceTracker::getVariance() {
//     return sumSampleSquares / nSamples;
// }


// void VarianceTracker::incrementSum(float sampleSquare) {
//     sumSampleSquares += sampleSquare;
// }


// void VarianceTracker::addSamples(unsigned int n) {
//     nSamples += n;
// }



template struct TreatmentPlan<std::mt19937>;
