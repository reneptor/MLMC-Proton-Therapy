#include "particles.hpp"


namespace cc = control_constants;
namespace gc = grid_constants;


// template void generateParticleHistory<MediumGrid, std::mt19937>(MediumGrid*, Particle*, std::mt19937&);
// template void generateParticleShadowHistory<MediumGrid, std::mt19937>(MediumGrid*, Particle*, ParticleShadow*, std::mt19937&);

// template void generateParticleHistory<MediumGrid, std::mt19937>(const MediumGrid*, Particle*, std::mt19937&);
// template void generateParticleShadowHistory<MediumGrid, std::mt19937>(const MediumGrid*, Particle*, ParticleShadow*, std::mt19937&, float);
// template void depositParticleEnergy<ScoringGrid, Particle>(ScoringGrid*, Particle*);
// template void depositParticleEnergy<ScoringGrid, ParticleShadow>(ScoringGrid*, ParticleShadow*);
// template void deleteParticleHistory<Particle>(Particle*);
// template void deleteParticleHistory<ParticleShadow>(ParticleShadow*);


Particle::Particle(double E, double x, double y, double z,
            double dir_x, double dir_y, double dir_z,
            double step,
            float effectiveRSPlength,
            unsigned int index, 
            unsigned int particleHistoryCount)
    : E(E), 
        x(x), y(y), z(z),
        dir_x(dir_x), dir_y(dir_y), dir_z(dir_z),
        step(step), effectiveRSPlength(effectiveRSPlength),
        index(index), particleHistoryCount(particleHistoryCount) {}


ParticleShadow::ParticleShadow(double E, unsigned int index) 
    : E(E), index(index) {}


template <typename P, typename Segs>
void addEnergyHistory(P* p, double dE, double totalLength, const Segs& segments){
    for (auto s : segments){
        p->dEhistory.push_back(std::make_tuple(dE*s.segmentLength()/totalLength, s.i, s.j, s.k));
    }
}


template <typename Grid, typename P>
void depositParticleEnergy(Grid* scoringGrid, P* p){
    while (p){
        for (auto [dE, i, j, k] : p->dEhistory){
            scoringGrid->depositEnergyLocal(i, j, k, dE);
        }
        p = p->next;
    }
}


template <typename P>
void deleteParticleHistory(P* p){
    while (p){
        P* p_old = p;
        p = p->next;
        delete p_old;
    }
}


template <typename Grid>
void depositParticleEnergiesMean(Grid* scoringGrid, Particle* p, ParticleShadow* ps){
    while (ps){

        unsigned int firstIndex = p->dEhistory.size();
        unsigned int secondIndex = p->next->dEhistory.size();

        for (int index = 0; index < firstIndex; index++){
            float dE = std::get<0>(p->dEhistory[index]); 
            const auto [dEs, i, j, k] = ps->dEhistory[index];
            scoringGrid->depositEnergyLocal(i, j, k, dE-dEs);
        }
        for (int index = firstIndex; index < firstIndex+secondIndex; index++){
            float dE = std::get<0>(p->next->dEhistory[index - firstIndex]);
            const auto [dEs, i, j, k] = ps->dEhistory[index];
            scoringGrid->depositEnergyLocal(i, j, k, dE-dEs);
        }
        ps = ps->next;
        p = p->next->next;
    }
}


template <typename VarTracker, typename Grid>
float computeParticleEnergiesVariance(const VarTracker* varTracker, const Grid* scoringGrid, Particle* p, ParticleShadow* ps){

    float sumSampleSquares = 0;

    while (ps){

        float sampleSquare = 0; // Every (particle + shadow) counts as a sample.

        unsigned int firstIndex = p->dEhistory.size();
        unsigned int secondIndex = p->next->dEhistory.size();

        for (int index = 0; index < firstIndex; index++){
            float dE = std::get<0>(p->dEhistory[index]); 
            const auto [dEs, i, j, k] = ps->dEhistory[index];

            float dE_diff_mean = scoringGrid->getEnergyLocal(i, j, k)/varTracker->nSamples;
            float dE_diff_sample = dE - dEs;
            sampleSquare += (dE_diff_sample - dE_diff_sample)*(dE_diff_sample - dE_diff_sample);
        }
        for (int index = firstIndex; index < firstIndex+secondIndex; index++){
            float dE = std::get<0>(p->next->dEhistory[index - firstIndex]);
            const auto [dEs, i, j, k] = ps->dEhistory[index];
            
            float dE_diff_mean = scoringGrid->getEnergyLocal(i, j, k)/varTracker->nSamples;
            float dE_diff_sample = dE - dEs;
            sampleSquare += (dE_diff_sample - dE_diff_sample)*(dE_diff_sample - dE_diff_sample);
        }
        sumSampleSquares += sampleSquare;
        ps = ps->next;
        p = p->next->next;
    }
    return sumSampleSquares;
}


template <typename VarTracker, typename Grid>
float computeParticleEnergiesVariance(const VarTracker* varTracker, const Grid* scoringGrid, Particle* p){

    float sumSampleSquares = 0;

    while (p){

        float sampleSquare = 0; // Every (particle + shadow) counts as a sample.

        unsigned int maxIndex = p->dEhistory.size();

        for (int index = 0; index < maxIndex; index++){
            const auto [dE, i, j, k] = p->dEhistory[index];

            float dE_diff_mean = scoringGrid->getEnergyLocal(i, j, k)/varTracker->nSamples;
            float dE_diff_sample = dE;
            sampleSquare += (dE_diff_sample - dE_diff_sample)*(dE_diff_sample - dE_diff_sample);
        }
        p = p->next;
    }
    return sumSampleSquares;
}




template <typename Grid, typename RandomGen>
void generateParticleHistory(const Grid* mediumGrid, Particle* initial, RandomGen& gen){
    Particle* p = initial;  // This particle responsible for complete history
    while (p){
        if (p->E < cc::MINIMUM_THRESHOLD_ENERGY && p->particleHistoryCount%2 == 0){
            p->next = nullptr;
            break; // Particle energy below threshold, stop history
        }
        Particle* p_next = particleStep<Grid, RandomGen>(mediumGrid, p, gen);
        p->next = p_next;
        p = p_next; 
    }
}


template <typename Grid, typename RandomGen>
void generateParticleShadowHistory(const Grid* mediumGrid, Particle* initial, ParticleShadow* initialShadow, RandomGen& gen, float alpha){
    Particle* p = initial;
    ParticleShadow* ps = initialShadow; // This particle responsible for complete history
    while (p->next && p->next->next){

        ParticleShadow* ps_next = particleShadowStep<Grid, RandomGen>(mediumGrid, p, p->next, p->next->next, ps, gen, alpha);
        ps->next = ps_next;
        ps = ps_next;
        p = p->next->next;
    }
}


template <typename Grid, typename RandomGen>
Particle* particleStep(const Grid* mediumGrid, const Particle* p, RandomGen& gen){

    double x = p->x;
    double y = p->y;
    double z = p->z;
    double x_next = x + p->dir_x * p->step;
    double y_next = y + p->dir_y * p->step;
    double z_next = z + p->dir_z * p->step;

    if (!((0 <= x_next && x_next <= gc::PHANTOM_LENGTH_X) && (0 <= y_next && y_next <= gc::PHANTOM_LENGTH_Y) && (0 <= z_next && z_next <= gc::PHANTOM_LENGTH_Z))){
        return nullptr; // Ended up outside phantom boundaries, no error ofc
    } 

    GridPath gridPath(x, y, z, x_next, y_next, z_next);
    gridPath.decomposePath();
    std::vector<LineSegment> segments = gridPath.traversalSegments;

    float effectiveRSPlength = effectiveRSPlengthFromMedium<Grid, std::vector<LineSegment>>(mediumGrid, segments);
    float effectiveRMSangle = effectiveRMSangleFromMedium<Grid, std::vector<LineSegment>>(mediumGrid, segments);

    double dE = sampleEnergyLoss<RandomGen>(effectiveRSPlength, p->E, gen);
    auto [dir_x, dir_y, dir_z] = sampleDir<RandomGen>(effectiveRMSangle, p->dir_x, p->dir_y, p->dir_z, p->E, gen);    
    // auto [dir_x, dir_y, dir_z] = sampleDirAlt<RandomGen>(effectiveRMSangle, p->dir_x, p->dir_y, p->dir_z, p->E, gen);    


    double E_next = p->E - dE;
    double step_next = p->step * E_next/p->E;

    Particle* p_next = new Particle(
        E_next,
        x_next, y_next, z_next,
        dir_x, dir_y, dir_z,
        step_next,              // 
        effectiveRSPlength,     // All used for shadow particle calculations  
        p->index,
        p->particleHistoryCount+1
    );

    addEnergyHistory<Particle, std::vector<LineSegment>>(p_next, dE, step_next, segments);
    return p_next;
}


template <typename Grid, typename RandomGen>
ParticleShadow* particleShadowStep(const Grid* mediumGrid, const Particle* p_prev, Particle* p, Particle* p_next, const ParticleShadow* ps_prev, RandomGen& gen,
    float alpha){

    float effectiveRSPlengthTotal = p->effectiveRSPlength + p_next->effectiveRSPlength;
    double dE = sampleEnergyLoss<RandomGen>(effectiveRSPlengthTotal, ps_prev->E, gen, p_prev->E, alpha);
    double E = ps_prev->E - dE;

    ParticleShadow* ps = new ParticleShadow(E, ps_prev->index);

    // Recaulculate these because can not save segments in memory with current implementation
    GridPath gridPathFirst(p_prev->x, p_prev->y, p_prev->z, p->x, p->y, p->z);
    gridPathFirst.decomposePath();
    std::vector<LineSegment> segmentsFirstHalf = gridPathFirst.traversalSegments;

    GridPath gridPathSecond(p->x, p->y, p->z, p_next->x, p_next->y, p_next->z);
    gridPathSecond.decomposePath();
    std::vector<LineSegment> segmentsSecondHalf = gridPathSecond.traversalSegments;

    float dEfirstHalf = dE * p->step / (p->step + p_next->step);
    float dEsecondHalf = dE * p_next->step / (p->step + p_next->step);

    addEnergyHistory<ParticleShadow, std::vector<LineSegment>>(ps, dEfirstHalf, p->step, segmentsFirstHalf);
    addEnergyHistory<ParticleShadow, std::vector<LineSegment>>(ps, dEsecondHalf, p_next->step, segmentsSecondHalf);

    return ps;
}

// template <typename Grid, typename RandomGen>
// ParticleShadow* particleShadowStep(const Grid* mediumGrid, const Particle* p, const ParticleShadow* ps, RandomGen& gen,
//     float alpha){
//     // ps corresponds to shadow of p (same step)

//     if (!p->next){
//         return nullptr; 
//     }

//     double x_start = p->x;
//     double y_start = p->y;
//     double z_start = p->z;

//     Particle* pFirst = p->next;
//     double x_midd = pFirst->x;
//     double y_midd = pFirst->y;
//     double z_midd = pFirst->z;

//     double firstStep = std::sqrt( // Necessary for energy history
//         (x_midd - x_start)*(x_midd - x_start) +
//         (y_midd - y_start)*(y_midd - y_start) +
//         (z_midd - z_start)*(z_midd - z_start)
//     );
//     double secondStep;

//     Particle* pSecond = pFirst->next;
//     double x_end = pSecond->x;
//     double y_end = pSecond->y;
//     double z_end = pSecond->z;

//     double secondStep = std::sqrt( // Necessary for energy history
//         (x_end - x_midd)*(x_end - x_midd) +
//         (y_end - y_midd)*(y_end - y_midd) +
//         (z_end - z_midd)*(z_end - z_midd)
//     );

//     GridPath gridPathFirst(x_start, y_start, z_start, x_midd, y_midd, z_midd);
//     gridPathFirst.decomposePath();
//     std::vector<LineSegment> segmentsFirst = gridPathFirst.traversalSegments;
//     std::vector<LineSegment> segmentsSecond;
//     float effectiveRSPlengthFirst = effectiveRSPlengthFromMedium<Grid, std::vector<LineSegment>>(mediumGrid, segmentsFirst);
//     float effectiveRSPlengthSecond = 0;

    
//     if (pSecond){
//         pSecond = pFirst->next; // MLMC principle for coarse layers - two fine for every coarse
        

//         GridPath gridPathSecond(x_midd, y_midd, z_midd, x_end, y_end, z_end);
//         gridPathSecond.decomposePath();
//         segmentsSecond = gridPathSecond.traversalSegments;
//         effectiveRSPlengthSecond = effectiveRSPlengthFromMedium<Grid, std::vector<LineSegment>>(mediumGrid, segmentsSecond);
//     }

//     float effectiveRSPlengthTotal = effectiveRSPlengthFirst + effectiveRSPlengthSecond;

//     double dE = sampleEnergyLoss<RandomGen>(effectiveRSPlengthTotal, ps->E, gen, p->E, alpha);

//     double E_next = ps->E - dE;

//     ParticleShadow* ps_next = new ParticleShadow(
//         E_next,
//         dE,
//         pFirst,
//         pSecond
//     );

//     float dEfirst = dE * firstStep / (firstStep + secondStep);
//     float dEsecond = dE * secondStep / (firstStep + secondStep);

//     addEnergyHistory<ParticleShadow, std::vector<LineSegment>>(ps_next, dEfirst, firstStep, segmentsFirst);
//     if (dEsecond != 0){
//         addEnergyHistory<ParticleShadow, std::vector<LineSegment>>(ps_next, dEsecond, secondStep, segmentsSecond);
//     }
    
//     return ps_next;
// }