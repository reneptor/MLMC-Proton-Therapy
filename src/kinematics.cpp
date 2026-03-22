#include "kinematics.hpp"


namespace pc = physics_constants;


template float effectiveRSPlengthFromMedium<MediumGrid, std::vector<LineSegment>>(const MediumGrid*, const std::vector<LineSegment>&);
template float effectiveRMSangleFromMedium<MediumGrid, std::vector<LineSegment>>(const MediumGrid*, const std::vector<LineSegment>&);
// template double sampleEnergyLoss<std::mt19937>(double, double, std::mt19937&, double, float);
// template std::array<double, 3> sampleDir<std::mt19937>(float, double, double, double, double, std::mt19937&);


template<typename Grid, typename Segs>
float effectiveRSPlengthFromMedium(const Grid* mediumGrid, const Segs& segments){
    float RSPlength = 0;
    for (auto s : segments){
        float RSP = mediumGrid->RSP(s.i, s.j, s.k);
        RSPlength += RSP*s.segmentLength();
    }
    return RSPlength;
}


template<typename Grid, typename Segs>
float effectiveRMSangleFromMedium(const Grid* mediumGrid, const Segs& segments){
    float MSangle = 0;
    for (auto s : segments){
        float X = mediumGrid->X(s.i, s.j, s.k);
        MSangle += pc::MSAngleHighland(s.segmentLength()/X);
    }
    return MSangle;
}


template <typename RandomGen>
std::array<double, 3> sampleDir(float MSangle, double dir_x, double dir_y, double dir_z, double E, RandomGen& gen){
    float angleEnergyDependance = pc::HighlandEnergyDependance(E);
    float RMSangle = 13.6*angleEnergyDependance*std::sqrt(MSangle);
    static thread_local std::normal_distribution<float> theta_normal(0.0, RMSangle);
    static thread_local std::normal_distribution<double> unit_normal(0.0, 1.0);

    float theta = theta_normal(gen);
    double s_x = unit_normal(gen);
    double s_y = unit_normal(gen);
    double s_z = unit_normal(gen);

    // Generate random unit vector and normalize
    double s_norm = std::sqrt(s_x*s_x + s_y*s_y + s_z*s_z);
    s_x /= s_norm;
    s_y /= s_norm;
    s_z /= s_norm;

    // Cross product to generate random perpendicular vector
    double k_x = dir_y*s_z - dir_z*s_y;
    double k_y = dir_z*s_x - dir_x*s_z;
    double k_z = dir_x*s_y - dir_y*s_x;

    // Normalize
    double k_norm = std::sqrt(k_x*k_x + k_y*k_y + k_z*k_z);
    k_x /= k_norm;
    k_y /= k_norm;
    k_z /= k_norm;

    // Precompute sin and cos
    // auto [sin_theta, cos_theta] = std::sincos(theta);
    float sin_theta = std::sin(theta);
    float cos_theta = std::cos(theta);

    // Apply Rodrigues formula for rotation
    double new_dir_x = dir_x*cos_theta + (k_y*dir_z - k_z*dir_y)*sin_theta;
    double new_dir_y = dir_y*cos_theta + (k_z*dir_x - k_x*dir_z)*sin_theta;
    double new_dir_z = dir_z*cos_theta + (k_x*dir_y - k_y*dir_x)*sin_theta;

    // Normalize the new direction vector (necessary?)
    // double new_dir_norm = std::sqrt(new_dir_x*new_dir_x + new_dir_y*new_dir_y + new_dir_z*new_dir_z);
    // new_dir_x /= new_dir_norm;
    // new_dir_y /= new_dir_norm;
    // new_dir_z /= new_dir_norm;

    std::array<double, 3> newDir = {
        new_dir_x,
        new_dir_y,
        new_dir_z 
    };

    return newDir;
}


template <typename RandomGen>
std::array<double, 3> sampleDirAlt(float MSangle,
                                double dir_x, double dir_y, double dir_z,
                                double E,
                                RandomGen& gen)
{
    // --- 1. Compute RMS scattering angle ---
    float angleEnergyDependance = pc::HighlandEnergyDependance(E);
    float RMSangle = 13.6f * angleEnergyDependance * std::sqrt(MSangle);

    // --- 2. thread_local RNG distributions ---
    static thread_local std::normal_distribution<float> theta_normal(0.0, RMSangle);
    static thread_local std::uniform_real_distribution<double> phi_uniform(0.0, 2.0 * M_PI);

    // --- 3. Sample scattering angle (Gaussian) and azimuth (uniform) ---
    double theta = theta_normal(gen);
    double phi   = uniform_phi(gen);

    // auto [sin_theta, cos_theta] = std::sincos(theta);
    // auto [sin_phi, cos_phi] = std::sincos(phi);
    float sin_theta = std::sin(theta);
    float cos_theta = std::cos(theta);
    float sin_phi = std::sin(phi);
    float cos_phi = std::cos(phi);

    // --- 4. Build an orthonormal basis (e1, e2, dir) ---

    // Choose a vector not parallel to (dir_x, dir_y, dir_z)
    double e1_x, e1_y, e1_z;
    if (std::fabs(dir_x) > std::fabs(dir_z)) {
        e1_x = -dir_y;
        e1_y =  dir_x;
        e1_z =  0.0;
    } else {
        e1_x =  0.0;
        e1_y = -dir_z;
        e1_z =  dir_y;
    }

    // Normalize e1
    double e1_norm = std::sqrt(e1_x*e1_x + e1_y*e1_y + e1_z*e1_z);
    e1_x /= e1_norm;
    e1_y /= e1_norm;
    e1_z /= e1_norm;

    // Compute e2 = dir × e1
    double e2_x = dir_y*e1_z - dir_z*e1_y;
    double e2_y = dir_z*e1_x - dir_x*e1_z;
    double e2_z = dir_x*e1_y - dir_y*e1_x;

    // --- 5. Apply scattering ---
    // u' = cosθ * dir + sinθ (cosφ * e1 + sinφ * e2)
    double new_x = dir_x * cos_theta + sin_theta * (cos_phi*e1_x + sin_phi*e2_x);
    double new_y = dir_y * cos_theta + sin_theta * (cos_phi*e1_y + sin_phi*e2_y);
    double new_z = dir_z * cos_theta + sin_theta * (cos_phi*e1_z + sin_phi*e2_z);

    return { new_x, new_y, new_z };
}



template <typename RandomGen>
double sampleEnergyLoss(double RSPlength, double E, RandomGen& gen, double E_real, float alpha){
    float gammaSquared = pc::gamma(E)*pc::gamma(E);
    float betaSquared = 1-1/gammaSquared;

    float stoppingPowerWater = pc::SP_CONSTANT_WATER/betaSquared * 
    (std::log(1e6*pc::ELECTRON_MASS*betaSquared*gammaSquared/pc::I_POT_WATER) - betaSquared); // MeV

    // Implement with ROOT
    double eta = pc::ETA_CONSTANT*RSPlength/betaSquared;
    double T_max = 2*pc::ELECTRON_MASS*betaSquared*gammaSquared/(1 + 2*std::sqrt(gammaSquared)*(pc::ELECTRON_MASS/pc::PROTON_MASS) + 
    (pc::ELECTRON_MASS/pc::PROTON_MASS)*(pc::ELECTRON_MASS/pc::PROTON_MASS));

    // float kappa = eta/T_max;

    // if (kappa > 10) {}... TODO WITH ROOT

    double dE_mean = RSPlength*stoppingPowerWater;
    if (E_real > 0) {                               // Sampling energy loss for shadow particle
        dE_mean *= (1 + alpha * (E - E_real)/E);    // Ensure correlation between particle and its shadow
    }
    double dE_std = std::sqrt(eta*T_max*(1-0.5*betaSquared));
    static thread_local std::normal_distribution<double> distEnergyLoss(dE_mean, dE_std);

    double dEsample = distEnergyLoss(gen);
    return dEsample;
}       