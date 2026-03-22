#include "pencilbeam.hpp"


namespace gc = grid_constants;
template struct PencilBeam<std::mt19937>;


template <typename RandomGen>
PencilBeam<RandomGen>::PencilBeam(float numPrimShare, double initialStep,
        RandomGen& gen, char entranceDir, float beamWidth, 
        double dir_x, double dir_y, double dir_z,
        double x_0, double y_0, double z_0,
        double E, double spread_E, float alpha)
    : numPrimShare(numPrimShare),
    initialStep(initialStep),
    gen(gen),
    entranceDir(entranceDir),
    beamWidth(beamWidth),
    dir_x(dir_x), dir_y(dir_y), dir_z(dir_z),
    x_0(x_0), y_0(y_0), z_0(z_0),
    E_0(E), spread_E(spread_E), alpha(alpha) {
    
    // Normalize direction vector
    double norm = std::sqrt(dir_x*dir_x + dir_y*dir_y + dir_z*dir_z);
    dir_x /= norm;
    dir_y /= norm;
    dir_z /= norm;

    assertWithMessage(validEntrance(x_0, y_0, z_0), "Invalid entrance location. Particles must enter at phantom boundary.");
    assertWithMessage(validEntranceDir(entranceDir), "Invalid entrance direction. Use 'x', 'y', or 'z");
    
    initBeamSpread();
}


template <typename RandomGen>
bool PencilBeam<RandomGen>::validEntrance(double x_0, double y_0, double z_0){
    return (x_0 == 0 || x_0 == gc::PHANTOM_LENGTH_X || y_0 == 0 || y_0 == gc::PHANTOM_LENGTH_Y || z_0 == 0 || z_0 == gc::PHANTOM_LENGTH_Z);
}


template <typename RandomGen>
bool PencilBeam<RandomGen>::validEntranceDir(char entranceDir){
    return (entranceDir == 'x' || entranceDir == 'y' || entranceDir == 'z');
}


template <typename RandomGen>
void PencilBeam<RandomGen>::initBeamSpread(){
    // 2D projection of gaussian beam spread on entrance inside phantom

    if (entranceDir == 'x'){
        spread_y = beamWidth*std::sqrt(dir_y*dir_y + dir_x*dir_x)/dir_x;
        spread_z = beamWidth*std::sqrt(dir_z*dir_z + dir_x*dir_x)/dir_x;
    }
    else if (entranceDir == 'y'){
        spread_x = beamWidth*std::sqrt(dir_x*dir_x + dir_y*dir_y)/dir_y;
        spread_z = beamWidth*std::sqrt(dir_z*dir_z + dir_y*dir_y)/dir_y;
    }
    else if (entranceDir == 'z'){
        spread_x = beamWidth*std::sqrt(dir_x*dir_x + dir_z*dir_z)/dir_z;
        spread_y = beamWidth*std::sqrt(dir_y*dir_y + dir_z*dir_z)/dir_z;
    }
    else{
        return;
    }
}



