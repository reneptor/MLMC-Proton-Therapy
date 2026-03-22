#include "voxeltraversal.hpp"


namespace gc = grid_constants;


Intersect::Intersect(double p_x, double p_y, double p_z, char c)
    : x(p_x), y(p_y), z(p_z), orientation(c) {}


LineSegment::LineSegment(double x1, double y1, double z1, double x2, double y2, double z2, 
    unsigned i, unsigned j, unsigned k)
    : x1(x1*gc::GRID_SPACING_X), y1(y1*gc::GRID_SPACING_Y), z1(z1*gc::GRID_SPACING_Z), 
    x2(x2*gc::GRID_SPACING_X), y2(y2*gc::GRID_SPACING_Y), z2(z2*gc::GRID_SPACING_Z), i(i), j(j), k(k) {}


float LineSegment::segmentLength(){
    return std::sqrt((x2-x1)*(x2-x1) + (y2-y1)*(y2-y1) + (z2-z1)*(z2-z1));
}


void GridPath::fillPlaneList(std::vector<unsigned int>& plane, double a1, double a2){
    unsigned int start = static_cast<unsigned int>(std::min(a1, a2)+1); // ceil
    unsigned int end = static_cast<unsigned int>(std::max(a1, a2)+1); // floor
    for (unsigned int i = start; i < end; ++i) {
        plane.push_back(i);
        n_intersects++;
    }
}


GridPath::GridPath(double x1_, double y1_, double z1_, double x2_, double y2_, double z2_) : 
    x1(x1_/gc::GRID_SPACING_X), y1(y1_/gc::GRID_SPACING_Y), z1(z1_/gc::GRID_SPACING_Z), 
    x2(x2_/gc::GRID_SPACING_X), y2(y2_/gc::GRID_SPACING_Y), z2(z2_/gc::GRID_SPACING_Z) {

    n_x = x2 - x1;
    n_y = y2 - y1;
    n_z = z2 - z1;

    planarConditionStart = n_x*x1 + n_y*y1 + n_z*z1;
    planarConditionEnd   = n_x*x2 + n_y*y2 + n_z*z2;

    fillPlaneList(x_planes, x1, x2);
    fillPlaneList(y_planes, y1, y2);
    fillPlaneList(z_planes, z1, z2);
}


std::array<unsigned int, 3> GridPath::findIntersectIndices(Intersect intersect){
    unsigned int i = static_cast<unsigned int>(intersect.x);
    unsigned int j = static_cast<unsigned int>(intersect.y);
    unsigned int k = static_cast<unsigned int>(intersect.z);

    if (intersect.orientation == 'x') {
        return (n_x > 0) ? std::array{i, j, k} : std::array{i - 1, j, k};
    } 
    else if (intersect.orientation == 'y') {
        return (n_y > 0) ? std::array{i, j, k} : std::array{i, j - 1, k};
    } 
    else if (intersect.orientation == 'z') {
        return (n_z > 0) ? std::array{i, j, k} : std::array{i, j, k - 1};
    }
    else {
        assertWithMessage(false, "Uninitialized orientation");
        return std::array<unsigned int, 3>{0, 0, 0}; // Should never get here based on prior logic
    }
}


std::array<double, 3> GridPath::pointAlongRayX(double x){
    double y_x = y1 + ((y2 - y1)/(x2 - x1)) * (x - x1);
    double z_x = z1 + ((z1 - z1)/(x2 - x1)) * (x - x1);
    return std::array{x, y_x, z_x};
}


std::array<double, 3> GridPath::pointAlongRayY(double y){
    double x_y = x1 + ((x2 - x1)/(y2 - y1)) * (y - y1);
    double z_y = z1 + ((z2 - z1)/(y2 - y1)) * (y - y1);
    return std::array{x_y, y, z_y};
}


std::array<double, 3> GridPath::pointAlongRayZ(double z){
    double x_z = x1 + ((x2 - x1)/(z2 - z1)) * (z - z1);
    double y_z = y1 + ((y2 - y1)/(z2 - z1)) * (z - z1);
    return std::array{x_z, y_z, z};
}


double GridPath::normalvecProduct(double p_x, double p_y, double p_z){
    return n_x*p_x + n_y*p_y + n_z*p_z;
}


bool GridPath::planarConditionSatisfied(double product){
    return ((planarConditionStart < product) && (product < planarConditionEnd));
}


void GridPath::decomposePath(){

    if (n_intersects == 0){ // Skipping expensive computaitons
        LineSegment segment(x1, y1, z1, x2, y2, z2,
        static_cast<unsigned int>(x1), static_cast<unsigned int>(y1), static_cast<unsigned int>(z1));
        traversalSegments.push_back(segment);
        return;
    }

    for (auto x : x_planes){
        auto [p_x, p_y, p_z] = pointAlongRayX(static_cast<double>(x));
        double product = normalvecProduct(p_x, p_y, p_z);
        if (planarConditionSatisfied(product)){
            Intersect intersect(p_x, p_y, p_z, 'x');
            planarIntersects.push_back(intersect);
        }
    }

    for (auto y : y_planes){
        auto [p_x, p_y, p_z] = pointAlongRayY(static_cast<double>(y));
        double product = normalvecProduct(p_x, p_y, p_z);
        if (planarConditionSatisfied(product)){
            Intersect intersect(p_x, p_y, p_z, 'y');
            planarIntersects.push_back(intersect);
        }
    }

    for (auto z : z_planes){
        auto [p_x, p_y, p_z] = pointAlongRayZ(static_cast<double>(z));
        double product = normalvecProduct(p_x, p_y, p_z);
        if (planarConditionSatisfied(product)){
            Intersect intersect(p_x, p_y, p_z, 'z');
            planarIntersects.push_back(intersect);
        }
    }

    if (planarIntersects.empty()){ // Band-aid solution. Sometimes unavoidable when points happen to be integers
        LineSegment segment(x1, y1, z1, x2, y2, z2,
        static_cast<unsigned int>((x2-x1)/2), static_cast<unsigned int>((y2-y1)/2), static_cast<unsigned int>((z2-z1)/2));
        traversalSegments.push_back(segment);
        return;
    }

    std::sort(planarIntersects.begin(), planarIntersects.end(), [this](Intersect int1, Intersect int2) {
        return (normalvecProduct(int1.x, int1.y, int1.z) < normalvecProduct(int2.x, int2.y, int2.z));  
    });
    
    Intersect firstIntersect = planarIntersects[0];
    LineSegment firstSegment(x1, y1, z1, firstIntersect.x, firstIntersect.y, firstIntersect.z,
    static_cast<unsigned int>(x1), static_cast<unsigned int>(y1), static_cast<unsigned int>(z1));
    traversalSegments.push_back(firstSegment);

    for (size_t t = 0; t < planarIntersects.size()-1; t++){
        Intersect int1 = planarIntersects[t];
        Intersect int2 = planarIntersects[t+1];
        auto [i, j, k] = findIntersectIndices(int1);
        LineSegment segment(int1.x, int1.y, int1.z, int2.x, int2.y, int2.z, i, j, k);
        traversalSegments.push_back(segment);
    }

    Intersect lastIntersect = planarIntersects[planarIntersects.size()-1];
    LineSegment lastSegment(lastIntersect.x, lastIntersect.y, lastIntersect.z, x2, y2, z2,
    static_cast<unsigned int>(x2), static_cast<unsigned int>(y2), static_cast<unsigned int>(z2));
    traversalSegments.push_back(lastSegment);
}


