#pragma once

#include "simulation_constants.hpp"
#include "error_logging.hpp"
#include <vector>
#include <cmath>
#include <algorithm>


struct Intersect {
    double x {}, y {}, z {};
    char orientation;

    Intersect(double p_x, double p_y, double p_z, char c);
};


struct LineSegment {
    double x1 {}, y1 {}, z1 {};
    double x2 {}, y2 {}, z2 {};
    unsigned int i {}, j {}, k {};

    LineSegment(double x1, double y1, double z1, double x2, double y2, double z2,
                unsigned i, unsigned j, unsigned k);

    float segmentLength();
};


struct GridPath {
    double x1 {}, y1 {}, z1 {};
    double x2 {}, y2 {}, z2 {};
    double n_x {}, n_y {}, n_z {};
    double planarConditionStart {}, planarConditionEnd {};
    unsigned int n_intersects = 0;

    std::vector<Intersect> planarIntersects;
    std::vector<LineSegment> traversalSegments;

    std::vector<unsigned int> x_planes;
    std::vector<unsigned int> y_planes;
    std::vector<unsigned int> z_planes;

    GridPath(double x1_, double y1_, double z1_, double x2_, double y2_, double z2_);

    void fillPlaneList(std::vector<unsigned int>& plane, double a1, double a2);
    std::array<unsigned int, 3> findIntersectIndices(Intersect intersect);
    std::array<double, 3> pointAlongRayX(double x);
    std::array<double, 3> pointAlongRayY(double y);
    std::array<double, 3> pointAlongRayZ(double z);
    double normalvecProduct(double p_x, double p_y, double p_z);
    bool planarConditionSatisfied(double product);
    void decomposePath();

};





