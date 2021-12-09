#include <Eigen/Core>
#include <iostream>

Eigen::Matrix3d computeRotation(const Eigen::RowVector3d& from,
                                const Eigen::RowVector3d& to){
    // There might already be something like this in Eigen? Couldn't find it
    // https://math.stackexchange.com/questions/180418/calculate-rotation-matrix-to-align-vector-a-to-vector-b-in-3d
    Eigen::RowVector3d a = from.normalized();
    Eigen::RowVector3d b = to.normalized();
    if (a==b || a == -b) std::cout << "ERROR: case not handled in computeRotation" << std::endl;
    Eigen::RowVector3d v = a.cross(b);
    double s = v.norm();
    double c = a.dot(b);
    Eigen::Matrix3d vs;
    vs <<     0, -v[2],  v[1], 
           v[2],     0, -v[0], 
          -v[1],  v[0],     0;

    return Eigen::Matrix3d::Identity() + vs + vs * vs * 1.0 / (1.0+c);
}