#pragma once
#include <Eigen/Core>
#include <vector>

class SurfaceNet {

public:
    SurfaceNet(const Eigen::MatrixXi& F,
               const Eigen::MatrixXf& V_3d,
               const Eigen::MatrixXf& V_2d);

    void computeNet();
    void vizBoundaryEdges(Eigen::MatrixXd& edge_begs, Eigen::MatrixXd& edge_ends) const;
    std::vector<std::vector<Eigen::MatrixXd>> vizNet() const {return {fiber_begs_list, fiber_ends_list};}

    std::vector<std::vector<int>> nearestFibers(const Eigen::RowVector2f& vertex) const; 

private:
    
    // Set by constructor
    const Eigen::MatrixXi F_;
    Eigen::MatrixXf V_3d_;
    Eigen::MatrixXf V_2d_; // Note: ideally we'd get rid of that and directly work on something closer to the rendering buffers
    Eigen::MatrixXi Eb_; // boundary edges
    Eigen::MatrixXi TT_; // triangle-triangle adjacency
    int n_tris_;
    
    // Net fiber variables
    std::vector<Eigen::MatrixXd> fiber_begs_list; // Fiber viz
    std::vector<Eigen::MatrixXd> fiber_ends_list; 
};