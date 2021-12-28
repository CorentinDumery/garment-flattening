#pragma once

#include <Eigen/Core>
#include <Eigen/Sparse>
#include <vector>

typedef Eigen::DiagonalMatrix<double, Eigen::Dynamic> DiagonalMatrixXd;

// Interesting discussion on Eigen performance for LSCM:
// https://forum.kde.org/viewtopic.php?f=74&t=125165

// benchmarking eigen solvers for |Ax - b|^2
// https://stackoverflow.com/questions/42116271/best-way-of-solving-sparse-linear-systems-in-c-gpu-possible

// WEIGHTED SOLVE: https://math.stackexchange.com/questions/709602/when-solving-an-overdetermined-linear-system-is-it-possible-to-weight-the-influ
// https://forum.kde.org/viewtopic.php?f=74&t=110784

class BaryOptimizer {
public:
    Eigen::MatrixXd localGlobal(const Eigen::MatrixXd& V_2d, const Eigen::MatrixXd& V_3d, 
                                const Eigen::MatrixXi& F);

    // Config parameters
    bool enable_stretch_eqs_ = true;
    double stretch_coeff_ = 10.0;
    bool enable_angle_eqs_ = false;
    double angle_coeff_ = 0.7;
    bool enable_set_seed_eqs_ = true;
    bool enable_edges_eqs_ = true;
    double edges_coeff_ = 0.001;

private:

    int next_equation_id_ = 0;
    
    void equationsFromTriangle(const Eigen::MatrixXd& V_2d, const Eigen::MatrixXd& V_3d,
                               const Eigen::MatrixXi& F, int f_id,
                               std::vector<Eigen::Triplet<double>>& triplet_list,
                               std::vector<double>& target_vector,
                               std::vector<double>& weight_vector);

    void makeSparseMatrix(const Eigen::MatrixXd& V_2d, const Eigen::MatrixXd& V_3d,
                      const Eigen::MatrixXi& F,
                      Eigen::SparseMatrix<double>& A, Eigen::VectorXd& b,
                      DiagonalMatrixXd& W);

    
};