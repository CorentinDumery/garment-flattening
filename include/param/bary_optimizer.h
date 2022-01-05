#pragma once

#include <Eigen/Core>
#include <Eigen/Sparse>
#include <vector>

#include "param/dart.h"

typedef Eigen::DiagonalMatrix<double, Eigen::Dynamic> DiagonalMatrixXd;

// Interesting discussion on Eigen performance for LSCM:
// https://forum.kde.org/viewtopic.php?f=74&t=125165

// benchmarking eigen solvers for |Ax - b|^2
// https://stackoverflow.com/questions/42116271/best-way-of-solving-sparse-linear-systems-in-c-gpu-possible

// WEIGHTED SOLVE: https://math.stackexchange.com/questions/709602/when-solving-an-overdetermined-linear-system-is-it-possible-to-weight-the-influ
// https://forum.kde.org/viewtopic.php?f=74&t=110784

class BaryOptimizer {
public:
    BaryOptimizer(int n_faces);

    Eigen::MatrixXd localGlobal(const Eigen::MatrixXd& V_2d, const Eigen::MatrixXd& V_3d, 
                                const Eigen::MatrixXi& F);

    void setSelectedVertices(std::vector<int> selected_vs) {selected_vs_ = selected_vs;};
    void setDarts(std::vector<SimpleDart> simple_darts) {
        for (SimpleDart d: simple_darts)
            simple_darts_.push_back(d);
    };
    void setDarts(std::vector<std::vector<int>> ordered_cuts) {
        std::vector<SimpleDart> simple_darts;
        for (int i=0; i<ordered_cuts.size(); i++){
            std::vector<int> cut = ordered_cuts[i]; 
            if (cut.size() % 2 == 0) continue;
            SimpleDart sd(cut);
            sd.print();
            simple_darts.push_back(sd);
        }

        setDarts(simple_darts);
    };

    // Config parameters
    bool enable_stretch_eqs_ = true;
    double stretch_coeff_ = 10.0;
    bool enable_angle_eqs_ = false;
    double angle_coeff_ = 0.7;
    bool enable_set_seed_eqs_ = true;
    bool enable_edges_eqs_ = true;
    double edges_coeff_ = 0.001;
    bool enable_selected_eqs_ = true;
    double selected_coeff_ = 1.0;
    bool enable_dart_sym_eqs_ = true;
    double dart_sym_coeff_ = 1.0;

private:

    int next_equation_id_ = 0;
    int n_equations_ = 0;
    int n_triplets_ = 0;
    std::vector<int> selected_vs_;
    std::vector<SimpleDart> simple_darts_;

    Eigen::VectorXd b;
    DiagonalMatrixXd W; 
    Eigen::MatrixXd V_tri_2d;
    Eigen::MatrixXd V_tri_3d;

    std::vector<Eigen::Triplet<double>> triplet_list; // Perf: get rid of std::vector
    std::vector<double> target_vector; // Perf: get rid of std::vector, replace with b
    std::vector<double> weight_vector; // Perf: get rid of std::vector, replace with W
    // TODO PERF: reserve for triplets?
    
    void equationsFromTriangle(const Eigen::MatrixXd& V_2d, const Eigen::MatrixXd& V_3d,
                               const Eigen::MatrixXi& F, int f_id,
                               std::vector<Eigen::Triplet<double>>& triplet_list,
                               std::vector<double>& target_vector, // TODO remove
                               std::vector<double>& weight_vector);

    void equationsFromDarts(const Eigen::MatrixXd& V_2d,
                            const Eigen::MatrixXi& F,
                            std::vector<Eigen::Triplet<double>>& triplet_list,
                            std::vector<double>& target_vector,
                            std::vector<double>& weight_vector);

    void makeSparseMatrix(const Eigen::MatrixXd& V_2d, const Eigen::MatrixXd& V_3d,
                      const Eigen::MatrixXi& F,
                      Eigen::SparseMatrix<double>& A, Eigen::VectorXd& b2,
                      DiagonalMatrixXd& W2);

    bool canUseSelectedEquation(){
        return enable_selected_eqs_ && selected_vs_.size() >= 2 && selected_vs_[0] >= 0 && selected_vs_[1] >= 0;
    }
    
    int edges_eq_time = 0;
    int stretch_shear_eq_time = 0;
    
};