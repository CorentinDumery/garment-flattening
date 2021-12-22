

#include <Eigen/Core>
#include <iostream>

#include "procustes.h"

// One of 
//#include <Eigen/SparseCholesky>
//#include <Eigen/SparseQR>
#include <Eigen/IterativeLinearSolvers> // https://forum.kde.org/viewtopic.php?f=74&t=125165

//#define LOCALGLOBAL_DEBUG
#define LOCALGLOBAL_TIMING
//#define LOCALGLOBAL_DEBUG_SMALL // print full matrices, should be small

#ifdef LOCALGLOBAL_TIMING
#include <chrono>
using std::chrono::steady_clock;
using std::chrono::duration_cast;
using std::chrono::microseconds;
#endif

// Interesting discussion on Eigen performance for LSCM:
// https://forum.kde.org/viewtopic.php?f=74&t=125165

int next_equation_id = 0; // TODO put somewhere else

Eigen::VectorXd vertices2dToVector(const Eigen::MatrixXd& V){
    Eigen::VectorXd res(2 * V.rows());
    for (int i=0; i<V.rows(); i++){
        res(2 * i) = V(i,0);
        res(2 * i + 1) = V(i,1);
    }
    return res;
}

void equationsFromTriangle(const Eigen::MatrixXd& V_2d, const Eigen::MatrixXd& V_3d,
                           const Eigen::MatrixXi& F, int f_id,
                           std::vector<Eigen::Triplet<double>>& triplet_list,
                           std::vector<double>& target_vector){
    // each triangle gives us 6 equations:
    // each edge gives us one eq for U and one for V (TODO DONT DO EDGES TWICE?)

    Eigen::MatrixXd V_tri_2d = makeTriPoints(V_2d, F, f_id); // Could be removed for perf
    Eigen::MatrixXd V_tri_3d = makeTriPoints(V_3d, F, f_id);
    
    V_tri_3d = move3Dto2D(V_tri_3d);
    Eigen::MatrixXd R_est;
    Eigen::VectorXd T_est;
    procustes(V_tri_2d, V_tri_3d, R_est, T_est);

    // TODO check assumptions

    Eigen::MatrixXd p1 = V_tri_2d; // TODO get rid of extra notation
    Eigen::MatrixXd p2 = V_tri_3d;

    Eigen::MatrixXd p2_rt, p2_r;
    Eigen::MatrixXd p2t = p2.transpose(); // TODO transposeInPlace ?
    p2_rt = p2t.colwise() - T_est;
    p2_rt = (R_est.transpose() * p2_rt);
    p2_r = p2_rt.transpose();
    

    double ABu = (p1.row(1) - p1.row(0))(0);
    double ApBpu = (p2_r.row(1) - p2_r.row(0))(0);
    double ACu = (p1.row(2) - p1.row(0))(0);
    double ApCpu = (p2_r.row(2) - p2_r.row(0))(0);

    double ABv = (p1.row(1) - p1.row(0))(1);
    double ApBpv = (p2_r.row(1) - p2_r.row(0))(1);
    double ACv = (p1.row(2) - p1.row(0))(1);
    double ApCpv = (p2_r.row(2) - p2_r.row(0))(1);

    double Eu = std::pow(ABu - ApBpu, 2) 
                + std::pow(ACu - ApCpu, 2);
    double Ev = std::pow(ABv - ApBpv, 2) 
                + std::pow(ACv - ApCpv, 2);

    std::vector<std::pair<int, int>> edges = {std::make_pair(0,1), std::make_pair(0,2), std::make_pair(1,2)}; 

    for (std::pair<int, int> edge : edges){
        Eigen::RowVectorXd Ap = p2_r.row(edge.first);
        Eigen::RowVectorXd Bp = p2_r.row(edge.second);

        double target_u = (Bp - Ap)(0);
        triplet_list.push_back(Eigen::Triplet<double>(next_equation_id, 2 * F(f_id, edge.second), 1.0));
        triplet_list.push_back(Eigen::Triplet<double>(next_equation_id, 2 * F(f_id, edge.first), -1.0));
        target_vector.push_back(target_u);
        next_equation_id ++;

        double target_v = (Bp - Ap)(1);
        triplet_list.push_back(Eigen::Triplet<double>(next_equation_id, 2 * F(f_id, edge.second) + 1, 1.0));
        triplet_list.push_back(Eigen::Triplet<double>(next_equation_id, 2 * F(f_id, edge.first) + 1, -1.0));
        target_vector.push_back(target_v);
        next_equation_id ++;
    }
}

void makeSparseMatrix(const Eigen::MatrixXd& V_2d, const Eigen::MatrixXd& V_3d,
                      const Eigen::MatrixXi& F, const Eigen::MatrixXi& E,
                      Eigen::SparseMatrix<double>& A, Eigen::VectorXd& b){

    // Sparse conventions:
    // we have n target equations
    // and vector x of V 2D vertices
    // vertex i has its u coord in x(2*i) and its v coord in x(2*i+1)
    // M(equation, 2*v_id);

    next_equation_id = 0; // TODO REMOVE GLOBAL
    int n_equations = 6 * F.rows();// + 2;

    std::vector<Eigen::Triplet<double>> triplet_list; // Perf: get rid of std::vector
    std::vector<double> target_vector; // Perf: get rid of std::vector
    triplet_list.reserve(n_equations);
    for (int f_id=0; f_id<F.rows(); f_id++) {
        equationsFromTriangle(V_2d, V_3d, F, f_id, triplet_list, target_vector);
    }

    /* constrain in 0,0, only needed for non iterative methods?
    triplet_list.push_back(Eigen::Triplet<double>(next_equation_id, 0, 1.0));
    next_equation_id ++;
    triplet_list.push_back(Eigen::Triplet<double>(next_equation_id, 1, 1.0));
    next_equation_id ++;
    target_vector.push_back(0);
    target_vector.push_back(0);*/

    //b = Eigen::VectorXd(target_vector.size(), target_vector.data()); // Perf: get rid of std::vector
    b.resize(target_vector.size());
    for (int i=0; i<target_vector.size(); i++){ // TODO DO ANOTHER WAY
        b(i) = target_vector[i];
    }

    #ifdef LOCALGLOBAL_DEBUG
    std::cout << "n_equations: " << n_equations << std::endl;
    std::cout << "next_equation_id: " << next_equation_id << std::endl;
    std::cout << "triplet_list.size(): " << triplet_list.size() << std::endl; 
    std::cout << "target_vector.size(): " << target_vector.size() << std::endl; 
    #endif

    if (n_equations != target_vector.size()){
        std::cout << "ERROR: n_equations != b.rows(): " << n_equations << " vs " << target_vector.size() << std::endl;
    }

    if (n_equations != triplet_list.size()/2){// + 1){
        std::cout << "ERROR: n_equations != triplet_list.size()/2 + 1: " << n_equations << " vs " << triplet_list.size() << std::endl;
    }

    A.resize(n_equations, 2*V_2d.rows());
    A.setFromTriplets(triplet_list.begin(), triplet_list.end());

    #ifdef LOCALGLOBAL_DEBUG
    std::cout << "Sparse matrix computed." << std::endl;
    #endif
}

Eigen::MatrixXd localGlobal(const Eigen::MatrixXd& V_2d, const Eigen::MatrixXd& V_3d, 
                            const Eigen::MatrixXi& F, const Eigen::MatrixXi& E){

    #ifdef LOCALGLOBAL_TIMING
    steady_clock::time_point pre_localglobal = steady_clock::now();
    #endif

    Eigen::SparseMatrix<double> A;
    Eigen::VectorXd b, x;
    makeSparseMatrix(V_2d, V_3d, F, E, A, b);

    // solve Ax = b
    #ifdef LOCALGLOBAL_DEBUG
    std::cout << "Optimizing 2d vertices of size: " << V_2d.rows() << std::endl;
    std::cout << "Solver init..." << std::endl;
    #endif

    //Eigen::SimplicialLDLT <Eigen::SparseMatrix<double>> solver;
    //Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver(A);
    Eigen::LeastSquaresConjugateGradient<Eigen::SparseMatrix<double>> solver(A);
    //Eigen::SparseQR<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>> solver;
    //solver.setPivotThreshold(0.0f); //better performance if matrix has full rank
    //solver.compute(A);
    if(solver.info() != Eigen::Success) {
        std::cout << "ERROR: decomposition failed" << std::endl;
        //return;
    }
    //x.setZero();
    x = vertices2dToVector(V_2d);

    #ifdef LOCALGLOBAL_TIMING
    steady_clock::time_point pre_solve = steady_clock::now();
    #endif

    x = solver.solve(b);
    if(solver.info() != Eigen::Success) {
        std::cout << "ERROR: solving failed" << std::endl;
        //return;
    }

    #ifdef LOCALGLOBAL_TIMING
    steady_clock::time_point post_solve = steady_clock::now();
    #endif

    #ifdef LOCALGLOBAL_DEBUG
    std::cout << "Linear system success: " << std::endl;
    std::cout << "x.rows(): " << x.rows() << std::endl;
    #endif

    #ifdef LOCALGLOBAL_DEBUG_SMALL
    std::cout << "A" << std::endl;
    std::cout << A << std::endl;
    std::cout << "b" << std::endl;
    std::cout << b << std::endl;
    std::cout << x.topRows(10) << std::endl;
    std::cout << x << std::endl;
    #endif

    Eigen::MatrixXd res = Eigen::MatrixXd::Zero(V_2d.rows(), 3);
    for (int i=0; i<x.rows()/2; i++){
        res(i, 0) = x(2 * i); 
        res(i, 1) = x(2 * i + 1);
    }

    #ifdef LOCALGLOBAL_TIMING
    steady_clock::time_point post_localglobal = steady_clock::now();
    int pre_time = duration_cast<microseconds>(pre_solve - pre_localglobal).count();
    int solve_time = duration_cast<microseconds>(post_solve - pre_solve).count();
    int post_time = duration_cast<microseconds>(post_localglobal - post_solve).count();
    std::cout << "Precomp time : " << pre_time << " [µs]" << std::endl;
    std::cout << "Solving time : " << solve_time << " [µs]" << std::endl;
    std::cout << "Postcomp time: " << post_time << " [µs]" << std::endl << std::endl;
    #endif

    return res;
}