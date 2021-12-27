

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

#define USE_WEIGTHS_IN_LINEAR_SYSTEM true

Eigen::VectorXd vertices2dToVector(const Eigen::MatrixXd& V){ // TODO put somewhere else
    Eigen::VectorXd res(2 * V.rows()); // TODO faster?
    for (int i=0; i<V.rows(); i++){
        res(2 * i) = V(i,0);
        res(2 * i + 1) = V(i,1);
    }
    return res;
}

// Compute barycentric coordinates (u, v, w) for
// point p with respect to triangle (a, b, c)
// credits https://gamedev.stackexchange.com/questions/23743/whats-the-most-efficient-way-to-find-barycentric-coordinates
Eigen::Vector3d barycentricCoords(const Eigen::RowVector3d& p, const Eigen::RowVector3d& a, 
                                         const Eigen::RowVector3d& b, const Eigen::RowVector3d& c){ // TODO put somewhere else
    Eigen::RowVector3d v0 = b - a;
    Eigen::RowVector3d v1 = c - a;
    Eigen::RowVector3d v2 = p - a;
    float d00 = v0.dot(v0);
    float d01 = v0.dot(v1);
    float d11 = v1.dot(v1);
    float d20 = v2.dot(v0);
    float d21 = v2.dot(v1);
    float denom = d00 * d11 - d01 * d01;
    double v = (d11 * d20 - d01 * d21) / denom;
    double w = (d00 * d21 - d01 * d20) / denom;
    double u = 1.0f - v - w;
    return Eigen::Vector3d(u, v, w);
}


void equationsFromTriangle(const Eigen::MatrixXd& V_2d, const Eigen::MatrixXd& V_3d,
                           const Eigen::MatrixXi& F, int f_id,
                           std::vector<Eigen::Triplet<double>>& triplet_list,
                           std::vector<double>& target_vector,
                           std::vector<Eigen::Triplet<double>>& weight_triplets){
    
    // each triangle gives us 2 equations:
    // one target u and one target v

    Eigen::MatrixXd V_tri_2d = makeTriPoints(V_2d, F, f_id); // Could be removed for perf
    Eigen::MatrixXd V_tri_3d = makeTriPoints(V_3d, F, f_id);

    Eigen::RowVectorXd D = (V_tri_2d.row(0) + V_tri_2d.row(1) + V_tri_2d.row(2))/3.0; // centroid
    Eigen::Vector3d D_bary(0.333333, 0.333333, 0.333333);
    Eigen::RowVectorXd DU = D;
    DU(0) += 1.0;
    Eigen::Vector3d DU_bary = barycentricCoords(DU, V_tri_2d.row(0), V_tri_2d.row(1), V_tri_2d.row(2));
    Eigen::RowVectorXd DV = D;
    DV(1) += 1.0;
    Eigen::Vector3d DV_bary = barycentricCoords(DV, V_tri_2d.row(0), V_tri_2d.row(1), V_tri_2d.row(2));
    Eigen::RowVectorXd DUV = D; // PERF: deduce DUV eq based on first two?
    DUV(0) += 1.0;
    DUV(1) += 1.0;
    Eigen::Vector3d DUV_bary = barycentricCoords(DUV, V_tri_2d.row(0), V_tri_2d.row(1), V_tri_2d.row(2));

    Eigen::RowVectorXd Dp = D_bary(0) * V_tri_3d.row(0) + D_bary(1) * V_tri_3d.row(1) + D_bary(2) * V_tri_3d.row(2);
    Eigen::RowVectorXd DUp = DU_bary(0) * V_tri_3d.row(0) + DU_bary(1) * V_tri_3d.row(1) + DU_bary(2) * V_tri_3d.row(2);
    Eigen::RowVectorXd DVp = DV_bary(0) * V_tri_3d.row(0) + DV_bary(1) * V_tri_3d.row(1) + DV_bary(2) * V_tri_3d.row(2);
    Eigen::RowVectorXd DUVp = DUV_bary(0) * V_tri_3d.row(0) + DUV_bary(1) * V_tri_3d.row(1) + DUV_bary(2) * V_tri_3d.row(2);

    double target_u = (DUp - Dp).norm();
    double target_v = (DVp - Dp).norm();
    //double target_uv_u = (DUVp - Dp)(0); No this doesnt work: you need to project onto 3D UV vectors
    //double target_uv_v = (DUVp - Dp)(1);
    double target_uv_u = (DUVp - Dp).dot(DUp - Dp)/(DUp - Dp).norm();
    double target_uv_v = (DUVp - Dp).dot(DVp - Dp)/(DVp - Dp).norm();;

    #ifdef LOCALGLOBAL_DEBUG
    std::cout << "bary centroid: " << barycentricCoords(D, V_tri_2d.row(0), V_tri_2d.row(1), V_tri_2d.row(2)) << std::endl;
    std::cout << "bary DU: " << DU_bary << std::endl;
    std::cout << "bary DV: " << DV_bary << std::endl;
    std::cout << "bary DUV: " << DUV_bary << std::endl;
    std::cout << "target_u: " << target_u << std::endl;
    std::cout << "target_v: " << target_v << std::endl;
    std::cout << "target_uv_u: " << target_uv_u << std::endl;
    std::cout << "target_uv_v: " << target_uv_v << std::endl;
    #endif

    // Au
    triplet_list.push_back(Eigen::Triplet<double>(next_equation_id, 2 * F(f_id, 0), DU_bary(0) - D_bary(0)));
    // Bu
    triplet_list.push_back(Eigen::Triplet<double>(next_equation_id, 2 * F(f_id, 1), DU_bary(1) - D_bary(1)));
    // Cu
    triplet_list.push_back(Eigen::Triplet<double>(next_equation_id, 2 * F(f_id, 2), DU_bary(2) - D_bary(2)));
    target_vector.push_back(target_u);
    if (USE_WEIGTHS_IN_LINEAR_SYSTEM)
        weight_triplets.push_back(Eigen::Triplet<double>(next_equation_id, next_equation_id, 1.0));
    next_equation_id ++;

    // Av
    triplet_list.push_back(Eigen::Triplet<double>(next_equation_id, 2 * F(f_id, 0) + 1, DV_bary(0) - D_bary(0)));
    // Bv
    triplet_list.push_back(Eigen::Triplet<double>(next_equation_id, 2 * F(f_id, 1) + 1, DV_bary(1) - D_bary(1)));
    // Cv
    triplet_list.push_back(Eigen::Triplet<double>(next_equation_id, 2 * F(f_id, 2) + 1, DV_bary(2) - D_bary(2)));
    target_vector.push_back(target_v);
    if (USE_WEIGTHS_IN_LINEAR_SYSTEM)
        weight_triplets.push_back(Eigen::Triplet<double>(next_equation_id, next_equation_id, 1.0));
    next_equation_id ++;

    // Shear: for shear we actually need two: check how much diagonal changes on U, and on V
    //*
    // A_uv_u
    triplet_list.push_back(Eigen::Triplet<double>(next_equation_id, 2 * F(f_id, 0), DUV_bary(0) - D_bary(0)));
    // B_uv_u
    triplet_list.push_back(Eigen::Triplet<double>(next_equation_id, 2 * F(f_id, 1), DUV_bary(1) - D_bary(1)));
    // C_uv_u
    triplet_list.push_back(Eigen::Triplet<double>(next_equation_id, 2 * F(f_id, 2), DUV_bary(2) - D_bary(2)));
    target_vector.push_back(target_uv_u);
    if (USE_WEIGTHS_IN_LINEAR_SYSTEM)
        weight_triplets.push_back(Eigen::Triplet<double>(next_equation_id, next_equation_id, 1.3));
    next_equation_id ++;

    // A_uv_v
    triplet_list.push_back(Eigen::Triplet<double>(next_equation_id, 2 * F(f_id, 0) + 1, DUV_bary(0) - D_bary(0)));
    // B_uv_v
    triplet_list.push_back(Eigen::Triplet<double>(next_equation_id, 2 * F(f_id, 1) + 1, DUV_bary(1) - D_bary(1)));
    // C_uv_v
    triplet_list.push_back(Eigen::Triplet<double>(next_equation_id, 2 * F(f_id, 2) + 1, DUV_bary(2) - D_bary(2)));
    target_vector.push_back(target_uv_v);
    if (USE_WEIGTHS_IN_LINEAR_SYSTEM)
        weight_triplets.push_back(Eigen::Triplet<double>(next_equation_id, next_equation_id, 1.3));
    next_equation_id ++;

    //*/
}

void makeSparseMatrix(const Eigen::MatrixXd& V_2d, const Eigen::MatrixXd& V_3d,
                      const Eigen::MatrixXi& F, const Eigen::MatrixXi& E,
                      Eigen::SparseMatrix<double>& A, Eigen::VectorXd& b,
                      Eigen::SparseMatrix<double>& W){

    // Sparse conventions:
    // we have n target equations
    // and vector x of V 2D vertices
    // vertex i has its u coord in x(2*i) and its v coord in x(2*i+1)
    // M(equation, 2*v_id);

    next_equation_id = 0; // TODO REMOVE GLOBAL
    int n_equations = 4 * F.rows();

    std::vector<Eigen::Triplet<double>> triplet_list; // Perf: get rid of std::vector
    std::vector<double> target_vector; // Perf: get rid of std::vector
    triplet_list.reserve(3*n_equations);
    std::vector<Eigen::Triplet<double>> weight_triplets; // Perf: get rid of std::vector
    // TODO PERF: reserve for triplets?
    for (int f_id=0; f_id<F.rows(); f_id++) {
        equationsFromTriangle(V_2d, V_3d, F, f_id, triplet_list, target_vector, weight_triplets);
    }

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

    if (n_equations != triplet_list.size()/3){// + 1){
        std::cout << "ERROR: n_equations != triplet_list.size()/3: " << n_equations << " vs " << triplet_list.size() << std::endl;
    }

    A.resize(n_equations, 2*V_2d.rows());
    A.setFromTriplets(triplet_list.begin(), triplet_list.end());

    if (USE_WEIGTHS_IN_LINEAR_SYSTEM){
        W.resize(n_equations, n_equations);
        W.setFromTriplets(weight_triplets.begin(), weight_triplets.end());
    }

    #ifdef LOCALGLOBAL_DEBUG
    std::cout << "Sparse matrix computed." << std::endl;
    #endif
}

Eigen::MatrixXd localGlobal(const Eigen::MatrixXd& V_2d, const Eigen::MatrixXd& V_3d, 
                            const Eigen::MatrixXi& F, const Eigen::MatrixXi& E){

    #ifdef LOCALGLOBAL_TIMING
    steady_clock::time_point pre_localglobal = steady_clock::now();
    #endif

    Eigen::SparseMatrix<double> A, W;
    Eigen::VectorXd b, x;
    makeSparseMatrix(V_2d, V_3d, F, E, A, b, W);
    Eigen::SparseMatrix<double> Wt = W;
    if (USE_WEIGTHS_IN_LINEAR_SYSTEM){
        Wt = Wt.transpose();
    }

    // WEIGHTED SOLVE: https://math.stackexchange.com/questions/709602/when-solving-an-overdetermined-linear-system-is-it-possible-to-weight-the-influ
    // https://forum.kde.org/viewtopic.php?f=74&t=110784

    #ifdef LOCALGLOBAL_DEBUG_SMALL
    if (USE_WEIGTHS_IN_LINEAR_SYSTEM){
        std::cout << "W" << std::endl;
        std::cout << W << std::endl;
    }
    #endif

    // solve Ax = b
    #ifdef LOCALGLOBAL_DEBUG
    std::cout << "Optimizing 2d vertices of size: " << V_2d.rows() << std::endl;
    std::cout << "Solver init..." << std::endl;
    #endif

    Eigen::SparseMatrix<double> Ap = A;
    if (USE_WEIGTHS_IN_LINEAR_SYSTEM){
        Eigen::SparseMatrix<double> At = A;
        At = At.transpose();
        Ap = At * W.transpose() * W * A;
        //Ap = W * A;
    }
    //Eigen::SimplicialLDLT <Eigen::SparseMatrix<double>> solver;
    //Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver(A);
    Eigen::LeastSquaresConjugateGradient<Eigen::SparseMatrix<double>> solver(Ap);
    //Eigen::SparseQR<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>> solver;
    //solver.setPivotThreshold(0.0f); //better performance if matrix has full rank
    //solver.compute(A); //for sparseQR
    if(solver.info() != Eigen::Success) {
        std::cout << "ERROR: decomposition failed" << std::endl;
        //return;
    }

    x = vertices2dToVector(V_2d); // Initial solution

    #ifdef LOCALGLOBAL_TIMING
    steady_clock::time_point pre_solve = steady_clock::now();
    #endif

    Eigen::VectorXd bp = b;
    if (USE_WEIGTHS_IN_LINEAR_SYSTEM){
        bp = A.transpose() * W.transpose() * W * b;
        //bp = W * b;
    }
    x = solver.solve(bp);
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