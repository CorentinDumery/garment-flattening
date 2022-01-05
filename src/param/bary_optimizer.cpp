#include "param/bary_optimizer.h"

#include <Eigen/Core>
#include <iostream>

#include <Eigen/IterativeLinearSolvers> // https://forum.kde.org/viewtopic.php?f=74&t=125165

#include "param/dart.h"
#include "param/param_utils.h"

//#define LOCALGLOBAL_DEBUG
#define LOCALGLOBAL_TIMING
//#define LOCALGLOBAL_DEBUG_SMALL // print full matrices, should be small

#ifdef LOCALGLOBAL_TIMING
#include <chrono>
using std::chrono::steady_clock;
using std::chrono::duration_cast;
using std::chrono::microseconds;
#endif

#define USE_WEIGTHS_IN_LINEAR_SYSTEM true

BaryOptimizer::BaryOptimizer(int n_faces){
    // Sparse conventions:
    // we have n target equations
    // and vector x of V 2D vertices
    // vertex i has its u coord in x(2*i) and its v coord in x(2*i+1)
    // M(equation, 2*v_id);

    n_equations_ = 0; // Predict size for memory allocation
    n_triplets_ = 0;

    if (enable_stretch_eqs_) {
        n_equations_ += 2 * n_faces;
        n_triplets_ += 3 * 2 * n_faces;

        #ifdef LOCALGLOBAL_DEBUG
        std::cout << "stretch eqs: " << 2 * n_faces << std::endl;
        #endif
    }

    if (enable_angle_eqs_) {
        n_equations_ += 2 * n_faces;
        n_triplets_ += 3 * 2 * n_faces; 

        #ifdef LOCALGLOBAL_DEBUG
        std::cout << "angle  eqs: " << 2 * n_faces << std::endl;
        #endif
    }

    if (enable_set_seed_eqs_) {
        n_equations_ += 2;
        n_triplets_ += 2;

        #ifdef LOCALGLOBAL_DEBUG
        std::cout << "seed   eqs: " << 2 << std::endl;
        #endif
    }

    if (enable_edges_eqs_){
        n_equations_ += 6 * n_faces;
        n_triplets_ += 2 * 2 * n_faces;

        #ifdef LOCALGLOBAL_DEBUG
        std::cout << "edges  eqs: " << 6 * n_faces << std::endl;
        #endif
    }

    if (enable_selected_eqs_){ //canUseSelectedEquation()){
        n_equations_ += 1;
        n_triplets_ += 2;

        #ifdef LOCALGLOBAL_DEBUG
        std::cout << "selec eqs: " << 6 * n_faces << std::endl;
        #endif
    }

    if (enable_dart_sym_eqs_){
        int dart_points = 0; // points with dart target position, i.e. all on dart except tip
        for (int i=0; i<simple_darts_.size(); i++){
            dart_points += simple_darts_[i].size(); // -1;
        }
        n_equations_ += dart_points * 2;
        n_triplets_ += dart_points * 2;

        #ifdef LOCALGLOBAL_DEBUG
        std::cout << "dart sym eqs: " << dart_points * 2 << std::endl;
        #endif
    }

    triplet_list.resize(n_triplets_);
    target_vector.resize(n_equations_);
    weight_vector.resize(n_equations_);
    b.resize(n_equations_);
    W.resize(n_equations_);
    V_tri_2d.resize(3,3);
    V_tri_3d.resize(3,3);
}

void BaryOptimizer::equationsFromTriangle(const Eigen::MatrixXd& V_2d, const Eigen::MatrixXd& V_3d,
                           const Eigen::MatrixXi& F, int f_id,
                           std::vector<Eigen::Triplet<double>>& triplet_list,
                           std::vector<double>& target_vector,
                           std::vector<double>& weight_vector){
    
    // each triangle gives us 2 equations:
    // one target u and one target v

    #ifdef LOCALGLOBAL_TIMING
    steady_clock::time_point pre_stretch_shear_eqs = steady_clock::now();
    #endif

    //Eigen::MatrixXd V_tri_2d = makeTriPoints(V_2d, F, f_id); // Could be removed for perf
    //Eigen::MatrixXd V_3d.row(F(f_id, ke)TriPoints(V_3d, F, f_id);

    Eigen::RowVectorXd D = (V_2d.row(F(f_id, 0)) + V_2d.row(F(f_id, 1)) + V_2d.row(F(f_id, 2)))/3.0; // centroid
    Eigen::Vector3d D_bary(0.333333, 0.333333, 0.333333);
    Eigen::RowVectorXd DU = D;
    DU(0) += 1.0;
    Eigen::Vector3d DU_bary = barycentricCoords(DU, V_2d.row(F(f_id, 0)), V_2d.row(F(f_id, 1)), V_2d.row(F(f_id, 2)));
    Eigen::RowVectorXd DV = D;
    DV(1) += 1.0;
    Eigen::Vector3d DV_bary = barycentricCoords(DV, V_2d.row(F(f_id, 0)), V_2d.row(F(f_id, 1)), V_2d.row(F(f_id, 2)));
    Eigen::RowVectorXd DUV = D; // PERF: deduce DUV eq based on first two?
    DUV(0) += 1.0;
    DUV(1) += 1.0;
    Eigen::Vector3d DUV_bary = barycentricCoords(DUV, V_2d.row(F(f_id, 0)), V_2d.row(F(f_id, 1)), V_2d.row(F(f_id, 2)));

    Eigen::RowVectorXd Dp = D_bary(0) * V_3d.row(F(f_id, 0)) + D_bary(1) * V_3d.row(F(f_id, 1)) + D_bary(2) * V_3d.row(F(f_id, 2));
    Eigen::RowVectorXd DUp = DU_bary(0) * V_3d.row(F(f_id, 0)) + DU_bary(1) * V_3d.row(F(f_id, 1)) + DU_bary(2) * V_3d.row(F(f_id, 2));
    Eigen::RowVectorXd DVp = DV_bary(0) * V_3d.row(F(f_id, 0)) + DV_bary(1) * V_3d.row(F(f_id, 1)) + DV_bary(2) * V_3d.row(F(f_id, 2));
    Eigen::RowVectorXd DUVp = DUV_bary(0) * V_3d.row(F(f_id, 0)) + DUV_bary(1) * V_3d.row(F(f_id, 1)) + DUV_bary(2) * V_3d.row(F(f_id, 2));

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

    if (enable_stretch_eqs_){
        // Au
        triplet_list.push_back(Eigen::Triplet<double>(next_equation_id_, 2 * F(f_id, 0), DU_bary(0) - D_bary(0)));
        // Bu
        triplet_list.push_back(Eigen::Triplet<double>(next_equation_id_, 2 * F(f_id, 1), DU_bary(1) - D_bary(1)));
        // Cu
        triplet_list.push_back(Eigen::Triplet<double>(next_equation_id_, 2 * F(f_id, 2), DU_bary(2) - D_bary(2)));
        //target_vector.push_back(target_u);
        b(next_equation_id_) = target_u;
        //weight_vector.push_back(stretch_coeff_);
        W.diagonal()[next_equation_id_] = stretch_coeff_;
        next_equation_id_ ++;

        // Av
        triplet_list.push_back(Eigen::Triplet<double>(next_equation_id_, 2 * F(f_id, 0) + 1, DV_bary(0) - D_bary(0)));
        // Bv
        triplet_list.push_back(Eigen::Triplet<double>(next_equation_id_, 2 * F(f_id, 1) + 1, DV_bary(1) - D_bary(1)));
        // Cv
        triplet_list.push_back(Eigen::Triplet<double>(next_equation_id_, 2 * F(f_id, 2) + 1, DV_bary(2) - D_bary(2)));
        //target_vector.push_back(target_v);
        b(next_equation_id_) = target_v;
        //weight_vector.push_back(stretch_coeff_);
        W.diagonal()[next_equation_id_] = stretch_coeff_;
        next_equation_id_ ++;
    }

    if (enable_angle_eqs_){
        // Shear: for shear we actually need two: check how much diagonal changes on U, and on V
        // A_uv_u
        triplet_list.push_back(Eigen::Triplet<double>(next_equation_id_, 2 * F(f_id, 0), DUV_bary(0) - D_bary(0)));
        // B_uv_u
        triplet_list.push_back(Eigen::Triplet<double>(next_equation_id_, 2 * F(f_id, 1), DUV_bary(1) - D_bary(1)));
        // C_uv_u
        triplet_list.push_back(Eigen::Triplet<double>(next_equation_id_, 2 * F(f_id, 2), DUV_bary(2) - D_bary(2)));
        //target_vector.push_back(target_uv_u);
        b(next_equation_id_) = target_uv_u;
        //weight_vector.push_back(angle_coeff_);
        W.diagonal()[next_equation_id_] = angle_coeff_;
        next_equation_id_ ++;

        // A_uv_v
        triplet_list.push_back(Eigen::Triplet<double>(next_equation_id_, 2 * F(f_id, 0) + 1, DUV_bary(0) - D_bary(0)));
        // B_uv_v
        triplet_list.push_back(Eigen::Triplet<double>(next_equation_id_, 2 * F(f_id, 1) + 1, DUV_bary(1) - D_bary(1)));
        // C_uv_v
        triplet_list.push_back(Eigen::Triplet<double>(next_equation_id_, 2 * F(f_id, 2) + 1, DUV_bary(2) - D_bary(2)));
        //target_vector.push_back(target_uv_v);
        b(next_equation_id_) = target_uv_v;
        //weight_vector.push_back(angle_coeff_);
        W.diagonal()[next_equation_id_] = angle_coeff_;
        next_equation_id_ ++;
    }

    #ifdef LOCALGLOBAL_TIMING
    steady_clock::time_point pre_edges_eq = steady_clock::now();
    #endif

    if (enable_edges_eqs_){
        //Eigen::MatrixXd V_tri_2d = makeTriPoints(V_2d, F, f_id);
        //Eigen::MatrixXd V_tri_3d = makeTriPoints(V_3d, F, f_id);
        makeTriPoints(V_2d, F, f_id, V_tri_2d); // 5 ms
        makeTriPoints(V_3d, F, f_id, V_tri_3d); // 5 ms


        V_tri_3d = move3Dto2D(V_tri_3d); // 1 ms
        Eigen::MatrixXd R_est;
        Eigen::VectorXd T_est;
        procustes(V_tri_2d, V_tri_3d, R_est, T_est); // 4 ms 
        
    

        //Eigen::MatrixXd p1 = V_tri_2d; // TODO get rid of extra notation
        Eigen::MatrixXd p2 = V_tri_3d;

        Eigen::MatrixXd p2_rt, p2_r;
        Eigen::MatrixXd p2t = p2.transpose(); // TODO transposeInPlace ?
        p2_rt = p2t.colwise() - T_est;
        p2_rt = (R_est.transpose() * p2_rt);
        p2_r = p2_rt.transpose();

        std::vector<std::pair<int, int>> edges = {std::make_pair(0,1), std::make_pair(0,2), std::make_pair(1,2)}; 

        for (std::pair<int, int> edge : edges){
            Eigen::RowVectorXd Ap = p2_r.row(edge.first);
            Eigen::RowVectorXd Bp = p2_r.row(edge.second);

            double target_u = (Bp - Ap)(0);
            triplet_list.push_back(Eigen::Triplet<double>(next_equation_id_, 2 * F(f_id, edge.second), 1.0));
            triplet_list.push_back(Eigen::Triplet<double>(next_equation_id_, 2 * F(f_id, edge.first), -1.0));
            //target_vector.push_back(target_u);
            b(next_equation_id_) = target_u;
            //weight_vector.push_back(edges_coeff_);
            W.diagonal()[next_equation_id_] = edges_coeff_;
            next_equation_id_ ++;

            double target_v = (Bp - Ap)(1);
            triplet_list.push_back(Eigen::Triplet<double>(next_equation_id_, 2 * F(f_id, edge.second) + 1, 1.0));
            triplet_list.push_back(Eigen::Triplet<double>(next_equation_id_, 2 * F(f_id, edge.first) + 1, -1.0));
            //target_vector.push_back(target_v);
            b(next_equation_id_) = target_v;
            //    weight_vector.push_back(edges_coeff_);
            W.diagonal()[next_equation_id_] = edges_coeff_;
            next_equation_id_ ++;
        }
    }

    #ifdef LOCALGLOBAL_TIMING
    steady_clock::time_point post_edges_eq = steady_clock::now();
    stretch_shear_eq_time += duration_cast<microseconds>(pre_edges_eq - pre_stretch_shear_eqs).count();
    edges_eq_time += duration_cast<microseconds>(post_edges_eq - pre_edges_eq).count();
    #endif
}

void BaryOptimizer::equationsFromDarts(const Eigen::MatrixXd& V_2d,
                                       const Eigen::MatrixXi& F,
                                       std::vector<Eigen::Triplet<double>>& triplet_list,
                                       std::vector<double>& target_vector,
                                       std::vector<double>& weight_vector){
    for (auto dart: simple_darts_){
        //if (dart.size() < 3) continue;
        Eigen::RowVector3d sym_axis = dart.computeSymmetryAxis(V_2d);
        Eigen::MatrixXd targets = dart.computeSymmetryTargets(V_2d, sym_axis);

        if (dart.size() != targets.rows()){ 
            std::cout << "Error in equationsFromDarts" << std::endl;
        }

        for (int i=0; i<dart.size(); i++){
            int v_id = dart.v_id(i);

            double target_u = targets(i, 0);
            triplet_list.push_back(Eigen::Triplet<double>(next_equation_id_, 2 * v_id, 1.0));
            //target_vector.push_back(target_u);
            b(next_equation_id_) = target_u;
            weight_vector.push_back(dart_sym_coeff_);
            W.diagonal()[next_equation_id_] = dart_sym_coeff_;
            next_equation_id_ ++;

            double target_v = targets(i, 1);
            triplet_list.push_back(Eigen::Triplet<double>(next_equation_id_, 2 * v_id + 1, 1.0));
            //target_vector.push_back(target_v);
            b(next_equation_id_) = target_v;
            //weight_vector.push_back(dart_sym_coeff_);
            W.diagonal()[next_equation_id_] = dart_sym_coeff_;
            next_equation_id_ ++;
        }
        
    }
}

void BaryOptimizer::makeSparseMatrix(const Eigen::MatrixXd& V_2d, const Eigen::MatrixXd& V_3d,
                                     const Eigen::MatrixXi& F,
                                     Eigen::SparseMatrix<double>& A, Eigen::VectorXd& b2,
                                     DiagonalMatrixXd& W2){

    next_equation_id_ = 0;
    #ifdef LOCALGLOBAL_TIMING
    edges_eq_time = 0;
    stretch_shear_eq_time = 0;
    #endif

    triplet_list.clear(); // empty vector, but keeping memory size
    //target_vector.clear(); // Perf: get rid of std::vector
    //weight_vector.clear(); // Perf: get rid of std::vector

    for (int f_id=0; f_id<F.rows(); f_id++) {
        equationsFromTriangle(V_2d, V_3d, F, f_id, triplet_list, target_vector, weight_vector);
    }

    if (enable_set_seed_eqs_){
        triplet_list.push_back(Eigen::Triplet<double>(next_equation_id_, 0, 1.0));
        //target_vector.push_back(V_2d(0,0));
        b(next_equation_id_) = V_2d(0,0);
        //weight_vector.push_back(1.0);
        W.diagonal()[next_equation_id_] = 1.0;
        next_equation_id_ ++;


        triplet_list.push_back(Eigen::Triplet<double>(next_equation_id_, 1, 1.0));
        //target_vector.push_back(V_2d(0,1));
        b(next_equation_id_) = V_2d(0,1);
        //weight_vector.push_back(1.0);
        W.diagonal()[next_equation_id_] = 1.0;
        next_equation_id_ ++;
    }

    if (enable_selected_eqs_){
        if (selected_vs_.size() >= 2 && selected_vs_[0] >= 0 && selected_vs_[1] >= 0){
            // selected vertices should have = V
            int v0 = selected_vs_[0];
            int v1 = selected_vs_[1];
            triplet_list.push_back(Eigen::Triplet<double>(next_equation_id_, 2 * v0 + 1, 1.0));
            triplet_list.push_back(Eigen::Triplet<double>(next_equation_id_, 2 * v1 + 1, -1.0));
            //target_vector.push_back(0);
            b(next_equation_id_) = 0;
            //weight_vector.push_back(selected_coeff_);
            W.diagonal()[next_equation_id_] = selected_coeff_;
            next_equation_id_ ++;
        }
        else { // not enough vertices selected
            triplet_list.push_back(Eigen::Triplet<double>(next_equation_id_, 2 * 0 + 1, 1.0));
            triplet_list.push_back(Eigen::Triplet<double>(next_equation_id_, 2 * 0 + 1, -1.0));
            //target_vector.push_back(0);
            b(next_equation_id_) = 0;
            //weight_vector.push_back(0);
            W.diagonal()[next_equation_id_] = 0;
            next_equation_id_ ++;
        }
    }

    if (enable_dart_sym_eqs_){
        equationsFromDarts(V_2d, F, triplet_list, target_vector, weight_vector);
    }


    //b = Eigen::VectorXd(target_vector.size(), target_vector.data()); // Perf: get rid of std::vector
    /*b.resize(target_vector.size());
    for (int i=0; i<target_vector.size(); i++){ // TODO DO ANOTHER WAY
        b(i) = target_vector[i];
    }*/

    #ifdef LOCALGLOBAL_DEBUG
    std::cout << "n_equations_: " << n_equations_ << std::endl;
    std::cout << "next_equation_id_: " << next_equation_id_ << std::endl;
    std::cout << "triplet_list.size(): " << triplet_list.size() << std::endl; 
    std::cout << "target_vector.size(): " << target_vector.size() << std::endl; 
    #endif

    //if (n_equations_ != target_vector.size()){
    //    std::cout << "ERROR: n_equations_ != b.rows(): " << n_equations_ << " vs " << target_vector.size() << std::endl;
    //}

    if (n_equations_ != next_equation_id_){
        std::cout << "ERROR: n_equations_ != next_equation_id_: " << n_equations_ << " vs " << next_equation_id_ << std::endl;
    }

    /*if (n_equations_ != triplet_list.size()/3){// + 1){
        std::cout << "ERROR: n_equations_ != triplet_list.size()/3: " << n_equations_ << " vs " << triplet_list.size() << std::endl;
    }*/

    A.resize(n_equations_, 2*V_2d.rows());
    A.setFromTriplets(triplet_list.begin(), triplet_list.end());

    /*if (USE_WEIGTHS_IN_LINEAR_SYSTEM){
        //W.resize(n_equations_, n_equations_);
        //W.setFromTriplets(weight_triplets.begin(), weight_triplets.end());
        //W.resize(n_equations_);
        Eigen::VectorXd temp(n_equations_);
        for (int i=0; i<n_equations_; i++){
            temp(i) = weight_vector[i]; // TODO perf
        }
        W = temp.asDiagonal();
    }*/

    #ifdef LOCALGLOBAL_DEBUG
    std::cout << "Sparse matrix computed." << std::endl;
    #endif
}

Eigen::MatrixXd BaryOptimizer::localGlobal(const Eigen::MatrixXd& V_2d, const Eigen::MatrixXd& V_3d, 
                                           const Eigen::MatrixXi& F){

    #ifdef LOCALGLOBAL_TIMING
    steady_clock::time_point pre_localglobal = steady_clock::now();
    #endif

    Eigen::SparseMatrix<double> A;
    DiagonalMatrixXd W2; // TODO remove
    Eigen::VectorXd b2, x;
    makeSparseMatrix(V_2d, V_3d, F, A, b2, W2);
    DiagonalMatrixXd Wt = W;
    if (USE_WEIGTHS_IN_LINEAR_SYSTEM){
        //Wt = Wt.transpose();
    }

    #ifdef LOCALGLOBAL_DEBUG_SMALL
    if (USE_WEIGTHS_IN_LINEAR_SYSTEM){
        std::cout << "W's diagonal coefficients: " << std::endl;
        std::cout << W.diagonal() << std::endl;
    }
    #endif

    // solve Ax = b
    #ifdef LOCALGLOBAL_DEBUG
    std::cout << "Optimizing 2d vertices of size: " << V_2d.rows() << std::endl;
    std::cout << "Solver init..." << std::endl;
    #endif

    Eigen::SparseMatrix<double> Ap = A;
    Eigen::SparseMatrix<double> At = A;
    At = At.transpose();
    if (USE_WEIGTHS_IN_LINEAR_SYSTEM){
        Ap = At * Wt * W * A;
        //Ap = W * A;
    }
    
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
    solver.compute(Ap);

    //Eigen::LeastSquaresConjugateGradient<Eigen::SparseMatrix<double>> solver(Ap);
    
    //Eigen::SparseQR<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>> solver;
    //solver.setPivotThreshold(0.0f); //better performance if matrix has full rank
    //solver.compute(A); //for sparseQR

    if(solver.info() != Eigen::Success) {
        std::cout << "ERROR: decomposition failed" << std::endl;
        //return;
    }

    x = vertices2dToVector(V_2d); // Initial solution
    x = Eigen::VectorXd::Zero(V_2d.rows() * 2);

    Eigen::VectorXd bp = b;
    if (USE_WEIGTHS_IN_LINEAR_SYSTEM){
        bp = At * Wt * W * b;
        //bp = W * b;
    }

    #ifdef LOCALGLOBAL_TIMING
    steady_clock::time_point pre_solve = steady_clock::now();
    #endif

    x = solver.solve(bp);

    if(solver.info() != Eigen::Success) {
        std::cout << "ERROR, solving failed: ";
        if(solver.info() == Eigen::NumericalIssue) 
            std::cout << "NumericalIssue" << std::endl;
        if(solver.info() == Eigen::NoConvergence) 
            std::cout << "NoConvergence" << std::endl;
        if(solver.info() == Eigen::InvalidInput) 
            std::cout << "InvalidInput" << std::endl;
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
    std::cout << "# of non-zero elements in linear system: " << Ap.nonZeros() << std::endl;
    steady_clock::time_point post_localglobal = steady_clock::now();
    int pre_time = duration_cast<microseconds>(pre_solve - pre_localglobal).count();
    int solve_time = duration_cast<microseconds>(post_solve - pre_solve).count();
    int post_time = duration_cast<microseconds>(post_localglobal - post_solve).count();
    std::cout << "Precomp time : " << pre_time << " [µs]" << std::endl;
    std::cout << "\t incl. stretch_shear_eq_time: " << stretch_shear_eq_time << " [µs]" << std::endl;
    std::cout << "\t incl. edges_eq_time: " << edges_eq_time << " [µs]" << std::endl;
    std::cout << "Solving time : " << solve_time << " [µs]" << std::endl;
    std::cout << "Postcomp time: " << post_time << " [µs]" << std::endl << std::endl;
    #endif

    return res;
}