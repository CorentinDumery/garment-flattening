#include <Eigen/Core>
#include <Eigen/Sparse>
#include <Eigen/IterativeLinearSolvers> // https://forum.kde.org/viewtopic.php?f=74&t=125165
#include <vector>

#include <igl/readOBJ.h>
#include <igl/writeOBJ.h>
#include <igl/boundary_loop.h>

#include "param/param_utils.h"

typedef Eigen::DiagonalMatrix<double, Eigen::Dynamic> DiagonalMatrixXd;

void makeSparseMatrix(const Eigen::MatrixXd& V_2d, const Eigen::MatrixXd& V_3d,
                      const Eigen::MatrixXi& F,
                      Eigen::SparseMatrix<double>& A,
                      Eigen::VectorXd& b, 
                      DiagonalMatrixXd& W,
                      Eigen::VectorXd& x,
                      std::vector<int> middle_ids,
                      std::vector<int> sleeve_ids,
                      const Eigen::MatrixXd& sleeve_pos){

    int next_equation_id = 0;
    
    std::vector<Eigen::Triplet<double>> triplet_list; 
    
    std::vector<double> b_vals;
    std::vector<double> w_vals;

    Eigen::MatrixXd V_tri_2d(3,3);
    Eigen::MatrixXd V_tri_3d(3,3);

    for (int f_id=0; f_id<F.rows(); f_id++) {
        bool enable_edges_eqs = true;
        double edges_coeff = 1.0;
        if (enable_edges_eqs){
            makeTriPoints(V_2d, F, f_id, V_tri_2d);
            makeTriPoints(V_3d, F, f_id, V_tri_3d);

            V_tri_3d = move3Dto2D(V_tri_3d); 
            Eigen::MatrixXd R_est;
            Eigen::VectorXd T_est;
            procrustes(V_tri_2d, V_tri_3d, R_est, T_est); 
            
            Eigen::MatrixXd p2 = V_tri_3d;
            Eigen::MatrixXd p2_rt, p2_r;
            Eigen::MatrixXd p2t = p2.transpose();
            p2_rt = p2t.colwise() - T_est;
            p2_rt = (R_est.transpose() * p2_rt);
            p2_r = p2_rt.transpose();

            std::vector<std::pair<int, int>> edges = {std::make_pair(0,1), std::make_pair(0,2), std::make_pair(1,2)}; 

            for (std::pair<int, int> edge : edges){
                Eigen::RowVectorXd Ap = p2_r.row(edge.first);
                Eigen::RowVectorXd Bp = p2_r.row(edge.second);

                double target_u = (Bp - Ap)(0);
                triplet_list.push_back(Eigen::Triplet<double>(next_equation_id, 2 * F(f_id, edge.second), 1.0));
                triplet_list.push_back(Eigen::Triplet<double>(next_equation_id, 2 * F(f_id, edge.first), -1.0));
                b_vals.push_back(target_u);
                w_vals.push_back(edges_coeff);
                next_equation_id ++;

                double target_v = (Bp - Ap)(1);
                triplet_list.push_back(Eigen::Triplet<double>(next_equation_id, 2 * F(f_id, edge.second) + 1, 1.0));
                triplet_list.push_back(Eigen::Triplet<double>(next_equation_id, 2 * F(f_id, edge.first) + 1, -1.0));
                b_vals.push_back(target_v);
                w_vals.push_back(edges_coeff);
                next_equation_id ++;
            }
        }
    }

    bool enable_set_seed_eqs = false;
    if (enable_set_seed_eqs){
        double seed_coeff = 1.0;
        triplet_list.push_back(Eigen::Triplet<double>(next_equation_id, 0, 1.0));
        b_vals.push_back(V_2d(0,0));
        w_vals.push_back(seed_coeff);
        next_equation_id ++;

        triplet_list.push_back(Eigen::Triplet<double>(next_equation_id, 1, 1.0));
        b_vals.push_back(V_2d(0,1));
        w_vals.push_back(seed_coeff);
        next_equation_id ++;
    }

    bool enable_middle_constraint = true;
    double middle_coeff = 1000.0;
    if (enable_middle_constraint){
        for (int i: middle_ids){
            triplet_list.push_back(Eigen::Triplet<double>(next_equation_id, 2 * i , 1.0));
            b_vals.push_back(0);
            w_vals.push_back(middle_coeff);
            next_equation_id ++;
        }
    }

    bool enable_sleeve_constraint = true;
    double sleeve_coeff = 1000.0;
    if (enable_sleeve_constraint){
        for (int i=0; i<sleeve_ids.size(); i++){
            triplet_list.push_back(Eigen::Triplet<double>(next_equation_id, 2 * sleeve_ids[i], 1.0));
            b_vals.push_back(sleeve_pos(i, 0));
            w_vals.push_back(sleeve_coeff);
            next_equation_id ++;

            triplet_list.push_back(Eigen::Triplet<double>(next_equation_id, 2 * sleeve_ids[i] + 1, 1.0));
            b_vals.push_back(sleeve_pos(i, 1));
            w_vals.push_back(sleeve_coeff);
            next_equation_id ++;
        }
    }

    A.resize(next_equation_id, 2 * V_2d.rows());
    A.setFromTriplets(triplet_list.begin(), triplet_list.end());

    if (b_vals.size() != next_equation_id) std::cout << "ERROR: wrong b size" << std::endl; 

    b.resize(next_equation_id);
    for (int i=0; i<next_equation_id; i++){ // TODO faster one liner
        b(i) = b_vals[i];
    }

    if (w_vals.size() != next_equation_id) std::cout << "ERROR: wrong W size" << std::endl; 
    W.resize(next_equation_id); 
    for (int i=0; i<next_equation_id; i++){ // TODO faster?
        W.diagonal()[i] = w_vals[i];
    }

    x = Eigen::VectorXd::Zero(V_2d.rows() * 2);
}

std::vector<int> readVectorIntFile(const std::string &filename) {
    std::vector<int> integers;
    std::ifstream inputFile(filename);

    if (!inputFile.is_open()) {
        std::cerr << "Failed to open the file: " << filename << std::endl;
        return integers;  // Return an empty vector
    }

    int value;
    while (inputFile >> value) {
        integers.push_back(value);
    }

    inputFile.close();
    return integers;
}

Eigen::MatrixXd readFloatsFromFile(const std::string& filename) {
    std::vector<std::pair<float, float>> floatPairs;
    std::ifstream inputFile(filename);

    if (!inputFile.is_open()) {
        std::cerr << "Failed to open the file: " << filename << std::endl;
        return Eigen::MatrixXd();  // Return an empty Eigen matrix
    }

    float value1, value2;
    while (inputFile >> value1 >> value2) {
        floatPairs.push_back(std::make_pair(value1, value2));
    }

    inputFile.close();

    // Convert floatPairs to an Eigen::MatrixXd
    int numPairs = floatPairs.size();
    Eigen::MatrixXd matrix(numPairs, 2);
    for (int i = 0; i < numPairs; ++i) {
        matrix(i, 0) = floatPairs[i].first;
        matrix(i, 1) = floatPairs[i].second;
    }

    return matrix;
}



int main(int argc, char *argv[]){
    std::string path_3d, path_output, path_middle_ids, path_sleeve_ids, path_sleeve_pos;
    std::vector<int> middle_ids, sleeve_ids;
    Eigen::MatrixXd sleeve_pos;
    if (argc >= 2){path_3d = argv[1];}
    if (argc >= 3){path_output = argv[2];}
    if (argc >= 4){
        path_middle_ids = argv[3];
        middle_ids = readVectorIntFile(path_middle_ids);
    }
    if (argc >= 5){
        path_sleeve_ids = argv[4];
        sleeve_ids = readVectorIntFile(path_sleeve_ids);
    }
    if (argc >= 6){
        path_sleeve_pos = argv[5];
        sleeve_pos = readFloatsFromFile(path_sleeve_pos);
    }

    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
    
    Eigen::MatrixXd V_2d, V_3d;
    Eigen::MatrixXi F, init_F;
    igl::readOBJ(path_3d, V_3d, F);

    Eigen::VectorXi bnd;
    igl::boundary_loop(F, bnd);

    std::vector<std::vector<int>> bnds;
    igl::boundary_loop(F, bnds);
    if (bnds.size() != 1) 
        std::cout << "ERROR: wrong topology, number of boundaries = " << bnds.size() << " (should be 1)" << std::endl;
    if (bnds[0].size() != bnd.rows())
        std::cout << "ERROR: wrong boundary size, " << bnds[0].size() << " vs " << bnd.rows() << std::endl;

    if (sleeve_ids.size() > 2) {
        int sv1 = 0;
        int sv2 = 1;
        V_2d = paramLSCMwithConstraint(V_3d, F, sleeve_ids[sv1], sleeve_pos(sv1,0), sleeve_pos(sv1,1), sleeve_ids[sv2], sleeve_pos(sv2,0), sleeve_pos(sv2,1));
    }
    else {
        V_2d = paramLSCM(V_3d, F, bnd); // TODO use the first vertex in sleeves for constraint?
    }
    
    // if LSCM flips the mesh, unflip it
    Eigen::Vector3d v1 = (V_2d.row(F(0, 1)) - V_2d.row(F(0, 0))).transpose();
    Eigen::Vector3d v2 = (V_2d.row(F(0, 2)) - V_2d.row(F(0, 0))).transpose();
    Eigen::Vector3d n0 = v1.cross(v2);
    if (n0(2) < 0) {
        std::cout << "Initial solution flipped, unflipping..." << std::endl;
        V_2d.col(0) *= - 1.0;
    }
    

    Eigen::SparseMatrix<double> A;
    Eigen::VectorXd b;
    DiagonalMatrixXd W;
    Eigen::VectorXd x;

    int n_iterations = 10;
    for (int it=0; it<n_iterations; it++){
        makeSparseMatrix(V_2d, V_3d, F, A, b, W, x, middle_ids, sleeve_ids, sleeve_pos);

        solver.compute(A.transpose() * W * W * A);

        if(solver.info() != Eigen::Success) {
            std::cout << "ERROR: decomposition failed" << std::endl;
        }

        x = solver.solve(A.transpose() * W * W * b);

        if (solver.info() != Eigen::Success) {
            std::cout << "ERROR, solving failed: ";
            if(solver.info() == Eigen::NumericalIssue) 
                std::cout << "NumericalIssue" << std::endl;
            if(solver.info() == Eigen::NoConvergence) 
                std::cout << "NoConvergence" << std::endl;
            if(solver.info() == Eigen::InvalidInput) 
                std::cout << "InvalidInput" << std::endl;
        }
    }

    Eigen::MatrixXd res = Eigen::MatrixXd::Zero(V_2d.rows(), 3);
    for (int i=0; i<x.rows()/2; i++){
        res(i, 0) = x(2 * i); 
        res(i, 1) = x(2 * i + 1);
    }

    igl::writeOBJ(path_output, res, F);
}