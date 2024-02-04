#include <Eigen/Core>
#include <Eigen/Sparse>
#include <Eigen/IterativeLinearSolvers> // https://forum.kde.org/viewtopic.php?f=74&t=125165
#include <vector>
#include <fstream>

#include <nlohmann/json.hpp>
#include <igl/readOBJ.h>
#include <igl/writeOBJ.h>
#include <igl/boundary_loop.h>
#include <igl/per_face_normals.h>

#include "param/param_utils.h"

#ifdef FLATTEN_WITH_UI
#include <igl/opengl/glfw/Viewer.h>
#endif

typedef Eigen::DiagonalMatrix<double, Eigen::Dynamic> DiagonalMatrixXd;

bool checkConstraint(const nlohmann::json& constraint) {
    std::string type = constraint["type"];    
    if (type == "equal_values") {
        return constraint.contains("energy_coeff") &&
               constraint.contains("axis") &&
               constraint.contains("opposite_sign") &&
               constraint.contains("ordered_vertex_ids1") &&
               constraint.contains("ordered_vertex_ids2");
    } else if (type == "shared_value") {
        return constraint.contains("axis") &&
               constraint.contains("vertex_ids");
    } else if (type == "predefined_value") {
        return constraint.contains("energy_coeff") &&
               constraint.contains("axis") &&
               constraint.contains("vertex_ids") &&
               constraint.contains("values");
    } else {
        return false; // Unsupported type
    }
}

Eigen::VectorXi readVectorXi(const nlohmann::json& array) {
    Eigen::VectorXi vectorXi(array.size());
    for (size_t i = 0; i < array.size(); ++i) {
        vectorXi[i] = array[i];
    }
    return vectorXi;
}

Eigen::VectorXd readVectorXd(const nlohmann::json& array) {
    Eigen::VectorXd vectorXd(array.size());
    for (size_t i = 0; i < array.size(); ++i) {
        vectorXd[i] = array[i];
    }
    return vectorXd;
}


void makeSparseMatrixFromJson(const Eigen::MatrixXd& V_2d, const Eigen::MatrixXd& V_3d,
                      const Eigen::MatrixXi& F,
                      Eigen::SparseMatrix<double>& A,
                      Eigen::VectorXd& b, 
                      DiagonalMatrixXd& W,
                      Eigen::VectorXd& x,
                      nlohmann::json constraints){

    int next_equation_id = 0;
    
    std::vector<Eigen::Triplet<double>> triplet_list; 
    
    std::vector<double> b_vals;
    std::vector<double> w_vals;

    Eigen::MatrixXd V_tri_2d(3,3);
    Eigen::MatrixXd V_tri_3d(3,3);

    for (int f_id=0; f_id<F.rows(); f_id++) { // The basic ARAP energy
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

    for (auto& [key, value] : constraints.items()) {
        std::string type = value["type"];
        if (type == "equal_values"){
            float energy = value["energy_coeff"];
            int axis = value["axis"];
            int opposite_sign = value["opposite_sign"];
            Eigen::VectorXi v1_ids = readVectorXi(value["ordered_vertex_ids1"]);
            Eigen::VectorXi v2_ids = readVectorXi(value["ordered_vertex_ids2"]);
            if (v1_ids.size() != v2_ids.size()) {
                std::cerr << "ERROR: Sizes of v_ids and v_vals don't match." << std::endl;
            }
            for (int i=1; i<v1_ids.size(); i++){
                triplet_list.push_back(Eigen::Triplet<double>(next_equation_id, 2 * v1_ids[i] + axis, 1.0));
                if (opposite_sign)
                    triplet_list.push_back(Eigen::Triplet<double>(next_equation_id, 2 * v2_ids[i] + axis, 1.0));
                else 
                    triplet_list.push_back(Eigen::Triplet<double>(next_equation_id, 2 * v2_ids[i] + axis, -1.0));
         
                b_vals.push_back(0);
                w_vals.push_back(energy);
                next_equation_id ++;
            }
        }
        else if (type == "shared_value"){
            float energy = value["energy_coeff"];
            int axis = value["axis"];
            Eigen::VectorXi v_ids = readVectorXi(value["vertex_ids"]);
            for (int i=1; i<v_ids.size(); i++){
                // NOTE: it would be better to have the values share a variable instead 
                // of just adding equations
                // Instead: add (vi - v0)^2 = 0 for i in 1,N
                triplet_list.push_back(Eigen::Triplet<double>(next_equation_id, 2 * v_ids[i] + axis, 1.0));
                triplet_list.push_back(Eigen::Triplet<double>(next_equation_id, 2 * v_ids[0] + axis, -1.0));
                b_vals.push_back(0);
                w_vals.push_back(energy);
                next_equation_id ++;
            }
        }
        else if (type == "predefined_value"){
            float energy = value["energy_coeff"];
            int axis = value["axis"];
            Eigen::VectorXi v_ids = readVectorXi(value["vertex_ids"]);
            Eigen::VectorXd v_vals = readVectorXd(value["values"]);
            if (v_ids.size() != v_vals.size()) {
                std::cerr << "ERROR: Sizes of v_ids and v_vals don't match." << std::endl;
            }
            for (int i=0; i<v_ids.size(); i++){
                triplet_list.push_back(Eigen::Triplet<double>(next_equation_id, 2 * v_ids[i] + axis, 1.0));
                b_vals.push_back(v_vals[i]);
                w_vals.push_back(energy);
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

    /*
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
    }*/

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

int main(int argc, char *argv[]){
    std::string path_3d, path_output, path_json, path_sleeve_ids, path_sleeve_pos;
    std::vector<int> middle_ids, sleeve_ids;
    Eigen::MatrixXd sleeve_pos;
    nlohmann::json constraints;
    if (argc >= 4){
        path_3d = argv[1];
        path_output = argv[2];
        path_json = argv[3];
        std::ifstream file(path_json);
    
        if (!file.is_open()) {
            std::cerr << "Failed to open file." << std::endl;
            return 1;
        }
        file >> constraints;

        for (auto& [key, value] : constraints.items()) {
            if (!value.contains("type")) {
                std::cerr << "Error: Missing 'type' field in constraint " << key << std::endl;
                return 1;
            }

            if (!checkConstraint(value)) {
                std::cerr << "Error: Invalid fields in constraint " << key << std::endl;
                return 1;
            }
            std::cout << value << std::endl;
        }

        std::cout << "All constraints are valid." << std::endl;

        //middle_ids = readVectorIntFile(path_middle_ids);
    }
    else {
        std::cerr << "Error: missing arguments" << std::endl;
        return 1;
    }
    
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
    Eigen::MatrixXd V_2d, V_3d;
    Eigen::MatrixXi F, init_F;
    igl::readOBJ(path_3d, V_3d, F);

    // ---------------------   Initial solution for V_2d   ---------------------
    Eigen::VectorXi bnd;
    igl::boundary_loop(F, bnd);

    std::vector<std::vector<int>> bnds;
    igl::boundary_loop(F, bnds);
    if (bnds.size() != 1) 
        std::cout << "ERROR: wrong topology, number of boundaries = " << bnds.size() << " (should be 1)" << std::endl;
    if (bnds[0].size() != bnd.rows())
        std::cout << "ERROR: wrong boundary size, " << bnds[0].size() << " vs " << bnd.rows() << std::endl;

    V_2d = paramLSCM(V_3d, F, bnd); // TODO use the first vertex in sleeves for constraint?

    bool check_degen = true;
    for (int i=0; i<F.rows(); i++){
        if (F(i,0)==F(i,1) || F(i,0)==F(i,2) || F(i,1)==F(i,2)){
            std::cout << "WARNING: degenerate face: " << F.row(i) << std::endl;
            check_degen = true;
        }
    }
    
    // if LSCM flips the mesh, unflip it
    ///*
    Eigen::MatrixXd normals;
    igl::per_face_normals(V_3d, F, normals);
    double mean_normal = normals.col(2).mean();
    Eigen::Vector3d v1 = (V_2d.row(F(0, 1)) - V_2d.row(F(0, 0))).transpose();
    Eigen::Vector3d v2 = (V_2d.row(F(0, 2)) - V_2d.row(F(0, 0))).transpose();
    Eigen::Vector3d n0 = v1.cross(v2);
    bool flip_output = false;
    if (n0(2) * mean_normal < 0) {
        std::cout << "Initial solution flipped, unflipping..." << std::endl;
        V_2d.col(0) *= - 1.0;
        flip_output = true;
    }//*/

    Eigen::MatrixXd R1 = rotationVote(V_3d, V_2d, F, Eigen::RowVector3d(0, 1.0,0.0), 
                                                     Eigen::RowVector3d(0.0,1.0,0.0));
    V_2d = (R1 * V_2d.transpose()).transpose();
    //V_2d = V_2d.rowwise() + (sleeve_pos.row(0) - V_2d.row(sleeve_ids[0]));
    std::cout << "Initial solution computed. " << std::endl;

    // ---------------------   Optimization   ---------------------
    Eigen::SparseMatrix<double> A;
    Eigen::VectorXd b;
    DiagonalMatrixXd W;
    Eigen::VectorXd x;

    #ifdef FLATTEN_WITH_UI
    std::vector<Eigen::MatrixXd> V_2d_list = {V_2d};
    #endif

    int n_iterations = 100;
    for (int it=0; it<n_iterations; it++){
        makeSparseMatrixFromJson(V_2d, V_3d, F, A, b, W, x, constraints);

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

        /*Eigen::MatrixXd res = Eigen::MatrixXd::Zero(V_2d.rows(), 3);
        for (int i=0; i<x.rows()/2; i++){
            res(i, 0) = x(2 * i); 
            res(i, 1) = x(2 * i + 1);
        }*/

        for (int i=0; i<x.rows()/2; i++){
            V_2d(i, 0) = x(2 * i); 
            V_2d(i, 1) = x(2 * i + 1);
        }

        #ifdef FLATTEN_WITH_UI
        V_2d_list.push_back(V_2d);
        #endif
    }

    //if (flip_output){
    //    V_2d.col(0) *= - 1.0;
    //}

    igl::writeOBJ(path_output, V_2d, F);

    #ifdef FLATTEN_WITH_UI
    igl::opengl::glfw::Viewer viewer;

    int current = 0;

    auto updateMesh = [&]() {
        viewer.data(0).clear();
        viewer.data(0).set_mesh(V_2d_list[current], F);
        return true;
    };

    // Bind the updateMesh function to the viewer's 'g' key
    viewer.callback_key_pressed = [&](igl::opengl::glfw::Viewer&, unsigned int key, int) -> bool {
        if (key == 'g' || key == 'G') {
            current = (current + 1) % V_2d_list.size();
            std::cout << "V_2d iteration: " << current << std::endl;
            updateMesh();
        }
        return true;
    };

    // Initialize the viewer with the first mesh
    updateMesh();


    //viewer.data().set_mesh(V_2d, F);
    viewer.append_mesh();
    V_3d = V_3d.rowwise() + Eigen::RowVector3d(1.0, 0, 0);
    viewer.data_list[1].set_mesh(V_3d, F);
    viewer.launch();
    #endif
}