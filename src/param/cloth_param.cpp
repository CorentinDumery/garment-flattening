#include "param/cloth_param.h"
#include <igl/avg_edge_length.h>
#include <nlohmann/json.hpp>
#include <fstream>
#include "param/self_intersect.h"

//#define DEBUG_CLOTH_PARAM
#ifdef DEBUG_CLOTH_PARAM
#include <igl/writeOBJ.h>
#endif

#define CHECK_TOPOLOGY_PARAM
//#define INPUT_CHECKS_PARAM

ClothParam::ClothParam(const Eigen::MatrixXd& V_3d, const Eigen::MatrixXi& F,
                       double max_stretch,
                       const std::vector<std::vector<std::pair<int, int>>>& dart_duplicates,
                       const std::vector<int>& dart_tips,
                       int seam_size,
                       CLOTH_INIT_TYPE init_type,
                       bool disable_stretch_penalty)
            : F_(F), V_3d_(V_3d), max_stretch_(max_stretch), seam_size_(seam_size), init_type_(init_type){
    
    igl::boundary_loop(F, bnd_);

    if (disable_stretch_penalty){
        bo_.enable_stretch_eqs_ = false;
    }

    #ifdef CHECK_TOPOLOGY_PARAM
    std::cout << "Boundary size: " << bnd_.rows() << std::endl;
    std::vector<std::vector<int>> bnds;
    igl::boundary_loop(F, bnds);
    if (bnds.size() != 1) 
        std::cout << "ERROR: wrong topology, number of boundaries = " << bnds.size() << " (should be 1)" << std::endl;
    if (bnds[0].size() != bnd_.rows())
        std::cout << "ERROR: wrong boundary size, " << bnds[0].size() << " vs " << bnd_.rows() << std::endl;
    #endif

    #ifdef INPUT_CHECKS_PARAM
    bool check_degen = false;
    for (int i=0; i<F.rows(); i++){
        if (F(i,0)==F(i,1) || F(i,0)==F(i,2) || F(i,1)==F(i,2)){
            std::cout << "ERROR: degenerate face: " << F.row(i) << std::endl;
            check_degen = true;
        }
    }
    std::cout << "Test degenerate face: " << check_degen << std::endl;

    double tol = 10e-6;
    for (int i=0; i<V_3d.rows(); i++){
        for (int j=i+1; j<V_3d.rows(); j++){
            if ((V_3d.row(i) - V_3d.row(j)).norm() < tol){
                std::cout << "ERROR: duplicate vertices: " << V_3d.row(i) << " vs " << V_3d.row(j) << std::endl;
            }
        }
    }
    #endif

    #ifdef DEBUG_CLOTH_PARAM
    std::cout << "Initializing param..." << std::endl;
    #endif

    if (init_type_ == CLOTH_INIT_ARAP){
        V_2d_ = paramARAP(V_3d_, F_, bnd_);
    }
    else if (init_type_ == CLOTH_INIT_LSCM){
        V_2d_ = paramLSCM(V_3d_, F_, bnd_);
        // if LSCM flips the mesh, unflip it
        Eigen::Vector3d v1 = (V_2d_.row(F_(0, 1)) - V_2d_.row(F_(0, 0))).transpose();
        Eigen::Vector3d v2 = (V_2d_.row(F_(0, 2)) - V_2d_.row(F_(0, 0))).transpose();
        Eigen::Vector3d n0 = v1.cross(v2);
        if (n0(2) < 0) V_2d_.col(1) *= - 1.0;
    }
    else if (init_type_ == CLOTH_INIT_SCAF){
        V_2d_ = paramSCAF(V_3d_, F_, bnd_);
    }

    #ifdef DEBUG_CLOTH_PARAM
    std::cout << "Param initialized." << std::endl;
    #endif

    bool rescale_init = true;
    if (rescale_init){
        V_2d_ *= igl::avg_edge_length(V_3d_, F_) / igl::avg_edge_length(V_2d_, F_);
    }

    setDartPairs(dart_duplicates, dart_tips);
    bo_.setSeamSize(seam_size_);

    #ifdef DEBUG_CLOTH_PARAM
    std::cout << "Allocating memory... (" << F.rows() << "," <<  V_3d.rows() << ")" << std::endl;
    #endif
    bo_.allocateMemory(F.rows(), V_3d.rows());
    
    if (rotate_each_iter_){
        Eigen::MatrixXd R1 = rotationVote(V_3d_, V_2d_, F_, Eigen::RowVector3d(0, 1.0,0.0), 
                                                            Eigen::RowVector3d(0.0,1.0,0.0));
        V_2d_ = (R1 * V_2d_.transpose()).transpose();
    }
    
};

void ClothParam::paramIter(int n_iter){
    for (int i=0; i<n_iter; i++){
        V_2d_ = bo_.localGlobal(V_2d_, V_3d_, F_);

        if (rotate_each_iter_){
            Eigen::MatrixXd R1 = rotationVote(V_3d_, V_2d_, F_, Eigen::RowVector3d(0, 1.0,0.0), 
                                                                Eigen::RowVector3d(0.0,1.0,0.0));
            V_2d_ = (R1 * V_2d_.transpose()).transpose();
        }
    }
}

bool ClothParam::paramAttempt(int max_iter){
    for (int current_iter = 0; current_iter < max_iter; current_iter++){
        bo_.measureScore(V_2d_, V_3d_, F_, stretch_u_, stretch_v_);

        #ifdef DEBUG_CLOTH_PARAM
        if (stretch_u_.maxCoeff() < -0.5){ // not supposed to happen if the initialization is ok
            igl::writeOBJ("../data/buggy/not_good.obj", V_3d_, F_);
            igl::writeOBJ("../data/buggy/not_good_uv.obj", V_2d_, F_);
        }
        if (std::isnan(stretch_u_.maxCoeff())){ // not supposed to happen if the initialization is ok
            igl::writeOBJ("./nanned.obj", V_3d_, F_);
            igl::writeOBJ("./nanned_uv.obj", V_2d_, F_);
        }
        #endif

        /*if (constraintSatisfied()){
            return true;
        }*/
        V_2d_ = bo_.localGlobal(V_2d_, V_3d_, F_);

        if (rotate_each_iter_){
            Eigen::MatrixXd R1 = rotationVote(V_3d_, V_2d_, F_, Eigen::RowVector3d(0, 1.0,0.0), 
                                                                Eigen::RowVector3d(0.0,1.0,0.0));
            V_2d_ = (R1 * V_2d_.transpose()).transpose();
        }
        
        if (enable_intersection_check_){
            if (checkSelfIntersect()){
                // interrupt process early
                // note: we do this after one iteration so we don't punish bad initialization
                return false;
            }
        }
    }
    bo_.measureScore(V_2d_, V_3d_, F_, stretch_u_, stretch_v_);
    return constraintSatisfied() && !checkSelfIntersect();
}

void ClothParam::printStretchStats() const {
    if (stretch_u_.rows() < 1) return;
    std::cout << "Stretch min -> max (avg): " << std::endl;
    printf("U: %f -> %f (%f)\n", stretch_u_.minCoeff(), stretch_u_.maxCoeff(), stretch_u_.mean());
    printf("V: %f -> %f (%f)\n", stretch_v_.minCoeff(), stretch_v_.maxCoeff(), stretch_v_.mean());
};

void vectorSanityCheck(Eigen::VectorXd& vec){
    if (std::isnan(vec.maxCoeff()) ||
        std::isnan(vec.minCoeff())){
        std::cout << "Sanitizing nan..." << std::endl;
        vec = Eigen::VectorXd::Constant(vec.rows(), 10e8);
    }
}

void ClothParam::getStretchStats(Eigen::VectorXd& stretch_u, Eigen::VectorXd& stretch_v) const {
    vectorSanityCheck(stretch_u);
    vectorSanityCheck(stretch_v);
    stretch_u = stretch_u_;
    stretch_v = stretch_v_;
}

void ClothParam::measureStretchStats(const Eigen::MatrixXd& V_2d, const Eigen::MatrixXd& V_3d, 
                                const Eigen::MatrixXi& F, Eigen::VectorXd& stretch_u, 
                                Eigen::VectorXd& stretch_v){
    BaryOptimizer bo;
    bo.measureScore(V_2d, V_3d, F, stretch_u, stretch_v);
    vectorSanityCheck(stretch_u);
    vectorSanityCheck(stretch_v);
}

void ClothParam::setAlignmentVertexPair(int v1_id, int v2_id){
    bo_.setSelectedVertices({v1_id, v2_id});
}

void ClothParam::setDartPairs(const std::vector<std::vector<std::pair<int, int>>>& dart_duplicates,
                    const std::vector<int>& dart_tips){
    bo_.setUnorderedDarts(dart_duplicates, dart_tips);
}

bool ClothParam::checkSelfIntersect() const {
    return selfIntersect(V_2d_, bnd_);
}

void ClothParam::loadConfig(std::string config_path){
    nlohmann::json j;
    std::ifstream i(config_path);
    i >> j;

    double stretch_coeff = j["stretch_coeff"];
    double edges_coeff = j["edges_coeff"];
    double selected_coeff = j["selected_coeff"];
    double tri_align_coeff = j["tri_align_coeff"];
    double dart_sym_coeff = j["dart_sym_coeff"];
    double seam_coeff = j["seam_coeff"];

    bo_.setCoeffs(stretch_coeff, edges_coeff, selected_coeff, 
                  tri_align_coeff, dart_sym_coeff, seam_coeff);
}