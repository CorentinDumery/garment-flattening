/**
 * @author Corentin Dumery
 * @brief Interface for woven textile parameterization.
 * @date 2022-02-04
 * 
 */

#pragma once

#include <Eigen/Core>
#include <igl/boundary_loop.h>
#include <memory>

#include "param/bary_optimizer.h"
#include "param/param_utils.h"

enum CLOTH_INIT_TYPE {CLOTH_INIT_ARAP, CLOTH_INIT_LSCM, CLOTH_INIT_SCAF};

class ClothParam {
public:

    // ----- MAIN FUNCTIONS ----- //

    // Initialize and allocate memory.
    // max_stretch: target stretch
    // Dart and seam information is optional.
    ClothParam(const Eigen::MatrixXd& V_3d, const Eigen::MatrixXi& F,
               double max_stretch = 0.05,
               const std::vector<std::vector<std::pair<int, int>>>& dart_duplicates = {},
               const std::vector<int>& dart_tips = {},
               int seam_size = 0,
               CLOTH_INIT_TYPE = CLOTH_INIT_LSCM, 
               bool disable_stretch_penalty = false);

    // Tries to parameterize. Returns whether max_stretch 
    // was satisfied within max number of iterations,
    // AND solution doesn't self intersect.
    bool paramAttempt(int max_iter = 10);

    // ----- OPTIONAL ----- //

    // Quick boundary self intersection test
    bool checkSelfIntersect() const;

    // Define seam targets as explained in paper
    void setSeamTargets(const std::vector<Eigen::MatrixXd>& targets_p,
                        const std::vector<Eigen::VectorXi>& p_ids){
        bo_.setSeamTargets(targets_p, p_ids);
    }

    // Tries to align V axis in parameterization with
    // direction v1 -> v2 
    // (only if enable_selected_eqs_ is enabled)
    void setAlignmentVertexPair(int v1_id, int v2_id);

    // ----- MISC. ----- //

    // Runs param for specified number of iterations, without stopping
    void paramIter(int n_iter);

    void printStretchStats() const;
    void getStretchStats(Eigen::VectorXd& stretch_u, Eigen::VectorXd& stretch_v) const;

    static void measureStretchStats(const Eigen::MatrixXd& V_2d, const Eigen::MatrixXd& V_3d, 
                                    const Eigen::MatrixXi& F, Eigen::VectorXd& stretch_u, 
                                    Eigen::VectorXd& stretch_v);

    Eigen::MatrixXd getV2d(){return V_2d_;}
    void setV2d(const Eigen::MatrixXd& V_2d){V_2d_ = V_2d;}

    bool constraintSatisfied() const {
        if (stretch_u_.rows() < 1) return false;
        return stretch_u_.minCoeff() > -max_stretch_ 
            && stretch_u_.maxCoeff() < max_stretch_ 
            && stretch_v_.minCoeff() > -max_stretch_ 
            && stretch_v_.maxCoeff() < max_stretch_;
    };

    void disableIntersectionCheck(){enable_intersection_check_ = false;};
    void enableIntersectionCheck(){enable_intersection_check_ = true;};

    void loadConfig(std::string config_path);

    Eigen::VectorXi getBnd(){return bnd_;}

    void setCoeffs(float stretch_f, float edges_f){
        bo_.stretch_coeff_ = stretch_f;
        bo_.edges_coeff_ = edges_f;
    }

    // -- Multiple poses -- //
    // We share b vectors across several instances to handle multiple poses
    Eigen::VectorXd getB() {return bo_.getB();};
    void stealSubBs(std::vector<std::unique_ptr<ClothParam>>& subs){
        std::vector<Eigen::VectorXd> other_bs;
        for (int i=0; i<subs.size(); i++){
            other_bs.push_back(subs[i]->getB());
        }
        bo_.other_bs = other_bs;
    };
    void sendMainV2d(std::vector<std::unique_ptr<ClothParam> >& subs){
        for (int i=0; i<subs.size(); i++){
            subs[i]->setV2d(V_2d_);
        }
    };

private:
    // Set during init
    const Eigen::MatrixXd V_3d_;
    const Eigen::MatrixXi F_;
    const double max_stretch_;
    int seam_size_;
    CLOTH_INIT_TYPE init_type_;

    Eigen::VectorXi bnd_;
    Eigen::MatrixXd V_2d_;
    BaryOptimizer bo_;

    Eigen::VectorXd stretch_u_, stretch_v_;

    // config
    bool enable_intersection_check_ = true;
    bool rotate_each_iter_ = true;

    // Sets pairs of vertex ids which should be symmetrical in a given dart.
    // first indexing: dart id
    // second indexing: pair id in dart    
    // pair: two vertex ids, which should be symmetrical w.r.t. dart
    // NOTE: needs to be called before memory allocation!
    void setDartPairs(const std::vector<std::vector<std::pair<int, int>>>& dart_duplicates,
                      const std::vector<int>& dart_tips);
};