#pragma once

/**
 * @file multi_patch_param.h
 * @author Corentin Dumery
 * @brief Functions that  
 * @date 2022-02-28
 * 
 */

#include <vector>
#include <memory>
#include <nlohmann/json.hpp>
#include <igl/writeOBJ.h>
#include <fstream>
#include "param/cloth_param.h"

/**
 * @brief Representation of a seam, i.e. a full cut
 * that will be sewn in the final garment.
 * 
 */
struct Seam {
    int patch1_id;
    int patch2_id;
    std::vector<std::pair<int, int>> corres;
    // where corres[i].first is in patch1_id
    //       corres[i].second is in patch2_id
    
    // (ordering of pairs not necessary)
    // (it is possible that patch1_id == patch2_id)
};

/**
 * @brief Final parameterization, where cuts and patches are assumed to be final.
 * This allows us to define seams and enforce seam reflectability as defined
 * in our paper. 
 * 
 */
bool finalParamMultiPatch(const std::vector<Eigen::MatrixXd>& vec_V_3d, 
                          const std::vector<Eigen::MatrixXi>& vec_F,
                          const std::vector<std::vector<std::vector<std::pair<int, int>>>>& vec_dart_duplicates,
                          const std::vector<std::vector<int>>& vec_dart_tips,
                          const std::vector<Seam>& seams,
                          std::vector<Eigen::MatrixXd>& vec_V_2d,
                          CLOTH_INIT_TYPE init_type = CLOTH_INIT_LSCM);

/**
 * @brief Computes target position for seam vertices,
 * following seam reflectability rules 
 *  
 */
void computeTargetPositions(const Eigen::MatrixXd& V1,
                            const Eigen::MatrixXd& V2,
                            const Seam& seam,
                            Eigen::MatrixXd& targets_p,
                            Eigen::MatrixXd& targets_q);