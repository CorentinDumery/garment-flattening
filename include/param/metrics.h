#pragma once
/**
 * @file metrics.h
 * @author Corentin Dumery
 * @brief Static functions measuring quality of a parameterization. 
 * @date 2022-02-28
 */
#include <Eigen/Core>

#include "param/multi_patch_param.h"
#include "param/param_utils.h" 

// TODO this is a duplicate of BaryOptimizer::measureScore
void measureStretchScore(const Eigen::MatrixXd& V_2d, const Eigen::MatrixXd& V_3d, 
                         const Eigen::MatrixXi& F, Eigen::VectorXd& stretch_u_vec,
                         Eigen::VectorXd& stretch_v_vec);

void measureSeamScore(const std::vector<Eigen::MatrixXd>& vec_V_2d,
                      const std::vector<Seam>& seams,
                      double& length_error, double& reflec_error);

void measureAlignmentScore(const Eigen::MatrixXd& V_2d,
                           const Eigen::MatrixXd& V_3d,
                           const Eigen::MatrixXi& F,
                           Eigen::VectorXd& align_error);