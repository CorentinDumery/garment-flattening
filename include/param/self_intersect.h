#pragma once
/**
 * @file self_intersect.h
 * @author Corentin Dumery
 * @brief Detect mesh self intersections 
 * @date 2022-02-28
 * 
 */
#include <Eigen/Geometry>

/**
 * @brief Checks for self intersection on the border edges of a 2D mesh 
 */
bool selfIntersect(const Eigen::MatrixXd& V_2d, const Eigen::VectorXi& bnd);