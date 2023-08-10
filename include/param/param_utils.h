#pragma once
/**
 * @author Corentin Dumery
 * @brief Some utilities for mesh parameterization and rotation
 * @date 2022-02-28
 */
#include <Eigen/Core>
#include <iostream>


void meshCleanup(Eigen::MatrixXd& V_3d, Eigen::MatrixXi& F);

/**
 * @brief Computes an ideal *reflection* between points1 and points2 using 
 * procrustean analysis.
 * Usage:
 * FROM 2 TO 1 
 * points3t = points2.transpose().colwise() - T_est;
 * points3t = (R_est.transpose() * points3);
 * FROM 1 TO 2 
 * points3t = (R_est * points1.transpose());
 * points3t = points3t.colwise() + T_est;
 * points3 = points3t.transpose(); 
 */
void procrustes(const Eigen::MatrixXd& points1, // to
               const Eigen::MatrixXd& points2, // from
               Eigen::MatrixXd& R_est,
               Eigen::VectorXd& T_est);

/**
 * @brief Make a vetrtex matrix from a triangle
 */
Eigen::MatrixXd makeTriPoints(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, int f_id);
void makeTriPoints(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, int f_id, Eigen::MatrixXd& out);

// Move triangle from 3D (Z != 0) to 2D (Z = 0), with
// same lengths but arbitrary orientation
Eigen::MatrixXd move3Dto2D(const Eigen::MatrixXd& V_tri);

/**
 * @brief Flatten vertex matrix 
 */
Eigen::VectorXd vertices2dToVector(const Eigen::MatrixXd& V);

// Compute barycentric coordinates (u, v, w) for
// point p with respect to triangle (a, b, c)
// credits https://gamedev.stackexchange.com/questions/23743/whats-the-most-efficient-way-to-find-barycentric-coordinates
Eigen::Vector3d barycentricCoords(const Eigen::RowVector3d& p, const Eigen::RowVector3d& a, 
                                  const Eigen::RowVector3d& b, const Eigen::RowVector3d& c);

// Parameterization methods relying on libigl
Eigen::MatrixXd paramARAP(const Eigen::MatrixXd& V_3d, const Eigen::MatrixXi& F, const Eigen::VectorXi& bnd);
Eigen::MatrixXd paramLSCM(const Eigen::MatrixXd& V_3d, const Eigen::MatrixXi& F, const Eigen::VectorXi& bnd);
Eigen::MatrixXd paramSCAF(const Eigen::MatrixXd& V_3d, const Eigen::MatrixXi& F, const Eigen::VectorXi& bnd);


Eigen::MatrixXd paramLSCMwithConstraint(const Eigen::MatrixXd& V_3d, const Eigen::MatrixXi& F,
                                        int v1_id, double v1_u, double v1_v,
                                        int v2_id, double v2_u, double v2_v);

/**
 * @brief Computes the rotation that aligns the first input vector with the second vector 
 */
Eigen::Matrix3d computeRotation(const Eigen::RowVector3d& from,
                                const Eigen::RowVector3d& to);

Eigen::Matrix3d rotationVote(const Eigen::MatrixXd& V_3d,
                             const Eigen::MatrixXd& V_2d,
                             const Eigen::MatrixXi& F,
                             const Eigen::RowVector3d& target_3d,
                             const Eigen::RowVector3d& target_2d);