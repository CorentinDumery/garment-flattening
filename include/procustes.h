#pragma once
#include <Eigen/Core>

/* FROM 2 TO 1 
line3 = line2t.colwise() - T_est;
line3 = (R_est.transpose() * line3);
*/

/* FROM 1 TO 2 
line3 = (R_est * line3);
line3 = line3.colwise() + T_est;
*/

void procustes(const Eigen::MatrixXd& points1, // to
               const Eigen::MatrixXd& points2, // from
               Eigen::MatrixXd& R_est,
               Eigen::VectorXd& T_est);

// Move triangle from 3D (Z != 0) to 2D (Z = 0), with
// same lengths but arbitrary orientation
Eigen::MatrixXd move3Dto2D(const Eigen::MatrixXd& V_tri);

// TODO Matrix3d ?
Eigen::MatrixXd makeTriPoints(const Eigen::MatrixXd V, const Eigen::MatrixXi F, int f_id);