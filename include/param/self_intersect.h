#pragma once

#include <Eigen/Geometry>

bool segmentIntersect(const Eigen::RowVector2d& A,
                      const Eigen::RowVector2d& B,
                      const Eigen::RowVector2d& C,
                      const Eigen::RowVector2d& D);

bool selfIntersect2D(const Eigen::MatrixXd& V_2d, const Eigen::VectorXi& bnd);

bool selfIntersect(const Eigen::MatrixXd& V_2d, const Eigen::VectorXi& bnd);