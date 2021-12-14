#pragma once
#include <cmath>

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
               Eigen::VectorXd& T_est){

    // https://igl.ethz.ch/projects/ARAP/svd_rot.pdf
    // https://math.stackexchange.com/questions/849217/estimate-rotation-and-translation-from-two-sets-of-points-in-different-coordinat
    // https://en.wikipedia.org/wiki/Procrustes_analysis

    Eigen::MatrixXd points1t = points1.transpose();
    Eigen::MatrixXd points2t = points2.transpose();

    Eigen::VectorXd pb = points1t.rowwise().mean();
    Eigen::VectorXd qb = points2t.rowwise().mean();

    Eigen::MatrixXd X = (points1t.colwise() - pb);
    Eigen::MatrixXd Y = (points2t.colwise() - qb);
    Eigen::MatrixXd S = X * Y.transpose();
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(S, Eigen::ComputeThinU | Eigen::ComputeThinV);
    Eigen::MatrixXd sigma = Eigen::MatrixXd::Identity(svd.matrixU().cols(), svd.matrixV().cols());
    sigma(sigma.rows() - 1, sigma.cols() - 1) = -(svd.matrixV() * svd.matrixU().transpose()).determinant();
    R_est = svd.matrixV() * sigma * svd.matrixU().transpose();
    T_est = qb - R_est * pb;
}

// Move triangle from 3D (Z != 0) to 2D (Z = 0), with
// same lengths but arbitrary orientation
Eigen::MatrixXd move3Dto2D(const Eigen::MatrixXd& V_tri){
    Eigen::MatrixXd V_2d(3,3); // TODO Matrix3d
    // First, move 3D triangle to 2D plane
    // V_2d: put A in (0,0), B in (0, |AB|), and find C 
    double r0 = (V_tri.row(1) - V_tri.row(0)).norm();
    double r1 = (V_tri.row(2) - V_tri.row(0)).norm();
    double r2 = (V_tri.row(2) - V_tri.row(1)).norm();
    V_2d.row(0) = Eigen::RowVector3d(0, 0, 0);
    V_2d.row(1) = Eigen::RowVector3d(r0, 0, 0);
    double CAB_angle = std::acos((r0*r0 + r1*r1 - r2*r2)/(2*r0*r1));
    double l1 = r1 * std::cos(CAB_angle);
    double h = l1 * std::tan(CAB_angle);
    V_2d.row(2) = Eigen::RowVector3d(l1, h, 0);

    r0 = (V_2d.row(1) - V_2d.row(0)).norm();
    r1 = (V_2d.row(2) - V_2d.row(0)).norm();
    r2 = (V_2d.row(2) - V_2d.row(1)).norm();
    if ( std::fabs(r0 - (V_tri.row(1) - V_tri.row(0)).norm()) // just checking...
        +std::fabs(r1 - (V_tri.row(2) - V_tri.row(0)).norm())
        +std::fabs(r2 - (V_tri.row(2) - V_tri.row(1)).norm()) > 0.0001){
        std::cout << "ERROR, flat triangle is different:" << std::endl;
        std::cout << r0 << " vs " << (V_tri.row(1) - V_tri.row(0)).norm() << std::endl;
    }

    return V_2d;
}