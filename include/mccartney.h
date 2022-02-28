#pragma once 
/**
 * @file mccartney.h
 * @author Corentin Dumery
 * @brief Measurements of McCartney's energies for anisotropic materials
 * @date 2022-02-28
 */
#include <Eigen/Core>
#include <igl/barycentric_coordinates.h>
#include "param/param_utils.h"

void computeFrameErrors(const Eigen::MatrixXd& V_2di, 
                        const Eigen::MatrixXd& V_3di,
                        double& Esu, double &Esv, double& Er){

    Eigen::MatrixXd V_ref2 = V_2di;
    Eigen::MatrixXd V_def2 = move3Dto2D(V_3di);

    
    for (int i=0; i<3; i++){
        V_ref2(i, 2) = 1.0;
        V_def2(i, 2) = 1.0;
    }

    V_ref2.transposeInPlace();
    V_def2.transposeInPlace();

    Eigen::MatrixXd T = V_def2 * V_ref2.inverse();

    double Su = T.col(0).topRows(2).norm();
    Su = T(0,0)*T(0,0) + T(1,0)*T(1,0);
    double Sv = T.col(1).topRows(2).norm();
    Sv = T(0,1)*T(0,1) + T(1,1)*T(1,1);

    // Nullify translation
    T(0,2) = 0.0;
    T(1,2) = 0.0;

    Eigen::VectorXd v1 = T * Eigen::Vector3d(1.0, 0.0, 0.0);
    Eigen::VectorXd v2 = T * Eigen::Vector3d(0.0, 1.0, 0.0);

    v1(2) = 0;
    v2(2) = 0;

    double angle = std::acos(v1.dot(v2)/(v1.topRows(2).norm() * v2.topRows(2).norm()));
    Esu = Su + 1.0/Su - 2.0;
    Esv = Sv + 1.0/Sv - 2.0;
    Er = std::pow((3.1415/2.0 - angle), 2);
}

Eigen::MatrixXd computeMcCartneyErrors(const Eigen::MatrixXd& V_2di, 
                                       const Eigen::MatrixXd& V_3di,
                                       double& Esu, double &Esv, double& Er){

    if (V_3di.rows() < 3 || V_2di.rows() < 3){
        std::cout << "ERROR: McCartney invalid inputs" << std::endl;
    } 

    Eigen::MatrixXd V_3d(3,3); // We'll move V_3di to align with V_2d
    Eigen::MatrixXd V_2d = V_2di;
    V_2d.row(1) -= V_2d.row(0);
    V_2d.row(2) -= V_2d.row(0);
    V_2d.row(0) -= V_2d.row(0);
    V_3d = move3Dto2D(V_3di);

    if (V_2d.row(0).maxCoeff() > 0 || V_2d.col(2).maxCoeff() > 0){
        std::cout << "ERROR: violated 2d assumptions" << std::endl;
    } 

    // Then, align triangles along weft axis

    Eigen::RowVector3d B = V_2d.row(1);
    Eigen::RowVector3d C = V_2d.row(2);
    Eigen::RowVector3d BC = C - B;

    double d = BC(1) / BC(0); 
    double alpha = - B(1) / (d );
    alpha = - B(1) / BC(1);

    Eigen::RowVector3d X = B + alpha * BC;
    Eigen::RowVector3d Bp = V_3d.row(1);
    Eigen::RowVector3d Cp = V_3d.row(2);
    Eigen::RowVector3d Xp = Bp + alpha * (Cp - Bp);
    Eigen::Matrix3d R = computeRotation(Xp, X);

    V_3d = (R * V_3d.transpose()).transpose();

    double ub = V_2d(1,0);
    double ubp = V_3d(1,0);
    double uc = V_2d(2,0);
    double ucp = V_3d(2,0);

    double vb = V_2d(1,1);
    double vbp = V_3d(1,1);
    double vc = V_2d(2,1);
    double vcp = V_3d(2,1);

    double Su = (vc * ubp - vb * ucp)/(ub * vc - uc * vb);
    double phiv = std::atan((ub * ucp - uc * ubp)/(ub * vcp - uc * vbp));
    double Sv = std::sqrt(
                std::pow(ub * ucp - uc * ubp, 2) + std::pow(ub*vcp - uc * vbp, 2))/(ub * vc - uc * vb);

    // triangle area
    double as = (V_2d.row(0) - V_2d.row(1)).norm();
    double bs = (V_2d.row(2) - V_2d.row(1)).norm();
    double cs = (V_2d.row(0) - V_2d.row(2)).norm();
    double s = (as + bs + cs)/2.0;
    double A_2d = std::sqrt(s * (s - as) * (s - bs) * (s - cs)); // Heron's formula

    double Ksu = 1.0;
    double Ksv = 1.0;
    double Kr = 1.0;
    Esu = 0.5 * A_2d * Ksu * std::pow(Su - 1.0, 2);
    Esv = 0.5 * A_2d * Ksv * std::pow(Sv - 1.0, 2);
    Er = 0.5 * A_2d * Kr * std::pow(phiv, 2);

    bool print_energies = false;
    if (print_energies){
        double total_E = Esu + Esv + Er;
        printf("Stretch U: %f (%f %%)\n", Esu, 100.0*Esu/total_E);
        printf("Stretch V: %f (%f %%)\n", Esv, 100.0*Esv/total_E);
        printf("Strain   : %f (%f %%)\n", Er, 100.0*Er/total_E);
    }

    return V_3d;
}

// Transport matrix "mat" of points in triangle V_2d to triangle V_3d 
Eigen::MatrixXd transportMatrix(const Eigen::MatrixXd& V_2d, const Eigen::MatrixXi& F, 
                                const Eigen::MatrixXd& mat, const Eigen::MatrixXd& V_3d){
    Eigen::MatrixXd res(mat.rows(), 3);
    for (int i=0; i<mat.rows(); i++){
        Eigen::RowVector3d bary;
        Eigen::RowVector3d p = mat.row(i);
        igl::barycentric_coordinates(p, V_2d.row(F(0,0)), V_2d.row(F(0,1)), V_2d.row(F(0,2)), bary);\
        res.row(i) = bary(0) * V_3d.row(F(0,0)) + bary(1) * V_3d.row(F(0,1)) + bary(2) * V_3d.row(F(0,2)); 
    }
    return res;
};