
#include <Eigen/Core>
#include <Eigen/LU>
#include <iostream>
#include <cmath>

#include "procustes.h"

double prop1(const Eigen::MatrixXd& V_ref, const Eigen::MatrixXd& V_def){
    Eigen::MatrixXd delta = V_ref - V_def;
    Eigen::RowVectorXd Ad = delta.row(0);
    Eigen::RowVectorXd Bd = delta.row(1);
    Eigen::RowVectorXd Cd = delta.row(2); 

    return std::pow((Bd-Ad)(0), 2) + std::pow((Cd-Ad)(0), 2) + std::pow((Bd-Cd)(0), 2);
}

double prop2(const Eigen::MatrixXd& V_ref, const Eigen::MatrixXd& V_def){
    Eigen::MatrixXd delta = V_ref - V_def;
    Eigen::RowVectorXd Ad = delta.row(0);
    Eigen::RowVectorXd Bd = delta.row(1);
    Eigen::RowVectorXd Cd = delta.row(2); 

    return std::pow((Bd - Ad + Cd - Ad + Cd - Bd)(0), 2);
}

double prop3(const Eigen::MatrixXd& V_ref, const Eigen::MatrixXd& V_def){
    Eigen::MatrixXd delta = V_ref - V_def;
    Eigen::RowVectorXd Ad = delta.row(0);
    Eigen::RowVectorXd Bd = delta.row(1);
    Eigen::RowVectorXd Cd = delta.row(2); 

    return std::pow((Bd - Ad - (Cd - Ad))(0), 2);
}

double prop4(const Eigen::MatrixXd& V_ref, const Eigen::MatrixXd& V_def){
    return prop1(V_ref, V_def) - prop3(V_ref, V_def);
}

//Find best rotation translation, then measure stretch
double prop5(const Eigen::MatrixXd& V_ref, const Eigen::MatrixXd& V_def){
    Eigen::MatrixXd p2_temp(3,3), p2(3,3);
    //p2_temp = makeTriPoints(V_3d, F, f_id);

    Eigen::MatrixXd p1 = V_ref;
    p2_temp = V_def;

    p2 = move3Dto2D(p2_temp);

    Eigen::MatrixXd R_est;
    Eigen::VectorXd T_est;
    procustes(p1, p2, R_est, T_est);

    /*Eigen::MatrixXd p2_r;
    p2_r = p2.transpose();
    p2_r = p2_r.colwise() - T_est;
    p2_r = (R_est.transpose() * p2_r);*/
    //p2_r.transpose();

    //p2_r = p2; // TODO REMOVE !!!!!!!!!!!!!!!

    Eigen::MatrixXd p2_rt, p2_r;
    Eigen::MatrixXd p2t = p2.transpose();
    p2_rt = p2t.colwise() - T_est;
    p2_rt = (R_est.transpose() * p2_rt);
    p2_r = p2_rt.transpose();
    
    //viewer.data().add_points(p1, Eigen::RowVector3d(1.0, 1.0, 0.0));
    //viewer.data().add_points(p2_r, Eigen::RowVector3d(0.0, 1.0, 1.0));

    double ABu = (p1.row(1) - p1.row(0))(0);
    double ApBpu = (p2_r.row(1) - p2_r.row(0))(0);
    double ACu = (p1.row(2) - p1.row(0))(0);
    double ApCpu = (p2_r.row(2) - p2_r.row(0))(0);

    double ABv = (p1.row(1) - p1.row(0))(1);
    double ApBpv = (p2_r.row(1) - p2_r.row(0))(1);
    double ACv = (p1.row(2) - p1.row(0))(1);
    double ApCpv = (p2_r.row(2) - p2_r.row(0))(1);

    double Eu = std::pow(ABu - ApBpu, 2) 
                + std::pow(ACu - ApCpu, 2);
    double Ev = std::pow(ABv - ApBpv, 2) 
                + std::pow(ACv - ApCpv, 2);
    return Eu;
}


double prop6(const Eigen::MatrixXd& V_ref, const Eigen::MatrixXd& V_def){
    //std::cout << V_ref << std::endl;
    //std::cout << V_def << std::endl;
    

    Eigen::MatrixXd V_ref2 = V_ref;
    Eigen::MatrixXd V_def2 = V_def;
    
    for (int i=0; i<3; i++){
        V_ref2(i, 2) = 1.0;
        V_def2(i, 2) = 1.0;
    }

    V_ref2.transposeInPlace();
    V_def2.transposeInPlace();

    Eigen::MatrixXd T = V_def2 * V_ref2.inverse();

    std::cout << std::endl << "Estimated T: " << std::endl << T << std::endl;

    double Su = T.col(0).norm();
    double Sv = T.col(1).norm();
    std::cout << std::endl << "Su: " << Su << std::endl;
    std::cout << "Sv: " << Sv << std::endl;

    Eigen::VectorXd v1 = T * Eigen::Vector3d(1.0, 0.0, 0.0);
    Eigen::VectorXd v2 = T * Eigen::Vector3d(0.0, 1.0, 0.0);

    double angle = std::acos(v1.dot(v2)/(v1.norm() * v2.norm()));
    std::cout << "Angle: " << angle << " (" << angle * 180.0 / 3.1415 <<" deg)"<< std::endl;

    return 0;
}

int main(int argc, char *argv[]){

    // What do we expect from our energies?
    //
    // Stretch:
    // - should be 0 for when there's only shear
    // - be proportional to stretch on the reference's UV axes

    double delta = 0.0;
    Eigen::MatrixXd V_ref(3,3);
    V_ref << 0.0, 0.0, 0.0,
             1.0, 0.0, 0.0,
             1.0, 1.0, 0.0;

    Eigen::MatrixXd V_def(3,3);
    /*V_def << 0.0, 0.0, 0.0,
             1.0, 0.0, 0.0,
             1.0 + delta, 1.0, 0.0; // only adding shear*/

    Eigen::MatrixXd shear(3,3);
    shear << 1.0, 0.1, 0.0,
             0.0, 1.0, 0.0,
             0.0, 0.0, 1.0;

    Eigen::MatrixXd stretch_u(3,3);
    stretch_u << 10.1, 0.0, 0.0,
                 0.0, 1.0, 0.0,
                 0.0, 0.0, 1.0;

    double d = std::sqrt(2.0)/2.0;
    double phi = 0.2;

    Eigen::MatrixXd angle_v(3,3);
    angle_v << 1.0, 0.0, 0.0,
               std::sin(phi), std::cos(phi), 0.0,
               0.0, 0.0, 1.0;

    Eigen::MatrixXd angle_u(3,3);
    angle_u << std::cos(phi), std::sin(phi), 0.0,
               0.0, 1.0, 0.0,
               0.0, 0.0, 1.0;

    angle_u << std::cos(phi), 0.0, 0.0,
               std::sin(phi), 1.0, 0.0,
               0.0, 0.0, 1.0;

    Eigen::MatrixXd transform = stretch_u;
    std::cout << "Transform applied: " << std::endl << transform << std::endl << std::endl; 
    V_def = (transform * V_ref.transpose()).transpose();

    //std::cout << "V_def: " << std::endl << V_def << std::endl << std::endl; 

    std::cout << "prop1: E= " << prop1(V_ref, V_def) << std::endl; 
    std::cout << "prop2: E= " << prop2(V_ref, V_def) << std::endl; 
    std::cout << "prop3: E= " << prop3(V_ref, V_def) << std::endl; 
    std::cout << "prop4: E= " << prop4(V_ref, V_def) << std::endl; 
    std::cout << "prop5: E= " << prop5(V_ref, V_def) << std::endl; 
    std::cout << "prop6: E= " << prop6(V_ref, V_def) << std::endl; 
}
