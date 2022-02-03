#pragma once
#include <Eigen/Core>
#include <vector>

Eigen::MatrixXi cutsMatrix(const std::vector<int>& selected, // TODO put somewhere else
                           const Eigen::MatrixXi& F){
    Eigen::MatrixXi cuts = Eigen::MatrixXi::Zero(F.rows(), 3);

    for (int i=0; i<F.rows(); i++){
        for (int j=0; j<selected.size(); j++){
            for (int k=0; k<selected.size(); k++){
                if (k==j) continue;
                for (int l=0; l<3; l++){
                    if (selected[j] == F(i,l)){
                        int m = (l+1)%3;
                        if (selected[k] == F(i,m)){
                            cuts(i, l) = 1;
                        } 
                    }
                }
            }
        }
    }
    return cuts;
}

// Example:
/*
Eigen::MatrixXd V;
Eigen::MatrixXi F;
Eigen::VectorXi I; // vertex correspondence pre/post cut
igl::readOBJ(..., V, F);

std::vector<int> selec = {...};
Eigen::MatrixXi cuts = cutsMatrix(selec, F);
igl::cut_mesh(V, F, cuts, I);
igl::writeOBJ("../../parafashion/data/semi_skirt_test/semi_skirt_quarter.obj", V, F);*/