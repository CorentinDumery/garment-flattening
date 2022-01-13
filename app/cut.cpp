

#include <igl/readOBJ.h>
#include <igl/writeOBJ.h>
#include <igl/cut_mesh.h>
#include <iostream>

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

int main(int argc, char *argv[]){

    Eigen::MatrixXd V0, V1, V2, V3;
    Eigen::MatrixXi F0, F1, F2, F3;
    Eigen::VectorXi I1, I2, I3; // vertex correspondence pre/post cut
    igl::readOBJ("../../parafashion/data/semi_skirt_test/semi_skirt_uncut.obj", V0, F0);
    std::cout << "F0 " << F0.rows() << std::endl;
    V1 = V0;
    V2 = V0;
    V3 = V0;
    F1 = F0;
    F2 = F0;
    F3 = F0;


    std::vector<int> selec_quarter = {116, 123, 134, 136, 138, 147, 148, 150, 244};
    std::vector<int> selec_semi = {70, 85, 86, 116, 119, 122, 123, 126, 128, 129, 130, 133, 134, 136, 138, 147, 148, 150, 244};
    std::vector<int> selec_three_quaters = {59, 60, 64, 67, 68, 70, 72, 78, 85, 86, 116, 119, 122, 123, 126, 128, 129, 130, 133, 134, 136, 138, 147, 148, 150, 160, 161, 216, 244};

    Eigen::MatrixXi cuts1 = cutsMatrix(selec_quarter, F1);
    Eigen::MatrixXi cuts2 = cutsMatrix(selec_semi, F2);
    Eigen::MatrixXi cuts3 = cutsMatrix(selec_three_quaters, F3);

    std::cout << "V1 " << V1.rows() << std::endl;
    std::cout << "F0 " << F0.rows() << std::endl;
    std::cout << "cuts1 " << cuts1.rows() << std::endl;
    igl::cut_mesh(V1, F1, cuts1, I1);
    igl::cut_mesh(V2, F2, cuts2, I2);
    igl::cut_mesh(V3, F3, cuts3, I3);

    igl::writeOBJ("../../parafashion/data/semi_skirt_test/semi_skirt_quarter.obj", V1, F1);
    igl::writeOBJ("../../parafashion/data/semi_skirt_test/semi_skirt_semi.obj", V2, F2);
    igl::writeOBJ("../../parafashion/data/semi_skirt_test/semi_skirt_three_quaters.obj", V3, F3);

    
}