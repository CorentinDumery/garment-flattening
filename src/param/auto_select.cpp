#include "param/auto_select.h"

std::vector<int> autoSelect(const Eigen::MatrixXd& V_3d, const Eigen::VectorXi& bnd){
    const int x_axis = 0;
    const int y_axis = 2;
    const int z_axis = 1;

    std::vector<int> best_pair = {-1, -1};
    double best_score = -1;
    for (int i=0; i<bnd.rows(); i++){
        for (int j=0; j<bnd.rows(); j++){
            if (i == j) continue;
            double delta_z = (V_3d(bnd(i), z_axis) - V_3d(bnd(j), z_axis)); // actually no fabs for z, we want to know which one is up
            double total_score = delta_z;
            if (total_score > best_score || best_score == -1){
                best_score = total_score;
                best_pair = {bnd(i), bnd(j)};
            }
        }
    } 

    return best_pair;
}