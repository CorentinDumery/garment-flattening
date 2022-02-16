/**
 * @author Corentin Dumery
 * @brief Measures alignment between parameterization and desired alignment axis. Also 
 * see scripts/angle_deviation_viz.ipynb for visualization.
 * @date 2022-02-04
 * 
 */

#include <igl/readOBJ.h>
#include <igl/doublearea.h>

#include "param/metrics.h"

int main(int argc, char *argv[]){

    std::string path_2d = "../data/skirt/skirt_UV_0.obj";
    std::string path_3d = "../data/skirt/skirt_3D_0.obj";

    if (argc >= 3){
        path_2d = argv[1];
        path_3d = argv[2];
    }

    Eigen::MatrixXd V_2db, V_2d, V_3d;
    Eigen::MatrixXi F;
    igl::readOBJ(path_2d, V_2db, F);
    igl::readOBJ(path_3d, V_3d, F);

    V_2d = Eigen::MatrixXd::Zero(V_2db.rows(), 3);
    V_2d.col(0) = V_2db.col(0);
    V_2d.col(1) = V_2db.col(1);

    Eigen::VectorXd align_error;
    measureAlignmentScore(V_2d, V_3d, F, align_error);

    align_error = align_error.array();

    Eigen::VectorXd A;
    igl::doublearea(V_3d, F, A);

    for (int i=0; i<A.rows(); i++){
        std::cout << A(i)/2.0 << ", ";
    }
    std::cout << std::endl;
    std::cout << "Total area: " << A.sum()/2.0 << std::endl;
}