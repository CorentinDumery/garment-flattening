

#include <igl/readOBJ.h>
#include <igl/doublearea.h>

#include "param/metrics.h"

int main(int argc, char *argv[]){

    std::string path_2d = "/media/corentin/A294-702A/fig2/redoFIG2/mark_skirt_uncut.objpatch_UV_0.obj";
    std::string path_3d = "/media/corentin/A294-702A/fig2/redoFIG2/mark_skirt_uncut.objpatch_3D_0.obj";

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

    /*for (int i=0; i<align_error.rows(); i++){
        std::cout << align_error(i) << ", ";
    }
    std::cout << std::endl;*/

    Eigen::VectorXd A;
    igl::doublearea(V_3d, F, A);

    for (int i=0; i<A.rows(); i++){
        std::cout << A(i)/2.0 << ", ";
    }
    std::cout << std::endl;
    std::cout << "Total area: " << A.sum()/2.0 << std::endl;
}