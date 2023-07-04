/**
 * @author Corentin Dumery
 * @brief Simple executable that flattens a mesh, optionally using 
 * some of the garment-flattening features.
 * @date 2023-07-01
 * 
 */

#include <igl/readOBJ.h>
#include <igl/doublearea.h>

#include "param/metrics.h"
#include "param/cloth_param.h"

int main(int argc, char *argv[]){

    std::string path_3d, path_output;

    if (argc >= 2){path_3d = argv[1];}
    if (argc >= 3){path_output = argv[2];}

    Eigen::MatrixXd V_2d, V_3d;
    Eigen::MatrixXi F;
    igl::readOBJ(path_3d, V_3d, F);

    ClothParam cp(V_3d, F, 0.00);
    for (int i=0; i<5; i++){
        cp.paramAttempt(10);
        cp.printStretchStats();
        std::cout << "Self intersect: " << cp.checkSelfIntersect() << std::endl;
    }
    V_2d = cp.getV2d();

    igl::writeOBJ(path_output, V_2d, F);
}