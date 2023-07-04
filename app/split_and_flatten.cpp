/**
 * @author Corentin Dumery
 * @brief Simple executable that splits a mesh along the XZ plane, 
 * and flattens each side, optionally using some of the garment-flattening features.
 * @date 2023-07-04
 * 
 */

#include <igl/readOBJ.h>
#include <igl/writeOBJ.h>
#include <igl/readPLY.h>
#include <igl/doublearea.h>

#include "param/cloth_param.h"
#include "cut_on_plane.h"

int main(int argc, char *argv[]){

    std::string path_3d;

    path_3d = "../data/dress.ply";

    if (argc >= 2){path_3d = argv[1];}
    if (argc >= 3){}

    Eigen::MatrixXd V_2d, V_3d;
    Eigen::MatrixXi F;
    if (path_3d.substr(path_3d.length() - 3) == "obj")
        igl::readOBJ(path_3d, V_3d, F);
    else if (path_3d.substr(path_3d.length() - 3) == "ply")
        igl::readPLY(path_3d, V_3d, F);
    else {
        std::cout << "Format not recognized" << std::endl;
        return -1;
    }

    std::vector<Eigen::MatrixXd> V_list;
    std::vector<Eigen::MatrixXi> F_list;
    cutMeshOnPlane(V_3d, F, V_list, F_list);
    
    for (int mesh=0; mesh<V_list.size(); mesh++){
        ClothParam cp(V_list[mesh], F_list[mesh], 0.00);
        for (int i=0; i<5; i++){
            cp.paramAttempt(10);
            cp.printStretchStats();
        }
        V_2d = cp.getV2d();

        std::string file_name = "../data/component" + std::to_string(mesh) + ".obj";
        igl::writeOBJ(file_name, V_2d, F_list[mesh]);
    }
}