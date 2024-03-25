#include <iostream>
#include <igl/readOBJ.h>
#include <igl/readPLY.h>
#include <igl/writeOBJ.h>
#include <igl/opengl/glfw/Viewer.h>

#include "split_ccs.h"


int main(int argc, char *argv[]) {
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " <input_mesh.obj>" << std::endl;
        return 1;
    }

    std::string path_3d = argv[1];

    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    std::cout << "Reading file..." << std::endl;
    if (path_3d.substr(path_3d.length() - 4) == ".ply"){
        igl::readPLY(path_3d, V, F);
    }
    else {
        igl::readOBJ(path_3d, V, F);
    }

    std::vector<Eigen::MatrixXd> V_comps;
    std::vector<Eigen::MatrixXi> F_comps;
    Eigen::VectorXi vertex_components;
    Eigen::VectorXi face_components;
    std::vector<Eigen::VectorXi> v_maps;
    std::vector<Eigen::VectorXi> f_maps;

    splitMeshIntoCCs(V, F, V_comps, F_comps, vertex_components,
                       face_components, v_maps, f_maps);

    Eigen::MatrixXd V_merged;
    Eigen::MatrixXi F_merged;
    mergeCCsBack(V_comps, vertex_components, v_maps, V.rows(), V_merged);
    F_merged = F;

    igl::opengl::glfw::Viewer viewer;
    //viewer.data(0).set_mesh(V_comps[0], F_comps[0]);
    viewer.data(0).set_mesh(V_merged, F_merged);
    viewer.launch();

    for (int i=0; i<V_comps.size(); i++){
        std::string output_path = "component_" + std::to_string(i) + ".obj";
        if (!igl::writeOBJ(output_path, V_comps[i], F_comps[i])) {
            std::cerr << "Failed to write component " << i << " to OBJ file." << std::endl;
            return 1;
        }
        std::cout << "Component " << i << " saved to " << output_path << std::endl;
    }

    return 0;
}
