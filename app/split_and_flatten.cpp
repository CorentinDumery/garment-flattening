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

#include "param/cloth_param.h"
#include "cut_on_plane.h"


#include <igl/doublearea.h>

double findBestCut(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, int axis){
    Eigen::VectorXd tri_area;
    igl::doublearea(V, F, tri_area);
    tri_area = tri_area / 2;
    double total_area = tri_area.sum();

    double lower_bound = V.col(axis).minCoeff();
    double upper_bound = V.col(axis).maxCoeff();

    double epsilon = 0.005 * (upper_bound - lower_bound);
    double z = 0.0;

    while (upper_bound - lower_bound > epsilon) {
        z = (lower_bound + upper_bound) / 2.0;

        double area_below = 0.0;
        double area_above = 0.0;

        for (int i = 0; i < F.rows(); i++) {
            if (V(F(i, 0), axis) < z && V(F(i, 1), axis) < z && V(F(i, 2), axis) < z) {
                area_below += tri_area(i);
            } else if (V(F(i, 0), axis) > z && V(F(i, 1), axis) > z && V(F(i, 2), axis) > z) {
                area_above += tri_area(i);
            } else { // Triangle intersects the plane
            /*
                double z_intercept = (z - V(F(i, 0), axis)) / (V(F(i, 1), axis) - V(F(i, 0), axis));
                Eigen::Vector3d intersection1 = V.row(F(i, 0)) + z_intercept * (V.row(F(i, 1)) - V.row(F(i, 0)));
                tri_area(i) *= (1.0 - z_intercept);
                area_below += tri_area(i) / 2.0;
                area_above += tri_area(i) / 2.0;*/
            }
        }

        if (area_below < area_above) {
            lower_bound = z;
        } else {
            upper_bound = z;
        }
    }

    std::cout << "The z value that divides the mesh into equal areas is: " << z << std::endl;
    return z;
}

void saveEigenVector(std::string savepath, const Eigen::VectorXi& vec) {
    std::ofstream file(savepath);
    if (file.is_open()) {
        Eigen::IOFormat fmt(Eigen::StreamPrecision, Eigen::DontAlignCols, " ", " ");
        file << vec.format(fmt) << std::endl;
        file.close();
    } else {
        std::cerr << "Unable to open the file!" << std::endl;
    }
}

int main(int argc, char *argv[]){

    std::string path_3d, output_dir, cut_axes;

    path_3d = "../data/dress.ply";
    output_dir = "../data/";
    cut_axes = "Z";

    if (argc >= 2){path_3d = argv[1];}
    if (argc >= 3){output_dir = argv[2];}
    if (argc >= 4){cut_axes = argv[3];}

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

    std::vector<double> cut_vals = {findBestCut(V_3d, F, 0), 
                                    findBestCut(V_3d, F, 1), 
                                    findBestCut(V_3d, F, 2)};
    Eigen::VectorXi cut0; 
    Eigen::VectorXi cut1;

    std::vector<Eigen::MatrixXd> V_list = {V_3d};
    std::vector<Eigen::MatrixXi> F_list = {F};
    std::vector<char> axes = {'X', 'Y', 'Z'};
    for (int id=0; id<3; id++){
        char ax = axes[id];
        if (cut_axes.find(ax) != std::string::npos){
            std::cout << "Cutting " << ax << std::endl;
            std::vector<Eigen::MatrixXd> V_list_n;
            std::vector<Eigen::MatrixXi> F_list_n; 
            
            for (int i=0; i<V_list.size(); i++){
                std::vector<Eigen::MatrixXd> V_list_n2;
                std::vector<Eigen::MatrixXi> F_list_n2;
                cutMeshOnPlane(V_list[i], F_list[i], V_list_n2, F_list_n2, cut0, cut1, id, cut_vals[id]);
                V_list_n.insert(V_list_n.end(), V_list_n2.begin(), V_list_n2.end());
                F_list_n.insert(F_list_n.end(), F_list_n2.begin(), F_list_n2.end());
            }

            V_list = V_list_n;
            F_list = F_list_n;
        }
    }
    
    if (cut_axes.size() > 1) std::cout << "WARNING: seam correspondences not (yet) supported for several cuts" << std::endl;
    
    saveEigenVector(output_dir + "/cut_vertex_ids_0.txt", cut0);
    saveEigenVector(output_dir + "/cut_vertex_ids_1.txt", cut1);
    
    for (int mesh=0; mesh<V_list.size(); mesh++){
        // Save 3D connected component for debugging
        std::string file_name_3d = output_dir + "/component3d_" + std::to_string(mesh) + ".obj";
        igl::writeOBJ(file_name_3d, V_list[mesh], F_list[mesh]);
    }
    
    for (int mesh=0; mesh<V_list.size(); mesh++){
        // Compute 2D parameterization
        ClothParam cp(V_list[mesh], F_list[mesh], 0.00);
        for (int i=0; i<5; i++){
            cp.paramAttempt(20);
            cp.printStretchStats();
        }
        V_2d = cp.getV2d();

        // Save 2D
        std::string file_name_2d = output_dir + "/component2d_" + std::to_string(mesh) + ".obj";
        igl::writeOBJ(file_name_2d, V_2d, F_list[mesh]);
    }

    for (int mesh=0; mesh<V_list.size(); mesh++){
        Eigen::VectorXd A;
        igl::doublearea(V_list[mesh], F_list[mesh], A);
        std::cout << "Component area:" << A.sum()  << std::endl;
    }
}