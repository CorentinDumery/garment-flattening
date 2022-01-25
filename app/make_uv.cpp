

#include <igl/readOBJ.h>
#include <igl/writeOBJ.h>

void join(std::string path_3d, std::string path_2d, std::string path_out){
    Eigen::MatrixXd V_3d, V_uv;
    Eigen::MatrixXi F;

    igl::readOBJ(path_3d, V_3d, F);
    igl::readOBJ(path_2d, V_uv, F);

    Eigen::MatrixXd CN;
    Eigen::MatrixXi FN, FTC;


    //   V  #V by 3 mesh vertex positions
    //   F  #F by 3|4 mesh indices into V
    //   CN #CN by 3 normal vectors
    //   FN  #F by 3|4 corner normal indices into CN
    //   TC  #TC by 2|3 texture coordinates
    //   FTC #F by 3|4 corner texture coord indices into TC

    //   igl::writeOBJ(path_out, V_3d, F, CN, FN, V_uv, FTC);

    Eigen::MatrixXd corner_normals;
    Eigen::MatrixXi fNormIndices;

    Eigen::MatrixXd UV_V = V_uv;
    Eigen::MatrixXi UV_F = F;

    igl::writeOBJ(path_out,
        V_3d,
        F,
        corner_normals, fNormIndices, UV_V, UV_F);
}

int main(int argc, char *argv[]){

    std::string path_3d = "../data/semisphere.obj";
    std::string path_2d = "../saved_uv.obj";
    std::string path_out = "../join_uv.obj";

    std::vector<std::string> paths_3d = {};
    std::vector<std::string> paths_2d = {};
    std::vector<std::string> paths_out = {};

    std::string base_name = "../../parafashion/data/perfect/attempt3/perfect.objpatch_";
    for (int i=0; i<7; i++){
        paths_3d.push_back(base_name + "3D_" + std::to_string(i) + ".obj");
        paths_2d.push_back(base_name + "UV_" + std::to_string(i) + ".obj");
        paths_out.push_back(base_name + "joined_" + std::to_string(i) + ".obj");

        join(paths_3d[i], paths_2d[i], paths_out[i]);
    }

    
}