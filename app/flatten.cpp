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
#include "param/param_utils.h"

int main(int argc, char *argv[]){

    std::string path_3d, path_output;

    if (argc >= 2){path_3d = argv[1];}
    if (argc >= 3){path_output = argv[2];}

    Eigen::MatrixXd V_2d, V_3d;
    Eigen::MatrixXi F, init_F;
    igl::readOBJ(path_3d, V_3d, F);

    init_F = F;

    meshCleanup(V_3d, F);

    Eigen::VectorXi bnd;
    igl::boundary_loop(F, bnd);

    enum PARAM_TYPE {PARAM_CLOTH, PARAM_ARAP, PARAM_SCAF};
    PARAM_TYPE param_type = PARAM_CLOTH;

    if (param_type == PARAM_ARAP){ // Faulty libigl implementation?
        V_2d = paramARAP(V_3d, F, bnd);
    }
    else if (param_type == PARAM_SCAF){
        V_2d = paramSCAF(V_3d, F, bnd);
    }
    else if (param_type == PARAM_CLOTH){
        float stretch_f = 0.0;
        float edges_f = 1.0;
        ClothParam cp(V_3d, F, 0.00, {}, {}, 0, CLOTH_INIT_LSCM, true);
        cp.setCoeffs(stretch_f, edges_f);
        for (int i=0; i<5; i++){
            cp.paramAttempt(10);
            cp.printStretchStats();
            std::cout << "Self intersect: " << cp.checkSelfIntersect() << std::endl;
        }
        V_2d = cp.getV2d();
    }

    igl::writeOBJ(path_output, V_2d, F);

    if (F.rows() != init_F.rows()) std::cout << "WARNING: output F has different number of rows from input" << std::endl;
}