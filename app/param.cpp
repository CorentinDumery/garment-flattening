
#include <igl/readOBJ.h>

#define ENABLE_PARAM_VIEWER
#ifdef ENABLE_PARAM_VIEWER
#include <igl/opengl/glfw/Viewer.h>
#endif

#include "param/cloth_param.h"

int main(int argc, char *argv[]){

    // example input: ./param ../data/mark_skirt/mark_skirt_back_left_cut.obj

    std::string input_path = "../data/buggy_patch.obj";
    if (argc > 1) input_path = std::string(argv[1]);

    std::vector<int> dart_tip = {296};
    std::vector<std::pair<int, int>> dart_pairs = {{617, 295}, 
                                                   {616, 294},
                                                   {615, 287},
                                                   {613, 286},
                                                   {614, 603}};
    std::vector<std::vector<std::pair<int, int>>> darts = {dart_pairs};

   

    Eigen::MatrixXd V_3d, V_2d;
    Eigen::MatrixXi F;
    igl::readOBJ(input_path, V_3d, F);

    ClothParam cloth(V_3d, F, 0.01, darts, dart_tip);

    //cloth.setDartPairs(darts, dart_tip);

    bool success = cloth.paramAttempt(20);

    std::cout << "Success: " << success << std::endl;
    cloth.printStretchStats();
    Eigen::VectorXd su, sv;
    cloth.getStretchStats(su, sv);
    std::cout << "Stretch start: " << sv.topRows(10) << std::endl;
    V_2d = cloth.getV2d();


     // --- TEST AREA: check UnorderedDart and SimpleDart are equivalent --- //

    /*SimpleDart sd({614, 613, 615, 616, 617, 296, 295, 294, 287, 286, 603});
    UnorderedDart ud(darts[0], dart_tip[0]);

    Eigen::RowVectorXd sd_axis = sd.computeSymmetryAxis(V_2d);
    Eigen::RowVector2d ud_axis = ud.computeSymmetryAxis(V_2d.leftCols(2));
    std::cout << "SD sym axis: " << sd_axis << std::endl;
    std::cout << "UD sym axis: " << ud_axis << std::endl;

    Eigen::MatrixXd sym_sd = sd.getSymmetricPoints(V_2d, sd_axis);
    std::vector<std::pair<Eigen::RowVector2d, Eigen::RowVector2d>> sym_ud = ud.getSymmetricPoints(V_2d.leftCols(2), ud_axis);
    std::cout << "sym_sd: " << sym_sd << std::endl;
    std::cout << "sym_ud: " << std::endl;
    for (auto i: sym_ud){
        std::cout << i.first << std::endl;
        std::cout << "-" << std::endl;
        std::cout << i.second << std::endl;
        std::cout << "--" << std::endl;
    }


    Eigen::MatrixXd targets_sd = sd.computeSymmetryTargets(V_2d, sd_axis);
    std::vector<std::pair<Eigen::RowVector2d, Eigen::RowVector2d>> targets_ud = ud.computeSymmetryTargets(V_2d.leftCols(2), ud_axis);
    std::cout << "targets_sd: " << targets_sd << std::endl;
    std::cout << "targets_ud: " << std::endl;
    for (auto i: targets_ud){
        std::cout << i.first << std::endl;
        std::cout << "-" << std::endl;
        std::cout << i.second << std::endl;
        std::cout << "--" << std::endl;
    }*/


    // --- END TEST AREA --- //

    // ------------- VIEWER ------------- //
    #ifdef ENABLE_PARAM_VIEWER
    bool show_uv = true;

    // Plot the mesh
    igl::opengl::glfw::Viewer viewer;
    viewer.data().set_mesh(V_3d, F);
    viewer.data().set_uv(10 * (V_2d.array() + V_2d.minCoeff()) / (V_2d.maxCoeff() - V_2d.minCoeff()));
    viewer.callback_key_down = [&](igl::opengl::glfw::Viewer& viewer, unsigned char key, int modifier){
        if (key == '1')
            show_uv = false;
        else if (key == '2')
            show_uv = true;

        if (show_uv){
            viewer.data().set_mesh(V_2d, F);
            viewer.core().align_camera_center(V_2d, F);
        }
        else
        {
            viewer.data().set_mesh(V_3d, F);
            viewer.core().align_camera_center(V_3d, F);
        }

        viewer.data().compute_normals();

        return false;
    };

    /*for (int i=0; i<dart_pairs.size(); i++){
        Eigen::RowVector3d colorf = Eigen::RowVector3d::Random(); 
        Eigen::RowVector3d colors = Eigen::RowVector3d::Random(); 
        viewer.data().add_points(V_2d.row(dart_pairs[i].first), colorf);
        viewer.data().add_points(V_2d.row(dart_pairs[i].second), colors);
        viewer.data().add_points(targets_ud[i].first, colorf);
        viewer.data().add_points(targets_ud[i].second, colors);
    }

    viewer.data().add_points(V_2d.row(dart_tip[0]), Eigen::RowVector3d(1.0, 1.0, 1.0));

    Eigen::RowVector2d temp = V_2d.row(dart_tip[0]).leftCols(2);
    Eigen::RowVector3d vvv; 
    vvv(0) = (temp + 10.0 * ud_axis)(0);
    vvv(1) = (temp + 10.0 * ud_axis)(1);
    vvv(2) = 0.0;
    viewer.data().add_edges(V_2d.row(dart_tip[0]), vvv, Eigen::RowVector3d(1.0, 1.0, 1.0));*/

    viewer.data().show_lines = true;
    viewer.data().show_texture = true;
    viewer.launch();
    #endif
}