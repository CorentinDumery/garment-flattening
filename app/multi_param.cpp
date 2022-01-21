
#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/readOBJ.h>

#include "param/multi_patch_param.h"
#include "param/metrics.h"

int main(int argc, char *argv[]){

    int n_patches = 3;

    Eigen::MatrixXd V1, V2, V3, V1_out, V2_out, V3_out;
    Eigen::MatrixXi F1, F2, F3;
    std::vector<Eigen::MatrixXd> out_vec;

    igl::readOBJ("../data/patches/patch_3D_1.obj", V1, F1);
    igl::readOBJ("../data/patches/patch_3D_3.obj", V2, F2);
    igl::readOBJ("../data/patches/patch_3D_4.obj", V3, F3);
    //igl::readOBJ("../data/seam_test/left_piece.obj", V1, F1);
    //igl::readOBJ("../data/seam_test/right_piece.obj", V2, F2);
    /*igl::readOBJ("../data/seam_test_bis/top_piece.obj", V1, F1);
    igl::readOBJ("../data/seam_test_bis/middle_piece.obj", V2, F2);
    igl::readOBJ("../data/seam_test_bis/right_piece.obj", V3, F3);*/


    /*Seam s; // for seam test
    s.patch1_id = 0;
    s.patch2_id = 1;
    s.corres.push_back(std::make_pair(14, 13));
    s.corres.push_back(std::make_pair(15, 16));
    s.corres.push_back(std::make_pair(6, 10));
    s.corres.push_back(std::make_pair(19, 20));
    s.corres.push_back(std::make_pair(11, 12));
    s.corres.push_back(std::make_pair(18, 17));
    s.corres.push_back(std::make_pair(7, 4));*/

    /*Seam s1;
    s1.patch1_id = 0;
    s1.patch2_id = 1;
    s1.corres.push_back(std::make_pair(33, 16));
    s1.corres.push_back(std::make_pair(13, 2));
    s1.corres.push_back(std::make_pair(31, 11));
    s1.corres.push_back(std::make_pair(21, 5));
    s1.corres.push_back(std::make_pair(32, 17));
    s1.corres.push_back(std::make_pair(10, 0));

    Seam s2;
    s2.patch1_id = 1;
    s2.patch2_id = 2;
    s2.corres.push_back(std::make_pair(0, 16));
    s2.corres.push_back(std::make_pair(9, 34));
    s2.corres.push_back(std::make_pair(4, 22));
    s2.corres.push_back(std::make_pair(15, 30));
    s2.corres.push_back(std::make_pair(1, 19));
    s2.corres.push_back(std::make_pair(10, 33));*/

    if (n_patches == 2){
        finalParamMultiPatch({V1, V2}, {F1, F2}, 
                          {{}, {}}, // dart dupl
                          {{}, {}}, // dart tips
                          {}, // seams
                          out_vec);
    }

    if (n_patches == 3){
        finalParamMultiPatch({V1, V2, V3}, {F1, F2, F3}, 
                          {{}, {}, {}}, // dart dupl
                          {{}, {}, {}}, // dart tips
                          {}, // seams
                          out_vec);
    }

    
    
    V1_out = out_vec[0];
    V2_out = out_vec[1];
    if (n_patches >=3)
        V3_out = out_vec[2];

    Eigen::MatrixXd R1 = rotationVote(V1, V1_out, F1, Eigen::RowVector3d(0,1.0,0.0), Eigen::RowVector3d(0.0,1.0,0.0));
    Eigen::MatrixXd R2 = rotationVote(V2, V2_out, F2, Eigen::RowVector3d(0,1.0,0.0), Eigen::RowVector3d(0.0,1.0,0.0));
    Eigen::MatrixXd R3 = rotationVote(V3, V3_out, F3, Eigen::RowVector3d(0,1.0,0.0), Eigen::RowVector3d(0.0,1.0,0.0));


    V1_out = (R1 * V1_out.transpose()).transpose();
    V2_out = (R2 * V2_out.transpose()).transpose();
    V3_out = (R3 * V3_out.transpose()).transpose();


    V1_out.col(0) = V1_out.col(0).array() + 50.0;
    V3_out.col(0) = V3_out.col(0).array() - 50.0;

    // --- VISUALIZATION ---

    igl::opengl::glfw::Viewer viewer;
    viewer.append_mesh();
    viewer.append_mesh();
    int mesh1_id = viewer.data_list[0].id;
    int mesh2_id = viewer.data_list[1].id;
    int mesh3_id = viewer.data_list[2].id;

    igl::opengl::glfw::imgui::ImGuiMenu menu;
    menu.callback_draw_viewer_window = []() {};
    viewer.plugins.push_back(&menu);

    auto updateViz = [&](){
        viewer.data().clear_edges();
        viewer.data().clear_points();
        
    };

    //helper function for menu
    auto make_checkbox = [&](const char *label, unsigned int &option) {
        return ImGui::Checkbox(
            label,
            [&]() { return viewer.core().is_set(option); },
            [&](bool value) { return viewer.core().set(option, value); });
    };

    menu.callback_draw_custom_window = [&]() {
        ImGui::SetNextWindowPos(ImVec2(10, 10), ImGuiCond_FirstUseEver);
        ImGui::SetNextWindowSize(ImVec2(350, -1), ImGuiCond_FirstUseEver);
        if (ImGui::Begin("Multiparam")) {
            
            
            ImGui::End();
        }
    
    };

    viewer.data(mesh1_id).set_mesh(V1_out, F1);
    viewer.data(mesh2_id).set_mesh(V2_out, F2);
    if (n_patches > 2)
        viewer.data(mesh3_id).set_mesh(V3_out, F3);

    updateViz();

    /*for (int i=0; i<s1.corres.size(); i++){
        Eigen::RowVector3d color = Eigen::RowVector3d::Random(); 
        viewer.data().add_points(V1_out.row(s1.corres[i].first), color);
        viewer.data().add_points(V2_out.row(s1.corres[i].second), color);
    }

    for (int i=0; i<s2.corres.size(); i++){
        Eigen::RowVector3d color = Eigen::RowVector3d::Random(); 
        viewer.data().add_points(V2_out.row(s2.corres[i].first), color);
        viewer.data().add_points(V3_out.row(s2.corres[i].second), color);
    }*/


    //viewer.data().line_width = 5;
    //viewer.data().point_size = 10;
    //viewer.core().orthographic = true;

    viewer.core().set_rotation_type(igl::opengl::ViewerCore::ROTATION_TYPE_TRACKBALL);
    viewer.core().background_color = Eigen::Vector4f(202.0/255.0, 190.0/255.0, 232.0/255.0, 1.0);
    viewer.launch();
}