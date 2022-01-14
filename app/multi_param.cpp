
#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/readOBJ.h>

#include "param/multi_patch_param.h"

int main(int argc, char *argv[]){

    Eigen::MatrixXd V1, V2, V3, V1_out, V2_out, V3_out;
    Eigen::MatrixXi F1, F2, F3;
    std::vector<Eigen::MatrixXd> out_vec;

    igl::readOBJ("../data/patches/patch_3D_1.obj", V1, F1);
    igl::readOBJ("../data/patches/patch_3D_3.obj", V2, F2);
    igl::readOBJ("../data/patches/patch_3D_4.obj", V3, F3);

    finalParamMultiPatch({V1, V2, V3}, {F1, F2, F3}, 
                          {{}, {}, {}}, // dart dupl
                          {{}, {}, {}}, // dart tips
                          {}, // seams
                          out_vec);
    
    V1_out = out_vec[0];
    V2_out = out_vec[1];
    V3_out = out_vec[2];

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
    viewer.data(mesh3_id).set_mesh(V3_out, F3);

    updateViz();


    //viewer.data().line_width = 5;
    //viewer.data().point_size = 10;
    //viewer.core().orthographic = true;

    viewer.core().set_rotation_type(igl::opengl::ViewerCore::ROTATION_TYPE_TRACKBALL);
    viewer.core().background_color = Eigen::Vector4f(202.0/255.0, 190.0/255.0, 232.0/255.0, 1.0);
    viewer.launch();
}