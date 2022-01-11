

#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/readOBJ.h>
#include <igl/boundary_loop.h>

#include "param/param_utils.h"
#include "param/auto_select.h"

int main(int argc, char *argv[]){

    Eigen::MatrixXd V_2d, V_3d;
    Eigen::MatrixXi F;

    std::string input = "../data/buggy_patch.obj";
    if (argc > 1) input = std::string(argv[1]);
    igl::readOBJ(input, V_3d, F);

    
    Eigen::VectorXi bnd;
    igl::boundary_loop(F, bnd);

    std::vector<int> selec = autoSelect(V_3d, bnd);

    V_2d = paramLSCM(V_3d, F);

    Eigen::RowVector3d from = V_2d.row(selec[0]) - V_2d.row(selec[1]);
    Eigen::RowVector3d to(1.0, 0, 0);  
    Eigen::Matrix3d R = computeRotation(from, to);
    V_2d = (R * V_2d.transpose()).transpose();


    // --- VISUALIZATION ---

    igl::opengl::glfw::Viewer viewer;

    igl::opengl::glfw::imgui::ImGuiMenu menu;
    menu.callback_draw_viewer_window = []() {};
    viewer.plugins.push_back(&menu);

    int curr_d = 3;

    auto updateViz = [&](){
        viewer.data().clear_edges();
        viewer.data().clear_points();

        for (int i=0; i<bnd.rows(); i++){
            //viewer.data().add_points(V.row(bnd(i)), Eigen::RowVector3d(0,0,0));
        }

        for (int i=0; i<selec.size(); i++){
            //viewer.data().add_points(V.row(bnd(i)), Eigen::RowVector3d(0,0,0));
            if (curr_d ==3)
                viewer.data().add_points(V_3d.row(selec[i]), Eigen::RowVector3d(1,0,0));
            if (curr_d ==2)
                viewer.data().add_points(V_2d.row(selec[i]), Eigen::RowVector3d(1,0,0));
        }
        
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
        if (ImGui::Begin("Net")) {

            float w = ImGui::GetContentRegionAvailWidth();
            float p = ImGui::GetStyle().FramePadding.x;
            
            if (ImGui::Button("3D",  ImVec2((w - p) / 2.f, 0))){
                viewer.data().set_mesh(V_3d, F);
                curr_d = 3;            
                updateViz();
            }
            
            ImGui::SameLine(0, p);
            if (ImGui::Button("2D",  ImVec2((w - p) / 2.f, 0))){
                viewer.data().set_mesh(V_2d, F);
                curr_d = 2;
                updateViz();
            }
            
            ImGui::End();
        }
    
    };

    viewer.data().set_mesh(V_3d, F);
    updateViz();


    //viewer.data().line_width = 5;
    //viewer.data().point_size = 10;
    //viewer.core().orthographic = true;

    viewer.core().set_rotation_type(igl::opengl::ViewerCore::ROTATION_TYPE_TRACKBALL);
    viewer.core().background_color = Eigen::Vector4f(202.0/255.0, 190.0/255.0, 232.0/255.0, 1.0);
    viewer.launch();
}