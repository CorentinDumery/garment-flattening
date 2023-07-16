/**
 * @author Corentin Dumery 
 * @brief Split a mesh in two. 
 * @date 2022-07-04
 * 
 */

#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/readOBJ.h>
#include <igl/avg_edge_length.h>

#include "cut_on_plane.h"


int main(int argc, char *argv[]){

    std::string path_3d;
    path_3d = "../data/dress.ply";

    if (argc >= 2){path_3d = argv[1];}

    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    if (path_3d.substr(path_3d.length() - 3) == "obj")
        igl::readOBJ(path_3d, V, F);
    else if (path_3d.substr(path_3d.length() - 3) == "ply")
        igl::readPLY(path_3d, V, F);
    else {
        std::cout << "Format not recognized" << std::endl;
        return -1;
    }

    Eigen::MatrixXd V_plane(4, 3);
    V_plane << -1, -1, 0,
                -1, 1, 0,
                1, 1, 0,
                1, -1, 0;

    V_plane *= 1000;

    Eigen::MatrixXi F_plane(2, 3);
    F_plane << 0, 1, 2,
               0, 2, 3;

    std::vector<Eigen::MatrixXd> V_list;
    std::vector<Eigen::MatrixXi> F_list;
    Eigen::VectorXi cut0, cut1; 
    cutMeshOnPlane(V, F, V_list, F_list, cut0, cut1);

    Eigen::MatrixXd V0, V1;
    Eigen::MatrixXi F0, F1;
    V0 = V_list[0];
    F0 = F_list[0];
    V1 = V_list[1];
    F1 = F_list[1];

    igl::writeOBJ("comp0.obj", V0, F0);
    igl::writeOBJ("comp1.obj", V1, F1);

    double scale = (V0.col(2).maxCoeff() - V0.col(2).minCoeff());

    // --- VISUALIZATION ---

    igl::opengl::glfw::Viewer viewer;
    viewer.append_mesh();
    viewer.append_mesh();
    int mesh1_id = viewer.data_list[1].id;
    int mesh2_id = viewer.data_list[2].id;
    int plane_id = viewer.data_list[0].id;

    // viz input mesh
    //viewer.data(mesh1_id).set_mesh(V, F);

    // --- Viz cut ---
    
    Eigen::MatrixXd colors = Eigen::MatrixXd::Random(cut0.rows(), 3);
    Eigen::MatrixXd points0(cut0.rows(), 3), points1(cut0.rows(), 3);
    
    for (int i=0; i<cut0.rows(); i++){
        points0.row(i) = V_list[0].row(cut0(i));
        points1.row(i) = V_list[1].row(cut1(i));
    }
    
    Eigen::MatrixXd points1_offset = points1.rowwise() + Eigen::RowVector3d(0, 0, scale/5);
    viewer.data(mesh1_id).add_points(points1_offset, colors);
    viewer.data(mesh1_id).add_points(points0, colors);
    viewer.data(mesh1_id).point_size = 12;

    // --- ---

    viewer.data(mesh1_id).set_mesh(V0, F0);
    Eigen::MatrixXd V1_offset = V1.rowwise() + Eigen::RowVector3d(0, 0, scale/5);
    viewer.data(mesh2_id).set_mesh(V1_offset, F1);

    viewer.data(plane_id).set_mesh(V_plane, F_plane);

    igl::opengl::glfw::imgui::ImGuiMenu menu;
    menu.callback_draw_viewer_window = []() {};
    viewer.plugins.push_back(&menu);

    auto updateViz = [&](){
        
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
        if (ImGui::Begin("Multiple poses")) {
            make_checkbox("Show plane", viewer.data(plane_id).is_visible);
            ImGui::End();
        }
    
    };

    updateViz();

    viewer.data(plane_id).is_visible = false;
    viewer.core().background_color = Eigen::Vector4f(1.0, 1.0, 1.0, 1.0);
    viewer.launch();
}
