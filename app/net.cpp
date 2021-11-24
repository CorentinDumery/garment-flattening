

#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/readOBJ.h>

#include "net_param.h"

int main(int argc, char *argv[]){

    Eigen::MatrixXd V_2d_d, V_3d_d;
    Eigen::MatrixXi F;
    igl::readOBJ("../data/dress_front_cut.obj", V_3d_d, F);
    igl::readOBJ("../data/flat_dress.obj", V_2d_d, F); // note: could also be encoded as UV in a single OBJ

    int n_vertices = V_2d_d.rows();
    Eigen::MatrixXf V_2d(n_vertices, 2), V_3d(n_vertices, 3);
    V_2d.col(0) = V_2d_d.col(0).cast<float>();
    V_2d.col(1) = V_2d_d.col(1).cast<float>();
    V_3d = V_3d_d.cast<float>();

    V_3d = V_3d.rowwise() - V_3d.colwise().minCoeff();

    int n_fibers = 50; // TODO

    V_3d *= static_cast<double>(n_fibers) / (V_3d.maxCoeff() - V_3d.minCoeff());

    std::cout << "Mesh range:" << std::endl;
    std::cout << "X: " << V_3d.col(0).minCoeff() << " -> " << V_3d.col(0).maxCoeff()  << std::endl;
    std::cout << "Y: " << V_3d.col(1).minCoeff() << " -> " << V_3d.col(1).maxCoeff()  << std::endl;
    std::cout << "Z: " << V_3d.col(2).minCoeff() << " -> " << V_3d.col(2).maxCoeff()  << std::endl;

    /*V_2d.resize(V_3d.rows(), 2);
    V_2d.col(0) = V_3d.col(0);
    V_2d.col(1) = V_3d.col(1);*/


    NetParam net_param(F, V_3d, V_2d);

    net_param.computeFibers();

    
    Eigen::MatrixXd edge_begs, edge_ends;
    net_param.vizBoundaryEdges(edge_begs, edge_ends);

    auto fib_viz = net_param.vizFibers();
    std::vector<Eigen::MatrixXd> fiber_begs_list = fib_viz[0];
    std::vector<Eigen::MatrixXd> fiber_ends_list = fib_viz[1];




    // --- VISUALIZATION ---

   
    igl::opengl::glfw::Viewer viewer;

    igl::opengl::glfw::imgui::ImGuiMenu menu;
    menu.callback_draw_viewer_window = []() {};
    viewer.plugins.push_back(&menu);

    viewer.data().add_edges(edge_begs, edge_ends, Eigen::RowVector3d(1.0, 0.5, 0.6));
    viewer.data().add_edges(fiber_begs_list[0], fiber_ends_list[0], Eigen::RowVector3d(0.5, 1.0, 0.6));
    viewer.data().add_edges(fiber_begs_list[1], fiber_ends_list[1], Eigen::RowVector3d(0.5, 0.5, 1.0));
    viewer.data().set_mesh(V_2d_d, F);

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
        if (ImGui::Begin("Net")) {
            
            make_checkbox("Show mesh", viewer.data().show_lines);
            ImGui::End();
        }
    
    };

    //updateViz();


    viewer.data().line_width = 5;

    viewer.data().show_lines = false;

    viewer.data().point_size = 10;
    viewer.core().orthographic = true;

    viewer.core().set_rotation_type(igl::opengl::ViewerCore::ROTATION_TYPE_NO_ROTATION);
    viewer.core().background_color = Eigen::Vector4f(202.0/255.0, 190.0/255.0, 232.0/255.0, 1.0);
    viewer.launch();
}