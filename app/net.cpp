

#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/readOBJ.h>
#include <igl/opengl/glfw/background_window.h>

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

    float start_u = 0, start_v=0, measure_length = 1.0;
    int measure_axis = 0;
    double measurement = -1;
    bool first_menu_call = true;
    Eigen::MatrixXf start_mat, end_mat;
    Eigen::RowVector2f start_point(start_u, start_v);
    Eigen::RowVector2f end_point;
   
    igl::opengl::glfw::Viewer viewer;

    viewer.append_mesh();

    int mesh_id = viewer.data_list[0].id;
    int line_id = viewer.data_list[1].id;

    igl::opengl::glfw::imgui::ImGuiMenu menu;
    menu.callback_draw_viewer_window = []() {};
    viewer.plugins.push_back(&menu);

    

    auto updateViz = [&](){
        // mesh_id
        viewer.data(mesh_id).clear_edges();
        viewer.data(mesh_id).clear_points();
        viewer.data(mesh_id).add_edges(edge_begs, edge_ends, Eigen::RowVector3d(1.0, 0.5, 0.6));
        viewer.data(mesh_id).add_edges(fiber_begs_list[0], fiber_ends_list[0], Eigen::RowVector3d(0.5, 1.0, 0.6));
        viewer.data(mesh_id).add_edges(fiber_begs_list[1], fiber_ends_list[1], Eigen::RowVector3d(0.5, 0.5, 1.0));
        viewer.data(mesh_id).set_mesh(V_2d_d, F);

        
        // line_id
        viewer.data(line_id).clear_edges();
        Eigen::RowVector3d start(start_u, start_v, 1.0);
        Eigen::RowVector3d end = start;
        end(measure_axis) += measure_length;
        viewer.data(line_id).add_edges(start, end, Eigen::RowVector3d(1.0, 0.0, 0.0));
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
            
            make_checkbox("Show mesh", viewer.data(mesh_id).show_lines);

            bool need_update = false;

            if (ImGui::SliderFloat("Starting U", &start_u, V_2d.col(0).minCoeff(), V_2d.col(0).maxCoeff())){
                need_update = true;
            }
            if (ImGui::SliderFloat("Starting V", &start_v, V_2d.col(1).minCoeff(), V_2d.col(1).maxCoeff())){
                need_update = true;
            }
            if (ImGui::SliderFloat("Length", &measure_length, 0, V_2d.maxCoeff()-V_2d.minCoeff())){
                need_update = true;
            }
            if (ImGui::Button("Change axis")){
                measure_axis = (measure_axis + 1) % 2;
                need_update = true;
            }

            if (need_update || first_menu_call){
                start_mat = start_point;
                net_param.fromInitToRenderCoords(start_mat);

                end_point = start_point;
                end_point(measure_axis) += measure_length;
                end_mat = end_point;
                net_param.fromInitToRenderCoords(end_mat);
                Eigen::RowVector2d temp_start, temp_end;
                temp_start.row(0) = start_mat.row(0).cast<double>();
                temp_end.row(0) = end_mat.row(0).cast<double>();
                
                updateViz();
            }

            ImGui::Text("Measurement: %f", measurement);
            first_menu_call = false;
            ImGui::End();
        }
    
    };


    updateViz();

    viewer.data(mesh_id).line_width = 5;
    viewer.data(line_id).line_width = 15;

    viewer.data(mesh_id).show_lines = false;

    viewer.data(mesh_id).point_size = 10;
    viewer.core().orthographic = true;

    viewer.core().set_rotation_type(igl::opengl::ViewerCore::ROTATION_TYPE_NO_ROTATION);
    viewer.core().background_color = Eigen::Vector4f(202.0/255.0, 190.0/255.0, 232.0/255.0, 1.0);
    //viewer.core().background_color = Eigen::Vector4f(0.0, 0.0, 0.0, 1.0);
    viewer.launch();
}