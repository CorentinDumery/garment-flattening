

#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/readOBJ.h>
#include <igl/opengl/glfw/background_window.h>
#include <igl/unproject_onto_mesh.h>

#include "draw_colormap.h" 
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

    int n_fibers = 20; // TODO
    V_2d *= static_cast<float>(n_fibers) / (V_3d.maxCoeff() - V_3d.minCoeff());
    V_2d_d *= static_cast<float>(n_fibers) / (V_3d.maxCoeff() - V_3d.minCoeff());
    V_3d *= static_cast<float>(n_fibers) / (V_3d.maxCoeff() - V_3d.minCoeff());

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
    std::vector<std::vector<int>> selected_fibers;
    Eigen::RowVector2d cursor_pos(-1,-1);

    std::vector<Eigen::MatrixXd> colors_per_axis;
   
    igl::opengl::glfw::Viewer viewer;

    viewer.append_mesh();
    viewer.append_mesh();

    int mesh_id = viewer.data_list[2].id;
    int line_id = viewer.data_list[1].id;
    int highlighted_id = viewer.data_list[0].id;

    igl::opengl::glfw::imgui::ImGuiMenu menu;
    menu.callback_draw_viewer_window = []() {};
    viewer.plugins.push_back(&menu);

    

    auto updateViz = [&](){

        // mesh_id
        viewer.data(mesh_id).clear_edges();
        viewer.data(mesh_id).clear_points();
        viewer.data(mesh_id).add_edges(edge_begs, edge_ends, Eigen::RowVector3d(1.0, 0.5, 0.6));

        if (colors_per_axis.size() == 0){
            viewer.data(mesh_id).add_edges(fiber_begs_list[0], fiber_ends_list[0], Eigen::RowVector3d(0.5, 1.0, 0.6));
            viewer.data(mesh_id).add_edges(fiber_begs_list[1], fiber_ends_list[1], Eigen::RowVector3d(0.5, 0.5, 1.0));
        }
        else {
            viewer.data(mesh_id).add_edges(fiber_begs_list[0], fiber_ends_list[0], colors_per_axis[0]);
            viewer.data(mesh_id).add_edges(fiber_begs_list[1], fiber_ends_list[1], colors_per_axis[1]);
        }
        viewer.data(mesh_id).set_mesh(V_2d_d, F);
        viewer.data(mesh_id).set_colors(Eigen::RowVector3d(1.0,1.0,1.0));

        
        // line_id
        viewer.data(line_id).clear_edges();
        Eigen::RowVector3d start(start_u, start_v, 1.0);
        Eigen::RowVector3d end = start;
        end(measure_axis) += measure_length;
        viewer.data(line_id).add_edges(start, end, Eigen::RowVector3d(1.0, 0.0, 0.0));


        // highlighted id
        viewer.data(highlighted_id).clear_edges();
        for (int i=0; i<selected_fibers.size(); i++){
            for (int j=0; j<selected_fibers[i].size(); j++){
                if (selected_fibers[i][j] >= 0){
                    viewer.data(highlighted_id).add_edges(fiber_begs_list[i].row(selected_fibers[i][j]), 
                                                          fiber_ends_list[i].row(selected_fibers[i][j]),
                                                          Eigen::RowVector3d(0.0, 0.0, 0.0));
                }
            }    
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
                start_point = Eigen::RowVector2f(start_u, start_v);
                start_mat = start_point;
                net_param.fromInitToRenderCoords(start_mat);

                end_point = start_point;
                end_point(measure_axis) += measure_length;
                end_mat = end_point;
                net_param.fromInitToRenderCoords(end_mat);
                Eigen::RowVector2f temp_start, temp_end;
                temp_start.row(0) = start_mat.row(0).cast<float>();
                temp_end.row(0) = end_mat.row(0).cast<float>();
                measurement = (float) net_param.simpleMeasureFiber(temp_start, temp_end);
                updateViz();
            }

            make_checkbox("Show measurement line", viewer.data(line_id).is_visible);
            make_checkbox("Show highlighted", viewer.data(highlighted_id).is_visible);

            ImGui::Text("Measurement: %f", measurement);
            first_menu_call = false;

            if (ImGui::Button("Update measurement")){
                //net_param.alternativeRenderingAttempt(start_point(0), start_point(1), measure_length);
                net_param.otherRenderingAttempt(start_point(0), start_point(1), measure_length);
            }

            if (ImGui::Button("Simple measurement")){
                //net_param.alternativeRenderingAttempt(start_point(0), start_point(1), measure_length);

                start_point = Eigen::RowVector2f(start_u, start_v);
                start_mat = start_point;
                net_param.fromInitToRenderCoords(start_mat);

                end_point = start_point;
                end_point(measure_axis) += measure_length;
                end_mat = end_point;
                net_param.fromInitToRenderCoords(end_mat);
                Eigen::RowVector2f temp_start, temp_end;
                temp_start.row(0) = start_mat.row(0).cast<float>();
                temp_end.row(0) = end_mat.row(0).cast<float>();
                net_param.simpleMeasureFiber(temp_start, temp_end);
            }

            if (ImGui::Button("Color fibers")){
                net_param.computeFibers();
                colors_per_axis = net_param.computeFibersColor();
                updateViz();
                
            }

            ImGui::Text(std::to_string(0.75).substr(0,5).c_str());
            std::string text = std::to_string(1.25).substr(0,5);
            auto textWidth   = ImGui::CalcTextSize(text.c_str()).x;
            ImGui::SameLine();
            ImGui::SetCursorPosX((ImGui::GetWindowSize().x - 1) * 0.5f);
            ImGui::Text("1");
            ImGui::SameLine();
            ImGui::SetCursorPosX(ImGui::GetWindowSize().x - 1.5*textWidth);
            ImGui::Text(text.c_str());
            drawColormap(igl::COLOR_MAP_TYPE_JET);
    
            if (ImGui::Button("Adjust UV")){
                net_param.adjustVertices();
                V_2d_d = net_param.getV2d().cast<double>();
                //net_param.f romRenderToInitCoords(V_2d_d);
                updateViz();
            }

            ImGui::Text("Cursor position on mesh: %f, %f", cursor_pos(0), cursor_pos(1));

            ImGui::End();
        }

        
    };

    viewer.callback_mouse_move = [&V_2d_d, &F, &mesh_id, &net_param, &selected_fibers, &updateViz, &cursor_pos](igl::opengl::glfw::Viewer& viewer, int mouse_x, int mouse_y){
        int fid = -2;
        Eigen::RowVector3d bc;
        if(igl::unproject_onto_mesh(Eigen::Vector2f(mouse_x, viewer.core().viewport(3) - mouse_y), viewer.core().view,
            viewer.core().proj, viewer.core().viewport, V_2d_d, F, fid, bc)){

            Eigen::RowVector2d p = V_2d_d.row(F(fid, 0)) * bc(0) + V_2d_d.row(F(fid, 1)) * bc(1) + V_2d_d.row(F(fid, 2)) * bc(2);

            cursor_pos = p;

            selected_fibers = net_param.nearestFibers(p.cast<float>()); 

            viewer.data(mesh_id).clear_points(); 
            viewer.data(mesh_id).add_points(p, Eigen::RowVector3d(0,0,0));
            updateViz();
            return true;
        }
        return false;
    };

    

    updateViz();

    viewer.data(mesh_id).line_width = 5;
    viewer.data(line_id).line_width = 15;
    viewer.data(highlighted_id).line_width = 20;

    viewer.data(mesh_id).show_lines = false;

    viewer.data(line_id).is_visible = false;

    viewer.data(mesh_id).point_size = 10;
    viewer.core().orthographic = true;

    viewer.core().set_rotation_type(igl::opengl::ViewerCore::ROTATION_TYPE_NO_ROTATION);
    viewer.core().background_color = Eigen::Vector4f(202.0/255.0, 190.0/255.0, 232.0/255.0, 1.0);
    //viewer.core().background_color = Eigen::Vector4f(0.0, 0.0, 0.0, 1.0);
    viewer.launch();
}