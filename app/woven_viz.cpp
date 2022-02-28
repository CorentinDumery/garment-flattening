/**
 * @author Corentin Dumery
 * @brief Anisotropic parameterization of a mesh with a disk topology. Visualizes
 * different stretch/shear coefficients.
 * @date 2022-02-16
 * 
 */
#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/readOBJ.h>
#include <igl/writeOBJ.h>
#include <igl/png/readPNG.h>
#include <igl/edges.h>
#include <igl/unproject_onto_mesh.h>
#include <igl/boundary_loop.h>
#include <igl/avg_edge_length.h>
#include <igl/doublearea.h>
#include <igl/jet.h>

#include "metrics.h"
#include <param/bary_optimizer.h>
#include <param/self_intersect.h>
#include "param/cloth_param.h"
#include "draw_colormap.h"

int main(int argc, char *argv[]){
    Eigen::MatrixXd V_3d, V_2d, V_2di;
    Eigen::MatrixXi F;

    if (argc > 1){
        igl::readOBJ(std::string(argv[1]), V_3d, F);
    }
    else {
        igl::readOBJ("../data/jumpsuit_multipose/front1.obj", V_3d, F);
    }

    double scale_f = 1.0;
    float scale_uv = 1.0;
    V_3d *= scale_f;

    ClothParam cp(V_3d, F, 0.00);
    V_2di = cp.getV2d();
    cp.paramAttempt(5);
    V_2d = cp.getV2d();

    V_3d.col(0) = V_3d.col(0).array() - (V_3d(0,0) - V_2d(0,0));
    V_3d.col(1) = V_3d.col(1).array() - (V_3d(0,1) - V_2d(0,1));

    std::cout << "F.rows(): " << F.rows() << std::endl;
    Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> R,G,B,A;
    igl::png::readPNG("../data/grid.png",R,G,B,A);


    // --- VISUALIZATION ---
    igl::opengl::glfw::Viewer viewer;
    int n_iterations = 5;
    int display_mode = 1;
    float stretch_f = 5.0;
    float edges_f = 1.0;
    bool show_debug_info = false;
    double anim_time = 0;
    double time_increment = 0.02;
    bool morphing = true;
    bool is_displaying_v2d = false;
    
    viewer.data().set_mesh(V_3d, F);

    viewer.data().set_uv(V_2d * scale_uv);
    viewer.data().show_texture = true;
    viewer.data().set_texture(R,G,B);

    igl::opengl::glfw::imgui::ImGuiMenu menu;
    menu.callback_draw_viewer_window = []() {};
    viewer.plugins.push_back(&menu);

    auto updateViz = [&](){
        viewer.data().clear_edges();
        viewer.data().clear_points();

        if (display_mode == 2){
            Eigen::VectorXd stretch_u_vec, stretch_v_vec;
            Eigen::MatrixXd stretch_u_colors, stretch_v_colors;
            measureStretchScore(V_2d, V_3d, F, stretch_u_vec, stretch_v_vec);
            igl::jet(stretch_u_vec, -0.5, 0.5, stretch_u_colors);
            viewer.data().set_colors(stretch_u_colors);
        }

        if (display_mode == 3){
            Eigen::VectorXd align_error;
            Eigen::MatrixXd align_colors;
            measureAlignmentScore(V_2d, V_3d, F, align_error);
            igl::jet(align_error, -0.5, 0.5, align_colors);
            viewer.data().set_colors(align_colors);
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
        if (ImGui::Begin("Woven Parameterization")) {
            float w = ImGui::GetContentRegionAvailWidth();
            float p = ImGui::GetStyle().FramePadding.x;

            if (ImGui::Button("Morphing",  ImVec2((w - p) / 3.f, 0))){
                morphing = true;
                is_displaying_v2d = false;
            }
            ImGui::SameLine(0, p);
            if (ImGui::Button("3D",  ImVec2((w - p) / 3.f, 0))){
                morphing = false;
                viewer.data().set_mesh(V_3d, F);
                viewer.core().align_camera_center(V_3d,F);
                viewer.core().lighting_factor = 1.0;
                is_displaying_v2d = false;
            }
            
            ImGui::SameLine(0, p);
            if (ImGui::Button("2D",  ImVec2((w - p) / 3.f, 0))){
                morphing = false;
                viewer.data().set_mesh(V_2d, F);
                viewer.core().align_camera_center(V_2d,F);
                viewer.core().lighting_factor = 0.0;
                is_displaying_v2d = true;
            }

            /*make_checkbox("Show mesh", viewer.data().show_lines);
            make_checkbox("Show texture", viewer.data().show_texture);*/

            const char* dir_items[4];
            dir_items[0] = "Mesh";
            dir_items[1] = "Woven texture";
            dir_items[2] = "Stretch";
            dir_items[3] = "Alignment";
            
            if (ImGui::ListBox("Display", &display_mode, dir_items, IM_ARRAYSIZE(dir_items), 4)){
                viewer.data().show_lines = display_mode == 0;
                viewer.data().show_texture = display_mode == 1;
                if (display_mode < 2) viewer.data().set_colors(Eigen::RowVector3d(1.0, 1.0, 1.0));
                updateViz();
            }


            ImGui::Text(std::to_string(-0.5).substr(0,5).c_str());
            std::string text = std::to_string(0.5).substr(0,5);
            auto textWidth   = ImGui::CalcTextSize(text.c_str()).x;
            ImGui::SameLine();
            ImGui::SetCursorPosX((ImGui::GetWindowSize().x - 1) * 0.5f);
            ImGui::Text("0");
            ImGui::SameLine();
            ImGui::SetCursorPosX(ImGui::GetWindowSize().x - 1.5*textWidth);
            ImGui::Text(text.c_str());
            drawColormap(igl::COLOR_MAP_TYPE_JET);
                                
            ImGui::SliderFloat("Light factor", &viewer.core().lighting_factor, 0.0f, 5.0f, "%.3f");
            
            if (ImGui::SliderFloat("Scale fibers", &scale_uv, 0.1f, 3.0f, "%.3f")){
                display_mode = 1;
                viewer.data().set_uv(V_2d * scale_uv);
            }

            ImGui::Separator();

            if (ImGui::SliderInt("# iterations", &n_iterations, 0, 40, "%.3f")){
                cp.setV2d(V_2di);
                cp.paramAttempt(n_iterations);
                V_2d = cp.getV2d();
                if (is_displaying_v2d) viewer.data().set_mesh(V_2d, F);
                viewer.data().set_uv(V_2d * scale_uv);
                updateViz();
            }

            if (ImGui::SliderFloat("Stretch penalty", &stretch_f, 0.0f, 50.0f, "%.3f")){
                cp.setCoeffs(stretch_f, edges_f);
                cp.setV2d(V_2di);
                cp.paramAttempt(n_iterations);
                V_2d = cp.getV2d();
                if (is_displaying_v2d) viewer.data().set_mesh(V_2d, F);
                viewer.data().set_uv(V_2d * scale_uv);
                updateViz();
            }

            if (ImGui::SliderFloat("Rigidity penalty", &edges_f, 0.0f, 10.0f, "%.3f")){
                cp.setCoeffs(stretch_f, edges_f);
                cp.setV2d(V_2di);
                cp.paramAttempt(n_iterations);
                V_2d = cp.getV2d();
                if (is_displaying_v2d) viewer.data().set_mesh(V_2d, F);
                viewer.data().set_uv(V_2d * scale_uv);
                updateViz();
            }

            ImGui::Separator();
            ImGui::Checkbox("Debug info", &show_debug_info);
            
            if (show_debug_info){
                ImGui::Text("Self intersects: %s", selfIntersect(V_2d, cp.getBnd()) ? "true" : "false");
                Eigen::VectorXd M;
                igl::doublearea(V_3d, F, M);
                ImGui::Text("3D Mesh area: %f", M.sum());
                ImGui::Text("3D Mesh range:");
                ImGui::Text("X: %f, (%f -> %f)", V_3d.col(0).maxCoeff(), V_3d.col(0).minCoeff(), V_3d.col(0).minCoeff(), V_3d.col(0).maxCoeff());
                ImGui::Text("Y: %f, (%f -> %f)", V_3d.col(1).maxCoeff(), V_3d.col(1).minCoeff(), V_3d.col(1).minCoeff(), V_3d.col(1).maxCoeff());
                ImGui::Text("Z: %f, (%f -> %f)", V_3d.col(2).maxCoeff(), V_3d.col(2).minCoeff(), V_3d.col(2).minCoeff(), V_3d.col(2).maxCoeff());
                
                ImGui::Separator();
                igl::doublearea(V_2d, F, M);
                ImGui::Text("2D Mesh area: %f", M.sum());
                ImGui::Text("2D Mesh range:");
                ImGui::Text("U: %f, (%f -> %f)", V_2d.col(0).maxCoeff(), V_2d.col(0).minCoeff(), V_2d.col(0).minCoeff(), V_2d.col(0).maxCoeff());
                ImGui::Text("V: %f, (%f -> %f)", V_2d.col(1).maxCoeff(), V_2d.col(1).minCoeff(), V_2d.col(1).minCoeff(), V_2d.col(1).maxCoeff());
            }
            ImGui::Separator();

            if (ImGui::Button("Save UV", ImVec2(-1, 0))){
                igl::writeOBJ("../saved_uv.obj", V_2d, F);
            }
            ImGui::End();
        }
    
    };

    viewer.callback_post_draw = [&](igl::opengl::glfw::Viewer &v) {
        if (morphing){
            anim_time += time_increment;
            double coeff = std::pow(std::sin(anim_time), 2);
            viewer.core().lighting_factor = coeff;
            Eigen::MatrixXd V = coeff * V_3d + (1-coeff) * V_2d;
            viewer.data().set_mesh(V, F);
            updateViz();
        }
        return true;
    };

    std::vector<int> selected_vs = {-1, -1};
    int next_select = 0;

    /*viewer.callback_key_down = [&V_3d, &F, &selected_vs, &next_select, &bo](igl::opengl::glfw::Viewer& viewer, unsigned char key, int modifiers)->bool{
        double mouse_x = viewer.current_mouse_x;
        double mouse_y = viewer.core().viewport(3) - viewer.current_mouse_y;
        int fid = -2;
        Eigen::RowVector3d bc;
        if(igl::unproject_onto_mesh(Eigen::Vector2f(mouse_x, mouse_y), viewer.core().view,
                                    viewer.core().proj, viewer.core().viewport, V_3d, F, fid, bc)){
            int argmax = -1;
            if (bc(0) >= bc(1) && bc(0) >= bc(2)) argmax = 0;
            if (bc(1) >= bc(0) && bc(1) >= bc(2)) argmax = 1;
            if (bc(2) >= bc(0) && bc(2) >= bc(1)) argmax = 2;
            std::cout << F(fid, argmax) << std::endl;

            selected_vs[next_select] = F(fid, argmax);
            next_select = 1 - next_select;

            bo.setSelectedVertices(selected_vs);

            viewer.data().clear_points();
            for (int i=0; i<selected_vs.size(); i++){
                if (selected_vs[i] == -1) continue;
                viewer.data().add_points(V_3d.row(selected_vs[i]), Eigen::RowVector3d(1.0, 0, 0));
            }
            return true;
        }
        return false;
    };*/

    updateViz();

    viewer.data().set_colors(Eigen::RowVector3d(1.0, 1.0, 1.0));
    viewer.data().line_width = 3;
    viewer.data().show_lines = false;
    viewer.core().is_animating = true;
    viewer.core().set_rotation_type(igl::opengl::ViewerCore::ROTATION_TYPE_TRACKBALL);
    viewer.core().background_color = Eigen::Vector4f(202.0/255.0, 190.0/255.0, 232.0/255.0, 1.0);
    viewer.launch();
}