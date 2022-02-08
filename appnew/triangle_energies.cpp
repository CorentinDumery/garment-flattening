/**
 * @author Corentin Dumery
 * @brief Given a reference and a deformed triangle, computes stretch/shear
 * using different methods. 
 * @date 2022-02-08
 * 
 */

#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/readOBJ.h>

#include "mccartney.h"
#include "param_utils.h"
#include "surface_net.h"

int main(int argc, char *argv[]){

    Eigen::MatrixXd V_2di(3,3), V_3di(3,3);
    Eigen::MatrixXi F(1,3);

    V_2di << -1.0, 1.0, 0.0,
            2.0,-1.0, 0.0,
            3.0, 2.0, 0.0;

    double scale = 2.0;
    V_2di *= scale;
    V_3di = V_2di;

    F << 0, 1, 2;    

    // --- VISUALIZATION ---

    double Esu, Esv, Er; // measured energies

    float Su_exp = 1.0;
    float Sv_exp = 1.0;
    float phiv_exp = 0.1;
    float shear_u_exp = 0.0;
    float shear_v_exp = 0.0;

    int align_mode = 1;
    int measure_mode = 1;

    igl::opengl::glfw::Viewer viewer;
    igl::opengl::glfw::imgui::ImGuiMenu menu;
    menu.callback_draw_viewer_window = []() {};
    viewer.plugins.push_back(&menu);

    viewer.append_mesh();
    viewer.append_mesh();
    int id_2d = viewer.data_list[0].id;
    int id_3d = viewer.data_list[1].id;

    auto updateViz = [&](){
        viewer.data(id_2d).clear_points();
        viewer.data(id_2d).clear_edges();
        viewer.data(id_3d).clear_edges();

        Eigen::MatrixXd V_2d(3,3);
        Eigen::MatrixXd V_3d(3,3);
        V_2d = V_2di;
        V_3d = V_3di;

        // -- Deform 3D triangle following input parameters --
        Eigen::Matrix3d stretch, shear, shear_bis;
        stretch << Su_exp, 0, 0,
                   0, Sv_exp, 0,
                   0,      0, 1;

        shear << 1, std::sin(phiv_exp), 0,
                 0, std::cos(phiv_exp), 0,
                   0,      0, 1;

        shear_bis << 1, shear_v_exp, 0,
                 shear_u_exp, 1, 0,
                   0,      0, 1;

        Eigen::Matrix3d Mcomp = shear * shear_bis * stretch;
        Eigen::Matrix3d V_2dt, V_3dt;
        V_2dt = V_2d;
        V_2dt.transposeInPlace();
        V_3dt = Mcomp * V_2dt;
        V_3dt.transposeInPlace();
        V_3d = V_3dt;

        if (align_mode == 0){ // procrustes
            Eigen::MatrixXd R_est;
            Eigen::VectorXd T_est;
            procrustes(V_2d, V_3d, R_est, T_est);

            V_3dt = Mcomp * V_2dt;
            Eigen::MatrixXd line3 = V_3dt.colwise() - T_est;
            V_3dt = (R_est.transpose() * line3);

            V_3dt.transposeInPlace();
            V_3d = V_3dt;
        }
        else { // align warp axis
            
            V_2d.row(1) -= V_2d.row(0);
            V_2d.row(2) -= V_2d.row(0);
            V_2d.row(0) -= V_2d.row(0);
            V_3d = move3Dto2D(V_3d);

            if (V_2d.row(0).maxCoeff() > 0 || V_2d.col(2).maxCoeff() > 0){
                std::cout << "ERROR: violated 2d assumptions" << std::endl;
            } 

            Eigen::RowVector3d B = V_2d.row(1);
            Eigen::RowVector3d C = V_2d.row(2);
            Eigen::RowVector3d BC = C - B;
            double d = BC(1) / BC(0); 
            double alpha = - B(1) / (d );
            alpha = - B(1) / BC(1);
            Eigen::RowVector3d X = B + alpha * BC;
            Eigen::RowVector3d Bp = V_3d.row(1);
            Eigen::RowVector3d Cp = V_3d.row(2);
            Eigen::RowVector3d Xp = Bp + alpha * (Cp - Bp);
            Eigen::Matrix3d R = computeRotation(Xp, X);

            V_3d = (R * V_3d.transpose()).transpose();
        }
        
        
        
        if (measure_mode == 0) computeFrameErrors(V_2d, V_3d, Esu, Esv, Er);
        else computeMcCartneyErrors(V_2d, V_3d, Esu, Esv, Er);


        // -- display matching nets on triangles (only for visualization) --
        SurfaceNet net_param(F, V_3d.cast<float>(), V_2d.cast<float>());
        net_param.computeNet();
        
        std::vector<std::vector<Eigen::MatrixXd>> fibs = net_param.vizNet();
        std::vector<Eigen::MatrixXd> fiber_begs_list = fibs[0];
        std::vector<Eigen::MatrixXd> fiber_ends_list = fibs[1];

        viewer.data(id_2d).add_edges(fiber_begs_list[0], fiber_ends_list[0], Eigen::RowVector3d(1.0, 0.0, 0.0));
        viewer.data(id_2d).add_edges(fiber_begs_list[1], fiber_ends_list[1], Eigen::RowVector3d(1.0, 0.0, 0.0));

        // Transport net to other triangle
        viewer.data(id_3d).add_edges(transportMatrix(V_2d, F, fiber_begs_list[0], V_3d), 
                                     transportMatrix(V_2d, F, fiber_ends_list[0], V_3d), 
                                     Eigen::RowVector3d(0.0, 0.0, 1.0));
        viewer.data(id_3d).add_edges(transportMatrix(V_2d, F, fiber_begs_list[1], V_3d), 
                                     transportMatrix(V_2d, F, fiber_ends_list[1], V_3d), 
                                     Eigen::RowVector3d(0.0, 0.0, 1.0));

        viewer.data(id_2d).set_mesh(V_2d, F);
        viewer.data(id_3d).set_mesh(V_3d, F);
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
        if (ImGui::Begin("Triangle energies")) {

            
            ImGui::Text("In this app, we consider the transformation from the red to");
            ImGui::Text("the blue triangle. We aim to measure stretch and shear");
            ImGui::Text("independently. We compare (1) frame transportation from red");
            ImGui::Text("to blue and (2) energies proposed by McCartney et al.");

            const char* dir_items2[2];
            dir_items2[0] = "Frame transportation";
            dir_items2[1] = "McCartney et al.";
            
            if (ImGui::ListBox("Measurement", &measure_mode, dir_items2, IM_ARRAYSIZE(dir_items2), 2)){
                updateViz();
            }

            ImGui::Text("The blue triangle is at an arbitrary location in 3D, but");
            ImGui::Text("we use either (1) the optimal rigid transformation or");
            ImGui::Text("(2) warp alignment to bring them together.");

            const char* dir_items[2];
            dir_items[0] = "Procrustes";
            dir_items[1] = "Align warp";
            
            if (ImGui::ListBox("Alignment", &align_mode, dir_items, IM_ARRAYSIZE(dir_items), 2)){
                updateViz();
            }

            ImGui::Separator();

            ImGui::Text("Transform blue triangle:");

            if (ImGui::SliderFloat("Stretch U", &Su_exp, 0.01, 2)){
                updateViz();
            }

            if (ImGui::SliderFloat("Stretch V", &Sv_exp, 0.01, 2)){
                updateViz();
            }

            if (ImGui::SliderFloat("Phi_v (Shear)", &phiv_exp, -3.14/2.0, 3.14/2.0)){
                updateViz();
            }

            if (ImGui::SliderFloat("Combined shear U", &shear_u_exp, 0.01, 2)){
                updateViz();
            }

            if (ImGui::SliderFloat("Combined shear V", &shear_v_exp, 0.01, 2)){
                updateViz();
            }

            ImGui::Separator();

            ImGui::Text("Measurements:");

            float w = ImGui::GetContentRegionAvailWidth();
            float p = ImGui::GetStyle().FramePadding.x;
            ImGui::SetCursorPosX((w - 2) * 0.06f);
            ImGui::Text("Stretch U",  ImVec2((w - p) / 3.f, 0));
            ImGui::SameLine(0, p);
            ImGui::SetCursorPosX((w - 2) * 0.36f);
            ImGui::Text("Stretch V",  ImVec2((w - p) / 3.f, 0));
            ImGui::SameLine(0, p);
            ImGui::SetCursorPosX((w - 2) * 0.72f);
            ImGui::Text("Shear",  ImVec2((w - p) / 3.f, 0));

            float draw_height = 200.0;
            ImVec2 size(ImGui::GetContentRegionAvailWidth(), draw_height);
            ImGui::InvisibleButton("canvas_energies", size);
            ImVec2 p0 = ImGui::GetItemRectMin();
            ImVec2 p1 = ImGui::GetItemRectMax();
            ImDrawList* draw_list = ImGui::GetWindowDrawList();
            draw_list->PushClipRect(p0, p1);

            std::vector<double> values = {Esu, Esv, Er};

            

            float max_val = 5.0;
            for (int i=0; i<3; i++){
                int err_viz = std::max(255 - static_cast<int>(255.0 * values[i] / max_val), 0);
                draw_list->AddRectFilled(
                        ImVec2(p0.x + ((float) (2*i)/6.0) * size.x, p0.y),
                        ImVec2(p0.x + (((float) (2*i + 1))/6.0) * size.x, p0.y + (values[i] / max_val) * size.y),
                        IM_COL32(255, err_viz, err_viz, 255));
            }

            draw_list->PopClipRect();

            ImGui::End();
        }
    };
    updateViz();

    viewer.data(id_2d).line_width = 5;
    viewer.data(id_2d).line_color = Eigen::RowVector4f(1.0, 0.0, 0.0, 1.0);
    viewer.data(id_3d).line_width = 5;
    viewer.data(id_3d).line_color = Eigen::RowVector4f(0.0, 0.0, 1.0, 1.0);
    viewer.core().orthographic = false;
    viewer.core().set_rotation_type(igl::opengl::ViewerCore::ROTATION_TYPE_NO_ROTATION);
    viewer.core().background_color = Eigen::Vector4f(202.0/255.0, 190.0/255.0, 232.0/255.0, 1.0);
    viewer.launch();
}