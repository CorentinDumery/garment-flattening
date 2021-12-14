

#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/readOBJ.h>

#include "net_param.h"
#include "mccartney.h"

int main(int argc, char *argv[]){

    double desired_fibers = 2.0;

    Eigen::MatrixXd V_2di(3,3), V_3di(3,3);
    Eigen::MatrixXi F(1,3);

    V_2di << -1.0, 1.0, 0.0,
            2.0,-1.0, 0.0,
            3.0, 2.0, 0.0;

    V_2di *= desired_fibers;

    V_3di = V_2di;

    //V_3di.row(1) = Eigen::RowVector3d(1.0, -2.0, 0.0);

    F << 0, 1, 2;    

    // --- VISUALIZATION ---

    double Esu=0.0, Esv=0.0, Er=0.0;

    float control_u=0, control_v=0;

    float Su_exp = 1.0;
    float Sv_exp = 1.0;
    float phiv_exp = 0.0;
    float shear_u_exp = 0.0;
    float shear_v_exp = 0.0;

    igl::opengl::glfw::Viewer viewer;

    igl::opengl::glfw::imgui::ImGuiMenu menu;
    menu.callback_draw_viewer_window = []() {};
    viewer.plugins.push_back(&menu);

    viewer.append_mesh();
    viewer.append_mesh();
    int id_2d = viewer.data_list[0].id;
    int id_3d_left = viewer.data_list[1].id;
    int id_3d_right = viewer.data_list[2].id;

    unsigned int left_view, right_view;
    viewer.callback_init = [&](igl::opengl::glfw::Viewer &){
        /*viewer.core().viewport = Eigen::Vector4f(0, 0, 640, 800);
        left_view = viewer.core_list[0].id;
        right_view = viewer.append_core(Eigen::Vector4f(640, 0, 640, 800));*/

        viewer.core().viewport = Eigen::Vector4f(0, 0, 640, 800);
        left_view = viewer.core_list[0].id;
        right_view = viewer.append_core(Eigen::Vector4f(640, 0, 640, 800));

        viewer.data(id_2d).set_visible(false, right_view);
        viewer.data(id_3d_left).set_visible(false, right_view);
        viewer.data(id_3d_right).set_visible(false, left_view);
        return false;
    };

    viewer.callback_post_resize = [&](igl::opengl::glfw::Viewer &v, int w, int h) {
        v.core( left_view).viewport = Eigen::Vector4f(0, 0, w / 2, h);
        v.core(right_view).viewport = Eigen::Vector4f(w / 2, 0, w - (w / 2), h);
        return true;
    };



    auto updateViz = [&](){
        viewer.data(id_2d).clear_points();
        viewer.data(id_2d).clear_edges();
        viewer.data(id_3d_left).clear_edges();
        viewer.data(id_3d_right).clear_edges();

        Eigen::MatrixXd V_2d(3,3);
        V_2d = V_2di;
        V_2d.rowwise() -= V_2di.row(0); 

        V_2d(1, 0) += control_u; // should be V_3d?
        V_2d(2, 0) += control_u;  
        V_2d(2, 1) += control_v;


        Eigen::Matrix3d stretch, shear, shear_bis;
        stretch << Su_exp, 0, 0,
                   0, Sv_exp, 0,
                   0,      0, 1;

        shear << 1, 0, 0,
                 std::sin(phiv_exp), std::cos(phiv_exp), 0,
                   0,      0, 1;

        shear_bis << 1, shear_u_exp, 0,
                 shear_v_exp, 1, 0,
                   0,      0, 1;

        Eigen::Matrix3d Mcomp = stretch * shear * shear_bis;
        V_2d = V_2d * Mcomp;




        Eigen::MatrixXd V_3d = computeMcCartneyErrors(V_2d, V_3di, Esu, Esv, Er);


        NetParam net_param(F, V_3d.cast<float>(), V_2d.cast<float>());
        net_param.computeFibers();
        
        std::vector<std::vector<Eigen::MatrixXd>> fibs = net_param.vizFibers();
        std::vector<Eigen::MatrixXd> fiber_begs_list = fibs[0];
        std::vector<Eigen::MatrixXd> fiber_ends_list = fibs[1];


        viewer.data(id_2d).add_edges(fiber_begs_list[0], fiber_ends_list[0], Eigen::RowVector3d(1.0, 0.0, 0.0));
        viewer.data(id_2d).add_edges(fiber_begs_list[1], fiber_ends_list[1], Eigen::RowVector3d(1.0, 0.0, 0.0));


        // Transport net to other triangle
        viewer.data(id_3d_left).add_edges(transportMatrix(V_2d, F, fiber_begs_list[0], V_3d), 
                                          transportMatrix(V_2d, F, fiber_ends_list[0], V_3d), 
                                          Eigen::RowVector3d(0.0, 0.0, 1.0));
        viewer.data(id_3d_left).add_edges(transportMatrix(V_2d, F, fiber_begs_list[1], V_3d), 
                                          transportMatrix(V_2d, F, fiber_ends_list[1], V_3d), 
                                          Eigen::RowVector3d(0.0, 0.0, 1.0));

        viewer.data(id_3d_right).add_edges(transportMatrix(V_2d, F, fiber_begs_list[0], V_3di), 
                                           transportMatrix(V_2d, F, fiber_ends_list[0], V_3di), 
                                           Eigen::RowVector3d(0.0, 0.0, 1.0));
        viewer.data(id_3d_right).add_edges(transportMatrix(V_2d, F, fiber_begs_list[1], V_3di), 
                                           transportMatrix(V_2d, F, fiber_ends_list[1], V_3di), 
                                           Eigen::RowVector3d(0.0, 0.0, 1.0));

        viewer.data(id_2d).set_mesh(V_2d, F);
        viewer.data(id_3d_left).set_mesh(V_3d, F);
        viewer.data(id_3d_right).set_mesh(V_3di, F);        
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

            if (ImGui::Button("Hide Left 3D")){
                viewer.data(id_3d_left).set_visible(false, left_view);
            }
            if (ImGui::Button("Show Left 3D")){
                viewer.data(id_3d_left).set_visible(true, left_view);
            }
            
            if (ImGui::Button("Randomize 3d tri")){
                V_3di = V_2di.array() + (Eigen::MatrixXd::Random(3,3).array() - 0.5) * 0.5 * desired_fibers;
                updateViz();
            }

            if (ImGui::SliderFloat("Stretch horizontally", &control_u, 0, 10*desired_fibers)){
                updateViz();
            }

            if (ImGui::SliderFloat("Stretch vertically", &control_v, 0, 10*desired_fibers)){
                updateViz();
            }

            if (ImGui::SliderFloat("Stretch U", &Su_exp, 0.01, 2)){
                updateViz();
            }

            if (ImGui::SliderFloat("Stretch V", &Sv_exp, 0.01, 2)){
                updateViz();
            }

            if (ImGui::SliderFloat("Phi_v (Shear)", &phiv_exp, -3.14, 3.14)){
                updateViz();
            }

            ImGui::Separator();

            if (ImGui::SliderFloat("Shear U", &shear_u_exp, -1.0, 1.0)){
                updateViz();
            }
            
            if (ImGui::SliderFloat("Shear V", &shear_v_exp, -1.0, 1.0)){
                updateViz();
            }


            ImGui::Separator();

            float draw_height = 200.0;
            ImVec2 size(ImGui::GetContentRegionAvailWidth(), draw_height);
            ImGui::InvisibleButton("canvas_energies", size);
            ImVec2 p0 = ImGui::GetItemRectMin();
            ImVec2 p1 = ImGui::GetItemRectMax();
            ImDrawList* draw_list = ImGui::GetWindowDrawList();
            draw_list->PushClipRect(p0, p1);

            float test1 = 25.0;
            float test2 = 70.0;
            float test3 = 20.0;
            float max_val = 5.0;

            std::vector<float> values = {Esu, Esv, Er};

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
    viewer.data(id_3d_left).line_width = 5;
    viewer.data(id_3d_left).line_color = Eigen::RowVector4f(0.0, 0.0, 1.0, 1.0);
    viewer.data(id_3d_right).line_width = 5;
    viewer.data(id_3d_right).line_color = Eigen::RowVector4f(0.0, 0.0, 1.0, 1.0);
    //viewer.data().point_size = 10;
    viewer.core(left_view).orthographic = false;

    
    //viewer.core(left_view).set_rotation_type(igl::opengl::ViewerCore::ROTATION_TYPE_NO_ROTATION); // libigl improvement? can only have one kind of rotation
    viewer.core(right_view).set_rotation_type(igl::opengl::ViewerCore::ROTATION_TYPE_TRACKBALL);
    viewer.core().background_color = Eigen::Vector4f(202.0/255.0, 190.0/255.0, 232.0/255.0, 1.0);
    viewer.launch();
}