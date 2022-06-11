/**
 * @author Corentin Dumery 
 * @brief Example and visualization of multiple poses parameterization on a jumpsuit. 
 * @date 2022-02-10
 * 
 */

#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/readOBJ.h>
#include <igl/avg_edge_length.h>
#include <igl/png/readPNG.h>

#include "param/multiple_poses_param.h"

int main(int argc, char *argv[]){

    Eigen::MatrixXd V1, V2, V3, Vf_solo, Vf_duo, Vf_trio;
    Eigen::MatrixXi F1, F2, F3;

    igl::readOBJ("../data/jumpsuit_multipose/front1.obj", V1, F1);
    igl::readOBJ("../data/jumpsuit_multipose/front2.obj", V2, F2);
    igl::readOBJ("../data/jumpsuit_multipose/front3.obj", V3, F3);

    Vf_solo = multiplePosesParam({V1}, F1, 0.0, {}, {});
    Vf_duo = multiplePosesParam({V1, V2}, F1, 0.0, {}, {});
    Vf_trio = multiplePosesParam({V1, V2, V3}, F1, 0.0, {}, {});

    // --- Show PNG image ---

    Eigen::MatrixXd V_im(4,3);
    V_im <<
        -0.5,-0.5, 0,
         0.5,-0.5, 0,
         0.5, 0.5, 0,
        -0.5, 0.5, 0;
    
    V_im.col(0) *= 1.75925925926;
    
    Eigen::MatrixXi F_im(2,3);
    F_im << 0,1,2,
            2,3,0;

    Eigen::MatrixXd UV_im(4,2);
    UV_im <<
        0,0,
        1,0,
        1,1,
        0,1;

    double l = igl::avg_edge_length(V1, F1);
    V_im *= 75 * l;
    V_im.col(0) = V_im.col(0).array() - 75 * l;

    Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> R(1080,1080);
    Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> G(1080,1080);
    Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> B(1080,1080);
    Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> A(1080,1080);
    igl::png::readPNG("../data/jumpsuit_multipose/pic_solo.png",R,G,B,A);

    // --- VISUALIZATION ---

    igl::opengl::glfw::Viewer viewer;
    viewer.append_mesh();
    int mesh1_id = viewer.data_list[1].id;
    int im_id = viewer.data_list[0].id;

    viewer.data(mesh1_id).set_mesh(Vf_solo, F1);
    viewer.data(im_id).set_mesh(V_im, F_im);
    viewer.data(im_id).set_uv(UV_im);
    viewer.data(im_id).set_colors(Eigen::RowVector3d(1.0, 1.0, 1.0));


    viewer.data(im_id).show_texture = true;
    viewer.data(im_id).set_texture(R,G,B,A);
    viewer.data(im_id).show_lines = false;



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
        if (ImGui::Begin("Multiple poses")) {
            
            if (ImGui::Button("Solo", ImVec2(-1, 0))){
                viewer.data(mesh1_id).set_mesh(Vf_solo, F1);
                igl::png::readPNG("../data/jumpsuit_multipose/pic_solo.png",R,G,B,A);
                viewer.data(im_id).set_texture(R,G,B,A);
            }
            if (ImGui::Button("Duo", ImVec2(-1, 0))){
                viewer.data(mesh1_id).set_mesh(Vf_duo, F1);
                igl::png::readPNG("../data/jumpsuit_multipose/pic_duo.png",R,G,B,A);
                viewer.data(im_id).set_texture(R,G,B,A);
            }
            if (ImGui::Button("Trio", ImVec2(-1, 0))){
                viewer.data(mesh1_id).set_mesh(Vf_trio, F1);
                igl::png::readPNG("../data/jumpsuit_multipose/pic_trio.png",R,G,B,A);
                viewer.data(im_id).set_texture(R,G,B,A);
            }

            ImGui::Text("Multiple pose parameterization, taking each pose into");
            ImGui::Text("account. The third input is not symmetric, leading to ");
            ImGui::Text("a non symmetric parameterization.");
            
            ImGui::End();
        }
    
    };

    updateViz();

    viewer.core().camera_translation = Eigen::Vector3f(8, 0, 0);
    viewer.core().camera_zoom = 1.3;
    viewer.core().set_rotation_type(igl::opengl::ViewerCore::ROTATION_TYPE_NO_ROTATION);
    //viewer.core().background_color = Eigen::Vector4f(202.0/255.0, 190.0/255.0, 232.0/255.0, 1.0);
    viewer.core().background_color = Eigen::Vector4f(1.0, 1.0, 1.0, 1.0);
    viewer.launch();
}
