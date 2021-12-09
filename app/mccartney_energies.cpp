

#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/readOBJ.h>
#include <igl/png/readPNG.h>

#include "net_param.h"
#include "mccartney.h"

void printMatStats(std::string name, const Eigen::VectorXd& mat){
    std::cout << name << ": " <<  mat.minCoeff() << " -> " << mat.maxCoeff() << " (avg " << mat.mean() << ")"  << std::endl;
}

Eigen::MatrixXd fromVectorToColors(const Eigen::VectorXd& vector){
    Eigen::VectorXd adjusted = vector;
    adjusted = adjusted.array() - adjusted.minCoeff();
    double mean = adjusted.mean();
    adjusted = adjusted.array() / (2 * mean); // center mean at 0.5



    adjusted = adjusted.cwiseMin(1.0);
    adjusted = adjusted.cwiseMax(0.0);

    printMatStats("adjusted", adjusted);

    Eigen::MatrixXd colors = Eigen::MatrixXd::Constant(vector.rows(), 3, 1);
    colors.col(1) = 1.0 - adjusted.array();
    colors.col(2) = 1.0 - adjusted.array();

    return colors;
} 

int main(int argc, char *argv[]){

    Eigen::MatrixXd V_3d, V_2d;
    Eigen::MatrixXi F, F0;

    igl::readOBJ("../data/dress_front_cut.obj", V_3d, F0);
    igl::readOBJ("../data/flat_dress.obj", V_2d, F);

    double scale_f = 2.0;
    V_3d *= scale_f;
    V_2d *= scale_f;

    Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> R,G,B,A;
    igl::png::readPNG("../data/grid.png",R,G,B,A);



    if (F.maxCoeff() != F0.maxCoeff()){
        std::cout << "ERROR: meshes don't match" << std::endl;
        return 0;
    }

    Eigen::VectorXd strain_u(F.rows());
    Eigen::VectorXd strain_v(F.rows());
    Eigen::VectorXd shear(F.rows());

    for (int i=0; i<F.rows(); i++){
        Eigen::MatrixXd V_2di(3,3);
        Eigen::MatrixXd V_3di(3,3);
        for (int j=0; j<3; j++){
            V_2di.row(j) = V_2d.row(F(i,j));
            V_3di.row(j) = V_3d.row(F(i,j));
        }
        double Esu, Esv, Er;
        computeMcCartneyErrors(V_2di, V_3di, Esu, Esv, Er);
        strain_u(i) = Esu;
        strain_v(i) = Esv;
        shear(i) = Er;
    }

    printMatStats("Strain U", strain_u);
    printMatStats("Strain V", strain_v);
    printMatStats("Shear", shear);

    

    // --- VISUALIZATION ---

    igl::opengl::glfw::Viewer viewer;
    viewer.data().set_mesh(V_3d, F);

    viewer.data().set_uv(V_2d);
    viewer.data().show_texture = true;
    viewer.data().set_texture(R,G,B);

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
        if (ImGui::Begin("McCartney energies")) {
            float w = ImGui::GetContentRegionAvailWidth();
            float p = ImGui::GetStyle().FramePadding.x;

            if (ImGui::Button("3D",  ImVec2((w - p) / 2.f, 0))){
                viewer.data().set_mesh(V_3d, F);            
            }
            
            ImGui::SameLine(0, p);
            if (ImGui::Button("2D",  ImVec2((w - p) / 2.f, 0))){
                viewer.data().set_mesh(V_2d, F);
            }
            
            if (ImGui::CollapsingHeader("Colors", ImGuiTreeNodeFlags_DefaultOpen)){
                if (ImGui::Button("Strain U colors", ImVec2(-1,0))){
                    viewer.data().set_colors(fromVectorToColors(strain_u));
                }

                if (ImGui::Button("Strain V colors", ImVec2(-1,0))){
                    viewer.data().set_colors(fromVectorToColors(strain_v));
                }

                if (ImGui::Button("Shear colors", ImVec2(-1,0))){
                    viewer.data().set_colors(fromVectorToColors(shear));
                }

                if (ImGui::Button("White", ImVec2(-1,0))){
                    viewer.data().set_colors(Eigen::RowVector3d(1.0, 1.0, 1.0));
                }
            }

            ImGui::Separator();

            if (ImGui::Button("Viz net", ImVec2(-1,0))){

                // Per triangle net transported from UV to 3D.
                // should be identical to texturing with a grid

                for (int f_id=0; f_id<F.rows(); f_id++){
                    Eigen::MatrixXd V_2di(3,3);
                    Eigen::MatrixXd V_3di(3,3);
                    for (int j=0; j<3; j++){
                        V_2di.row(j) = V_2d.row(F(f_id,j));
                        V_3di.row(j) = V_3d.row(F(f_id,j));
                    }

                    Eigen::MatrixXi F2(1,3);
                    F2 << 0, 1, 2;

                    NetParam net_param(F2, V_3di.cast<float>(), V_2di.cast<float>());
                    net_param.computeFibers();
                    
                    std::vector<std::vector<Eigen::MatrixXd>> fibs = net_param.vizFibers();
                    std::vector<Eigen::MatrixXd> fiber_begs_list = fibs[0];
                    std::vector<Eigen::MatrixXd> fiber_ends_list = fibs[1];

                    // Transport net to other triangle
                    viewer.data().add_edges(transportMatrix(V_2di, F2, fiber_begs_list[0], V_3di), 
                                            transportMatrix(V_2di, F2, fiber_ends_list[0], V_3di), 
                                            Eigen::RowVector3d(0.0, 0.0, 1.0));
                    viewer.data().add_edges(transportMatrix(V_2di, F2, fiber_begs_list[1], V_3di), 
                                            transportMatrix(V_2di, F2, fiber_ends_list[1], V_3di), 
                                            Eigen::RowVector3d(0.0, 0.0, 1.0));
                }
            }
            
            ImGui::End();
        }
    
    };

    //updateViz();

    viewer.data().set_colors(Eigen::RowVector3d(1.0, 1.0, 1.0));
    viewer.data().line_width = 5;
    viewer.data().show_lines = false;
    //viewer.data().point_size = 10;
    //viewer.core().orthographic = true;

    viewer.core().set_rotation_type(igl::opengl::ViewerCore::ROTATION_TYPE_TRACKBALL);
    viewer.core().background_color = Eigen::Vector4f(202.0/255.0, 190.0/255.0, 232.0/255.0, 1.0);
    viewer.launch();
}