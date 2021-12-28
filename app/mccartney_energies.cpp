

#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/readOBJ.h>
#include <igl/png/readPNG.h>
#include <igl/edges.h>

//#define COMP_WITH_NET_PARAM
#ifdef COMP_WITH_NET_PARAM
#include "net_param.h"
#endif

#include "mccartney.h"
#include "procustes.h"
#include "bary_optimizer.h"

igl::opengl::glfw::Viewer viewer; // TODO MOVE

void printMatStats(std::string name, const Eigen::VectorXd& mat){
    std::cout << name << ": " <<  mat.minCoeff() << " -> " << mat.maxCoeff() << " (avg " << mat.mean() << ")"  << std::endl;
}

Eigen::MatrixXd fromVectorToColors(const Eigen::VectorXd& vector){
    Eigen::VectorXd adjusted = vector;
    adjusted = adjusted.array() - adjusted.minCoeff();
    double mean = adjusted.mean();
    //adjusted = adjusted.array() / (2 * mean); // center mean at 0.5
    
    adjusted = adjusted.array() * 100.0;

    printMatStats("vector", vector);


    adjusted = adjusted.cwiseMin(1.0);
    adjusted = adjusted.cwiseMax(0.0);

    printMatStats("adjusted", adjusted);

    Eigen::MatrixXd colors = Eigen::MatrixXd::Constant(vector.rows(), 3, 1);
    colors.col(1) = 1.0 - adjusted.array();
    colors.col(2) = 1.0 - adjusted.array();

    return colors;
}

Eigen::VectorXd paramLocalGlobal(const Eigen::MatrixXd& V_3d, const Eigen::MatrixXi& F,
                      Eigen::MatrixXd& V_2d, int viz_axis_error=2){

    
    
    // Compute ideal rotations per triangle
    /*for (int f_id=0; f_id<F.rows(); f_id++){
        Eigen::MatrixXd p1 = makeTriPoints(V_2d, F, f_id);
        Eigen::MatrixXd p2_temp(3,3), p2(3,3);
        p2_temp = makeTriPoints(V_3d, F, f_id);
        p2 = move3Dto2D(p2_temp);
        
        Eigen::MatrixXd R_est;
        Eigen::VectorXd T_est;
        procustes(p1, p2, R_est, T_est);

        std::cout << R_est << std::endl << std::endl;
    }*/

    Eigen::VectorXd energy_u, energy_v;
    energy_u = Eigen::VectorXd::Zero(F.rows());
    energy_v = Eigen::VectorXd::Zero(F.rows());

    for (int f_id=0; f_id<F.rows(); f_id++){
        Eigen::MatrixXd p1 = makeTriPoints(V_2d, F, f_id);
        Eigen::MatrixXd p2_temp(3,3), p2(3,3);
        p2_temp = makeTriPoints(V_3d, F, f_id);
        p2 = move3Dto2D(p2_temp);
        

        Eigen::MatrixXd R_est;
        Eigen::VectorXd T_est;
        procustes(p1, p2, R_est, T_est);

        /*Eigen::MatrixXd p2_r;
        p2_r = p2.transpose();
        p2_r = p2_r.colwise() - T_est;
        p2_r = (R_est.transpose() * p2_r);*/
        //p2_r.transpose();

        //p2_r = p2; // TODO REMOVE !!!!!!!!!!!!!!!

        Eigen::MatrixXd p2_rt, p2_r;
        Eigen::MatrixXd p2t = p2.transpose();
        p2_rt = p2t.colwise() - T_est;
        p2_rt = (R_est.transpose() * p2_rt);
        p2_r = p2_rt.transpose();
        
        //viewer.data().add_points(p1, Eigen::RowVector3d(1.0, 1.0, 0.0));
        //viewer.data().add_points(p2_r, Eigen::RowVector3d(0.0, 1.0, 1.0));

        double ABu = (p1.row(1) - p1.row(0))(0);
        double ApBpu = (p2_r.row(1) - p2_r.row(0))(0);
        double ACu = (p1.row(2) - p1.row(0))(0);
        double ApCpu = (p2_r.row(2) - p2_r.row(0))(0);

        double ABv = (p1.row(1) - p1.row(0))(1);
        double ApBpv = (p2_r.row(1) - p2_r.row(0))(1);
        double ACv = (p1.row(2) - p1.row(0))(1);
        double ApCpv = (p2_r.row(2) - p2_r.row(0))(1);

        double Eu = std::pow(ABu - ApBpu, 2) 
                  + std::pow(ACu - ApCpu, 2);
        double Ev = std::pow(ABv - ApBpv, 2) 
                  + std::pow(ACv - ApCpv, 2);

        energy_v(f_id) = Eu;
        energy_u(f_id) = Ev;

        //break;
    }

    if (viz_axis_error == 0){
        return energy_u;
    }
    else if (viz_axis_error == 1){
        return energy_v;
    }
    else {
        return energy_u + energy_v;
    }    
}

int main(int argc, char *argv[]){

    //Eigen::setNbThreads(1);
    std::cout << "Eigen is using " << Eigen::nbThreads() << " threads." << std::endl;

    Eigen::MatrixXd V_3d, V_2d, V_2di;
    Eigen::MatrixXi F, F0;

    /*
    igl::readOBJ("../data/dress_front_cut.obj", V_3d, F0);
    igl::readOBJ("../data/flat_dress.obj", V_2d, F);
    //*/

    //*
    igl::readOBJ("../data/semisphere_uncut.obj", V_3d, F0);
    igl::readOBJ("../data/semisphere_uncut_flat.obj", V_2d, F);
    //*/


    /*
    igl::readOBJ("../data/semisphere0.obj", V_3d, F0);
    igl::readOBJ("../data/semisphere0_flat_halfbad.obj", V_2d, F);//*/


    /*
    igl::readOBJ("../data/cross.obj", V_3d, F0);
    igl::readOBJ("../data/cross_flat_bad.obj", V_2d, F);
    Eigen::VectorXd temp = V_3d.col(2);
    V_3d.col(2) = V_3d.col(1);
    V_3d.col(1) = temp;
    temp = V_2d.col(2);
    V_2d.col(2) = V_2d.col(1);
    V_2d.col(1) = temp; //*/


    /*
    V_2d.resize(4, 3);
    V_3d.resize(4, 3);
    F.resize(2, 3);

    V_2d <<  0,   0, 0,
           1.0,   0, 0,
             0, 1.0, 0,
           1.0, 1.5, 0;

    V_3d <<  0,   0, 0,
           1.0,   0, 0,
             0, 1.0, 0,
           1.0, 1.0, 0;

    F << 0, 1, 2,
         1, 3, 2;

    F0 = F;
    //*/

    /*
    V_2d.resize(4, 3);
    V_3d.resize(4, 3);
    F.resize(2, 3);

    double d = std::sqrt(2.0)/2.0;
    V_2d <<  0,   0, 0,
           1.0,   0, 0,
             d,     d, 0,
           1.0 + d, d, 0;

    V_3d <<  0,   0, 0,
           1.0,   0, 0,
             0, 1.0, 0,
           1.0, 1.0, 0;

    F << 0, 1, 2,
         1, 3, 2;

    F0 = F;
    //*/

    /*
    V_2d.resize(4, 3);
    V_3d.resize(4, 3);
    F.resize(2, 3);

    double d = std::sqrt(2.0)/2.0;
    V_3d <<  0,   0, 0,
           1.0,   0, 0,
             d,     d, 0,
           1.0 + d, d, 0;

    V_2d <<  0,   0, 0,
           1.0,   0, 0,
             0, 1.0, 0,
           1.0, 1.0, 0;

    F << 0, 1, 2,
         1, 3, 2;

    F0 = F;
    //*/

    double scale_f = 1.0;
    V_3d *= scale_f;
    V_2d *= scale_f;

    V_2di = V_2d;

    //Eigen::MatrixXi E;
    //igl::edges(F, E);

    std::cout << "F.rows(): " << F.rows() << std::endl;

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

    BaryOptimizer bo;    

    // --- VISUALIZATION ---

    bool update_v2d = false;
    
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

            make_checkbox("Show mesh", viewer.data().show_lines);
            make_checkbox("Show texture", viewer.data().show_texture);
            
            if (ImGui::CollapsingHeader("Colors", ImGuiTreeNodeFlags_DefaultOpen)){
                /*if (ImGui::Button("Strain U colors", ImVec2(-1,0))){
                    viewer.data().set_colors(fromVectorToColors(strain_u));
                }

                if (ImGui::Button("Strain V colors", ImVec2(-1,0))){
                    viewer.data().set_colors(fromVectorToColors(strain_v));
                }

                if (ImGui::Button("Shear colors", ImVec2(-1,0))){
                    viewer.data().set_colors(fromVectorToColors(shear));
                }*/

                if (ImGui::Button("White", ImVec2(-1,0))){
                    viewer.data().set_colors(Eigen::RowVector3d(1.0, 1.0, 1.0));
                }

                if (ImGui::Button("My Stretch U", ImVec2(-1, 0))){
                    Eigen::VectorXd E = paramLocalGlobal(V_3d, F, V_2d, 0);
                    viewer.data().set_colors(fromVectorToColors(E));
                }

                if (ImGui::Button("My Stretch V", ImVec2(-1, 0))){
                    Eigen::VectorXd E = paramLocalGlobal(V_3d, F, V_2d, 1);
                    viewer.data().set_colors(fromVectorToColors(E));
                }

                if (ImGui::Button("Sum", ImVec2(-1, 0))){
                    Eigen::VectorXd E = paramLocalGlobal(V_3d, F, V_2d, 2);
                    viewer.data().set_colors(fromVectorToColors(E));
                }

                ImGui::SliderFloat("Light factor", &viewer.core().lighting_factor, 0.0f, 5.0f, "%.3f");
                
            }

            ImGui::Separator();

            #ifdef COMP_WITH_NET_PARAM
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
            #endif

            
            ImGui::Separator();

            if (ImGui::Button("Local global test", ImVec2(-1, 0))){
                for (int i=0; i<30; i++){
                    V_2d = bo.localGlobal(V_2d, V_3d, F);
                }
                if (update_v2d) viewer.data().set_mesh(V_2d, F);
                viewer.data().set_uv(V_2d);
                //Eigen::VectorXd E = paramLocalGlobal(V_3d, F, V_2d, 0);
                //viewer.data().set_colors(fromVectorToColors(E));
            }

            float strech_f = bo.stretch_coeff_;
            if (ImGui::SliderFloat("Stretch factor", &strech_f, 0.0f, 50.0f, "%.3f")){
                V_2d = V_2di;
                bo.stretch_coeff_ = strech_f;
                for (int i=0; i<1; i++){
                    V_2d = bo.localGlobal(V_2d, V_3d, F);
                }
                if (update_v2d) viewer.data().set_mesh(V_2d, F);
                viewer.data().set_uv(V_2d);
            }

            float angle_f = bo.angle_coeff_;
            if (ImGui::SliderFloat("Angles factor", &angle_f, 0.0f, 15.0f, "%.3f")){
                V_2d = V_2di;
                bo.angle_coeff_ = angle_f;
                for (int i=0; i<1; i++){
                    V_2d = bo.localGlobal(V_2d, V_3d, F);
                }
                if (update_v2d) viewer.data().set_mesh(V_2d, F);
                viewer.data().set_uv(V_2d);
            }

            float edges_f = bo.edges_coeff_;
            if (ImGui::SliderFloat("Edges factor", &edges_f, 0.0f, 10.0f, "%.3f")){
                V_2d = V_2di;
                bo.edges_coeff_ = edges_f;
                for (int i=0; i<1; i++){
                    V_2d = bo.localGlobal(V_2d, V_3d, F);
                }
                if (update_v2d) viewer.data().set_mesh(V_2d, F);
                viewer.data().set_uv(V_2d);
            }

            ImGui::Checkbox("Update V_2d", &update_v2d);

            ImGui::End();
        }
    
    };

    //updateViz();

    viewer.data().set_colors(Eigen::RowVector3d(1.0, 1.0, 1.0));
    viewer.data().line_width = 3;
    viewer.data().show_lines = false;
    //viewer.data().point_size = 10;
    //viewer.core().orthographic = true;

    viewer.core().set_rotation_type(igl::opengl::ViewerCore::ROTATION_TYPE_TRACKBALL);
    viewer.core().background_color = Eigen::Vector4f(202.0/255.0, 190.0/255.0, 232.0/255.0, 1.0);
    viewer.launch();
}