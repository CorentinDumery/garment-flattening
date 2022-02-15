/**
 * @author Corentin Dumery
 * @brief Visualize parameterization with/without seam/dart reflectability. 
 * @date 2022-02-15
 * 
 */
#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/readOBJ.h>
#include <iostream>
#include <fstream>

#include "param/multi_patch_param.h"
#include "param/metrics.h"

void loadInputFilesSkirt(std::vector<std::vector<std::vector<std::pair<int, int>>>>& vec_dart_duplicates,
                         std::vector<std::vector<int>>& vec_dart_tips,
                         std::vector<Seam>& seams){
    auto loadDart = [](std::string dart_path, std::vector<std::pair<int, int>>& dupls, int& tip_v){
        std::ifstream input_file(dart_path);
        if (!input_file.is_open()) {
            std::cout << "Could not open the file - '"
                << dart_path << "'" << std::endl;
        }
        int dart_tip, dupl1, dupl2;
        dupls.clear();
        input_file >> tip_v;
        while (input_file >> dupl1) {
            input_file >> dupl2;
            dupls.push_back(std::make_pair(dupl1, dupl2));
        }
        input_file.close();
    };

    std::vector<std::vector<std::pair<int, int>>> darts0, darts1;
    std::vector<int> tips0, tips1;

    std::vector<std::pair<int, int>> dupls;
    int tip_v;
    loadDart("../data/skirt_seams/patch_0_dart_0.txt", dupls, tip_v);
    darts0.push_back(dupls);
    tips0.push_back(tip_v);
    loadDart("../data/skirt_seams/patch_0_dart_1.txt", dupls, tip_v);
    darts0.push_back(dupls);
    tips0.push_back(tip_v);

    vec_dart_duplicates.push_back(darts0);
    vec_dart_tips.push_back(tips0);

    loadDart("../data/skirt_seams/patch_1_dart_0.txt", dupls, tip_v);
    darts1.push_back(dupls);
    tips1.push_back(tip_v);
    loadDart("../data/skirt_seams/patch_1_dart_1.txt", dupls, tip_v);
    darts1.push_back(dupls);
    tips1.push_back(tip_v);

    vec_dart_duplicates.push_back(darts1);
    vec_dart_tips.push_back(tips1);

    auto loadSeam = [](std::string seam_path, Seam& s){
        std::ifstream input_file(seam_path);
        if (!input_file.is_open()) {
            std::cout << "Could not open the file - '"
                << seam_path << "'" << std::endl;
        }
        int patch1, patch2, dupl1, dupl2;
        s.corres.clear();
        input_file >> s.patch1_id;
        input_file >> s.patch2_id;
        while (input_file >> dupl1) {
            input_file >> dupl2;
            s.corres.push_back(std::make_pair(dupl1, dupl2));
        }
        input_file.close();
    };

    seams.clear();
    Seam s;
    loadSeam("../data/skirt_seams/seam0.txt", s);
    seams.push_back(s);
    loadSeam("../data/skirt_seams/seam1.txt", s);
    seams.push_back(s);
}

int main(int argc, char *argv[]){

    Eigen::MatrixXd V1, V2, V1_out, V2_out, V1_dart, V2_dart, V1_seam, V2_seam, V1_both, V2_both;
    Eigen::MatrixXi F1, F2;
    std::vector<Eigen::MatrixXd> out_vec;

    igl::readOBJ("../data/skirt_seams/skirt_0.obj", V1, F1);
    igl::readOBJ("../data/skirt_seams/skirt_1.obj", V2, F2);

    std::vector<std::vector<std::vector<std::pair<int, int>>>> vec_dart_duplicates;
    std::vector<std::vector<int>> vec_dart_tips;
    std::vector<Seam> seams;
    loadInputFilesSkirt(vec_dart_duplicates, vec_dart_tips, seams);

    finalParamMultiPatch({V1, V2}, {F1, F2}, 
                        {{}, {}}, // dart dupl
                        {{}, {}}, // dart tips
                        {}, // seams
                        out_vec);
    
    V1_out = out_vec[0];
    V2_out = out_vec[1];

    finalParamMultiPatch({V1, V2}, {F1, F2}, 
                        vec_dart_duplicates, 
                        vec_dart_tips,
                        {}, // seams
                        out_vec);

    
    V1_dart = out_vec[0];
    V2_dart = out_vec[1];

    finalParamMultiPatch({V1, V2}, {F1, F2}, 
                         {{}, {}}, // dart dupl
                         {{}, {}}, // dart tips
                         seams, // seams
                         out_vec);

    V1_seam = out_vec[0];
    V2_seam = out_vec[1];

    finalParamMultiPatch({V1, V2}, {F1, F2}, 
                        vec_dart_duplicates, 
                        vec_dart_tips,
                        seams, // seams
                        out_vec);
    
    V1_both = out_vec[0];
    V2_both = out_vec[1];

    V1_out.col(0) = V1_out.col(0).array() + 25.0;
    V1_dart.col(0) = V1_dart.col(0).array() + 25.0;
    V1_seam.col(0) = V1_seam.col(0).array() + 25.0;
    V1_both.col(0) = V1_both.col(0).array() + 25.0;

    // --- VISUALIZATION ---

    bool dart_ref = true;
    bool seam_ref = true;
    bool show_seam_corres = false;

    igl::opengl::glfw::Viewer viewer;
    viewer.append_mesh();
    viewer.append_mesh();
    int mesh1_id = viewer.data_list[0].id;
    int mesh2_id = viewer.data_list[1].id;
    int hud_id = viewer.data_list[2].id;

    igl::opengl::glfw::imgui::ImGuiMenu menu;
    menu.callback_draw_viewer_window = []() {};
    viewer.plugins.push_back(&menu);

    auto updateViz = [&](){
        viewer.data().clear_edges();
        viewer.data().clear_points();

        Eigen::MatrixXd sel_V1, sel_V2;
        Eigen::MatrixXi sel_F;
        if (!dart_ref && !seam_ref){ 
            sel_V1 = V1_out;
            sel_V2 = V2_out;
        }
        else if (dart_ref && !seam_ref){  
            sel_V1 = V1_dart;
            sel_V2 = V2_dart;
        }
        else if (!dart_ref && seam_ref){ 
            sel_V1 = V1_seam;
            sel_V2 = V2_seam;
        }
        else if (dart_ref && seam_ref){  
            sel_V1 = V1_both;
            sel_V2 = V2_both;
        }

        viewer.data(mesh1_id).set_mesh(sel_V1, F1);
        viewer.data(mesh2_id).set_mesh(sel_V2, F2);

        if (show_seam_corres){
            for (Seam s1: seams){
                for (int i=0; i<s1.corres.size(); i++){
                    Eigen::RowVector3d color = Eigen::RowVector3d::Random(); 
                    viewer.data(hud_id).add_points(sel_V1.row(s1.corres[i].first), color);
                    viewer.data(hud_id).add_points(sel_V2.row(s1.corres[i].second), color);
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
        if (ImGui::Begin("Multiparam")) {
            
            if (ImGui::Checkbox("Dart reflectability", &dart_ref)){
                updateViz();
            }
            if (ImGui::Checkbox("Seam reflectability", &seam_ref)){
                updateViz();
            }

            ImGui::Separator();
            if (ImGui::Checkbox("Display seam correspondences", &show_seam_corres)){
                updateViz();
            }   
            ImGui::End();
        }
    };

    updateViz();

    viewer.core().set_rotation_type(igl::opengl::ViewerCore::ROTATION_TYPE_TRACKBALL);
    viewer.core().background_color = Eigen::Vector4f(202.0/255.0, 190.0/255.0, 232.0/255.0, 1.0);
    viewer.launch();
}