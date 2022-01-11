
#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/readOBJ.h>
#include <igl/arap.h>
#include <igl/cut_mesh.h>
#include <imgui.h>
#include <igl/per_vertex_normals.h>
#include <math.h>
#include <igl/doublearea.h>
#include <igl/jet.h>
#include <igl/internal_angles.h>
#include <igl/writeOBJ.h>

//#define USE_IMPLOT
#ifdef USE_IMPLOT
#include "implot.h"
#endif

#include "param/dart.h"
#include "param/param_utils.h"

// SCAF
#include <igl/triangle/scaf.h>
#include <igl/arap.h>
#include <igl/boundary_loop.h>
#include <igl/harmonic.h>
#include <igl/map_vertices_to_circle.h>
#include <igl/MappingEnergyType.h>
#include <igl/doublearea.h>
#include <igl/PI.h>
#include <igl/flipped_triangles.h>
#include <igl/topological_hole_fill.h>

#include "draw_colormap.h"

Eigen::MatrixXi cutsMatrix(const std::vector<int>& selected, 
                           const Eigen::MatrixXi& F){
    Eigen::MatrixXi cuts = Eigen::MatrixXi::Zero(F.rows(), 3);

    for (int i=0; i<F.rows(); i++){
        for (int j=0; j<selected.size(); j++){
            for (int k=0; k<selected.size(); k++){
                if (k==j) continue;
                for (int l=0; l<3; l++){
                    if (selected[j] == F(i,l)){
                        int m = (l+1)%3;
                        if (selected[k] == F(i,m)){
                            cuts(i, l) = 1;
                        } 
                    }
                }
            }
        }
    }
    return cuts;
}

std::vector<std::vector<int>> computeConnectivity(const std::vector<int>& selected,
                                                  const Eigen::MatrixXi& F){
    
    std::vector<std::vector<int>> neighbors;
    for (int i=0; i<selected.size(); i++){
        neighbors.push_back({});
    }
    for (int i=0; i<F.rows(); i++){
        for (int j=0; j<3; j++){
            for (int k=0; k<selected.size(); k++){
                if (F(i,j) == selected[k]){
                    for (int l=0; l<selected.size(); l++){
                        if (l==k) continue;
                        if (selected[l] == F(i,(j+1)%3)) neighbors[k].push_back(l);
                        if (selected[l] == F(i,(j+2)%3)) neighbors[k].push_back(l);
                    }
                }
            }    
        }
    }

    for (int i=0; i<neighbors.size(); i++){
        std::sort(neighbors[i].begin(), neighbors[i].end());
        neighbors[i].erase(std::unique(neighbors[i].begin(), neighbors[i].end()), neighbors[i].end());

    }

    /*std::cout << "neighbors: " << std::endl;
    for (auto i: neighbors){
        std::cout << "x: ";
        for (int j: i){
            std::cout << selected[j] << " ";
        }
        std::cout << std::endl;
    }*/

    return neighbors;
}

std::vector<int> restoreCutOrder(const std::vector<int>& selected,
                                 const Eigen::MatrixXi& F){

    std::vector<std::vector<int>> neigbhors = computeConnectivity(selected, F);

    int curr_v = -1;
    for (int i=0; i<neigbhors.size(); i++){
        if (neigbhors[i].size()==1){
            curr_v = i;
            break;
        }
    }

    std::cout << "Cut starts in " << curr_v << std::endl;

    if (curr_v < 0){
        std::cout << "ERROR: cut starting point not found" << std::endl;
        return {};
    }
    
    std::vector<int> new_selected = {selected[curr_v]};
    int next_v = neigbhors[curr_v][0];
    while (neigbhors[next_v].size() > 1){
        if (neigbhors[next_v][0] == curr_v){
            curr_v = next_v;
            next_v = neigbhors[next_v][1];
        }
        else if (neigbhors[next_v][1] == curr_v){
            curr_v = next_v;
            next_v = neigbhors[next_v][0];
        }
        else {
            std::cout << "Failed to order cutting path." << std::endl;
            break;
        }
        new_selected.push_back(selected[curr_v]);
    }
    new_selected.push_back(selected[next_v]);

    /*std::cout << "new_selected: ";
    for (auto i: new_selected){
        std::cout << i << " ";
    }
    std::cout << std::endl;*/
    return new_selected;
}

std::vector<int> identifyCut(const Eigen::VectorXi& corres, 
                             const Eigen::MatrixXi& F,
                             const std::vector<int>& selected){ // assumes pre-cut ids are still used in cut (on either side)

    std::vector<int> selected_dupl = selected;
    for (int i=0; i<corres.size(); i++){
        if (corres[i] != i){
            for (int j=0; j<selected.size(); j++){
                if (corres[i] == selected[j]){
                    selected_dupl.push_back(i);
                }
            }
        }
    }

    std::cout << "Restoring cut order..." << std::endl;
    std::vector<int> cut = restoreCutOrder(selected_dupl, F);

    /*std::cout << "cut: ";
    for (auto i: cut){
        std::cout << i << " ";
    }
    std::cout << std::endl;*/  

    return cut;
}

void symmetrizeDarts(Eigen::MatrixXd& V, const std::vector<SimpleDart>& simple_darts){

}

int main(int argc, char *argv[])
{

    std::string input_name;
    std::vector<std::vector<int>> darts;
    std::vector<int> reference_orientation = {0, 1};

    int input_id = 0;

    if (argc > 1) {
        input_name = argv[1];
        input_id = -1;
    }

    if (input_id == 0){
        input_name = "../data/dress_front.obj";
        std::vector<int> dart1 = {123, 132, 133, 316, 370, 372, 446, 470, 593, 594, 597, 746, 1033, 1037, 1191};
        std::vector<int> dart2 = {786, 787, 792, 906, 907}; // small dart in necklace location
        std::vector<int> dart3 = {71, 74, 265, 423, 424, 442, 443}; 
        std::vector<int> dart4 = {340, 341, 371, 372, 446, 746, 773, 1191};

        std::vector<int> diamond1 = {221, 222, 223, 316, 470, 593, 594, 597, 746, 1033, 1191}; // diamond below chest

        darts = {dart2, dart3, dart4};
        //darts = {dart2, diamond1};
        //darts = {dart3};

        reference_orientation = {792, 1700};
    }

    else if (input_id == 1){
        input_name = "../data/dress2_front.obj";

        std::vector<int> dart_rightarm = {360, 361, 652, 5050, 5052, 5367, 5369, 5552, 5600, 5601, 5804, 6092, 6093, 6153, 6161, 6622, 6648};
        std::vector<int> dart_rightarm_other = {5068, 5069, 5149, 5672, 5820, 5821, 5932, 6064, 6065};
        darts = {dart_rightarm};
        darts = {dart_rightarm_other};
        reference_orientation = {1169, 5168};
    }

    else if (input_id == 2){
        input_name = "../data/semisphere.obj";

        std::vector<int> dart = {448, 451, 455, 460, 466, 473, 481, 490, 500, 511, 523, 536, 549, 957, 958, 967, 1186, 1188, 1191, 1195};
        darts = {dart};
        reference_orientation = {429, 1017}; // one on side one on top
    }

    else if (input_id == 3){
        input_name = "../data/semisphere_straight.obj";

        std::vector<int> dart = {5, 73, 211, 212, 213, 359, 367, 1145, 1241, 1242, 1243, 1244, 1245, 1246, 1247, 1248, 
                                 1249, 1250, 1251, 1252, 1253, 1254, 1255, 1256, 1257, 1258, 1259, 1260, 1261, 1262, 
                                 1263, 1264, 1265, 1266, 1267, 1268};
        darts = {dart};
        reference_orientation = {5, 1038}; // one on side one on top
    }

    //selected = restoreCutOrder(selected, F);
    

    Eigen::MatrixXd V, V_3D, V_uv;
    Eigen::MatrixXi F;
    igl::readOBJ(input_name, V, F);
    std::cout << "Mesh size:" << std::endl;
    std::cout << "V: " << V.rows() << " " << V.cols() << std::endl;
    std::cout << "V min max: " << V.minCoeff() << " -> " << V.maxCoeff() << std::endl;
    std::cout << "F: " << F.rows() << " " << F.cols() << std::endl;

    V = V.array() * 20.0 / (V.maxCoeff() - V.minCoeff());
    std::cout << "V min max: " << V.minCoeff() << " -> " << V.maxCoeff() << std::endl;

    // --- Find selected edges --- 

    std::vector<std::vector<int>> ordered_cuts;

    for (auto selection: darts){
        Eigen::MatrixXi cuts = cutsMatrix(selection, F);

        std::cout << "Cutting..." << std::endl;
        Eigen::VectorXi I; // vertex correspondence pre/post cut
        igl::cut_mesh(V, F, cuts, I);
        std::vector<int> ordered_cut = identifyCut(I, F, selection);

        // TODO here we PRAY that one cut doesn't change ids for future/previous cuts....
        // someone should make sure of that

        ordered_cuts.push_back(ordered_cut);
    }

    std::vector<SimpleDart> simple_darts;
    for (int i=0; i<ordered_cuts.size(); i++){
        std::vector<int> cut = ordered_cuts[i]; 
        if (cut.size() % 2 == 0) continue;
        SimpleDart sd(cut);
        sd.print();
        simple_darts.push_back(sd);
    }

    V_3D = V;

    igl::writeOBJ("../data/cut_mesh.obj", V_3D, F);
    std::cout << "Cut." << std::endl;

    // --- SCAF --- 

    Eigen::MatrixXd V_uv_2d;
    bool enable_scaf = true;
    std::vector<double> flattening_energies;

    if (enable_scaf){
        Eigen::MatrixXd bnd_uv, uv_init;
        igl::triangle::SCAFData scaf_data;

        Eigen::VectorXd M;
        igl::doublearea(V, F, M);
        std::vector<std::vector<int>> all_bnds;
        igl::boundary_loop(F, all_bnds);

        std::cout << "Boundaries: " << all_bnds.size() << std::endl;

        // Heuristic primary boundary choice: longest
        auto primary_bnd = std::max_element(all_bnds.begin(), all_bnds.end(), [](const std::vector<int> &a, const std::vector<int> &b) { return a.size()<b.size(); });

        Eigen::VectorXi bnd = Eigen::Map<Eigen::VectorXi>(primary_bnd->data(), primary_bnd->size());

        igl::map_vertices_to_circle(V, bnd, bnd_uv);
        bnd_uv *= sqrt(M.sum() / (2 * igl::PI));
        if (all_bnds.size() == 1)
        {
            if (bnd.rows() == V.rows()) // case: all vertex on boundary
            {
                uv_init.resize(V.rows(), 2);
                for (int i = 0; i < bnd.rows(); i++)
                uv_init.row(bnd(i)) = bnd_uv.row(i);
            }
            else
            {
                igl::harmonic(V, F, bnd, bnd_uv, 1, uv_init);
                if (igl::flipped_triangles(uv_init, F).size() != 0)
                igl::harmonic(F, bnd, bnd_uv, 1, uv_init); // fallback uniform laplacian
            }
        }
        else
        {
            // if there is a hole, fill it and erase additional vertices.
            all_bnds.erase(primary_bnd);
            Eigen::MatrixXi F_filled;
            igl::topological_hole_fill(F, bnd, all_bnds, F_filled);
            igl::harmonic(F_filled, bnd, bnd_uv ,1, uv_init);
            uv_init.conservativeResize(V.rows(), 2);
        }

        Eigen::VectorXi b; Eigen::MatrixXd bc;

        igl::triangle::scaf_precompute(V, F, uv_init, scaf_data, igl::MappingEnergyType::SYMMETRIC_DIRICHLET, b, bc, 0);

        for (int i=0; i<10; i++){
            igl::triangle::scaf_solve(scaf_data, 1);
            flattening_energies.push_back(scaf_data.energy);
        }
        V_uv_2d = scaf_data.w_uv.topRows(V.rows());
        V_uv = Eigen::MatrixXd::Zero(V.rows(), 3);
        V_uv.col(0) = V_uv_2d.col(0);
        V_uv.col(1) = V_uv_2d.col(1);

        
    }

    Eigen::Matrix3d rot = computeRotation(V_3D.row(reference_orientation[0]) - V_3D.row(reference_orientation[1]),
                                          V_uv.row(reference_orientation[0]) - V_uv.row(reference_orientation[1]));
    Eigen::RowVector3d translation;
    
    bool match_rotation = true;
    if (match_rotation){
        V_3D = (rot * V_3D.transpose()).transpose();
        translation = V_uv.row(reference_orientation[0]) - V_3D.row(reference_orientation[0]);
        V_3D = V_3D.rowwise() + translation;
        V = V_3D;
    }

    igl::writeOBJ("../data/cut_mesh_flat.obj", V_uv, F);

    Eigen::MatrixXd N, N_uv, N_3D;
    igl::per_vertex_normals(V_uv, F, N_uv);
    igl::per_vertex_normals(V_3D, F, N_3D);
    N = N_uv;

    Eigen::VectorXd A_uv, A_3D;
    igl::doublearea(V_uv, F, A_uv);
    igl::doublearea(V_3D, F, A_3D);

    Eigen::VectorXd area_dist = A_uv.array() / A_3D.array();

    Eigen::MatrixXd area_dist_colors(area_dist.rows(), 3);
    double min_area_ratio = 0.75;
    double max_area_ratio = 1.25;
    igl::jet(area_dist, min_area_ratio, max_area_ratio, area_dist_colors);

    Eigen::MatrixXd Ang_uv, Ang_3D;
    igl::internal_angles(V_uv, F, Ang_uv);
    igl::internal_angles(V_3D, F, Ang_3D);

    Eigen::VectorXd ang_dist = ((Ang_uv.array() - Ang_3D.array()) / Ang_3D.array()).abs();
    Eigen::MatrixXd ang_dist_colors(ang_dist.rows(), 3);
    igl::colormap(igl::COLOR_MAP_TYPE_MAGMA, ang_dist, 0.0, 0.1, ang_dist_colors);


    // --- Viewer ---

    float anim_time = 3.14/2.0;
    bool display_x_pattern = true;

    enum DIST_TYPE {AREA_DIST, ANGLE_DIST};
    DIST_TYPE dist_type = AREA_DIST;

    igl::opengl::glfw::Viewer viewer;
    viewer.data().set_mesh(V, F);
    viewer.data().set_uv(V_uv_2d);

    igl::opengl::glfw::imgui::ImGuiMenu menu;
    menu.callback_draw_viewer_window = []() {};
    viewer.plugins.push_back(&menu);

    auto updateViz = [&](){
        viewer.data().clear_edges();
        viewer.data().clear_points();
        //viewer.data().clear();
        viewer.data().set_mesh(V, F);
        //viewer.data().set_normals(N); // needs a data().clear() ?? (so basically doesn't work?)

        if (!viewer.data().show_texture){
            if (dist_type == AREA_DIST) viewer.data().set_colors(area_dist_colors);
            if (dist_type == ANGLE_DIST) viewer.data().set_colors(ang_dist_colors);
        }

        // --- VIZ CUT ---
        
        // cross X pattern
        for (auto dart: simple_darts)
            if (display_x_pattern){
                Eigen::MatrixXd edge_begs, edge_ends;
                dart.xPattern(V, edge_begs, edge_ends);
                viewer.data().add_edges(edge_begs, edge_ends, Eigen::RowVector3d(0., 0., 0.));
            }
        
        //std::cout << "# of darts: " << simple_darts.size() << std::endl;
        for (auto dart: simple_darts){
            if (dart.size() < 3) continue;
            Eigen::RowVector3d axis = dart.computeSymmetryAxis(V_uv_2d);
        
            Eigen::MatrixXd edge_begs(1, 3);
            Eigen::MatrixXd edge_ends(1, 3);
            edge_begs.row(0) = V_uv_2d.row(dart.tip_id());
            edge_ends.row(0) = V_uv_2d.row(dart.tip_id());
            edge_begs(2) = 0;
            edge_ends(2) = 0;
            edge_ends(0) = edge_ends(0) + axis(0);
            edge_ends(1) = edge_ends(1) + axis(1);

            if (false && match_rotation) {
                edge_begs = (rot * edge_begs.transpose()).transpose();
                edge_ends = (rot * edge_ends.transpose()).transpose();
                edge_begs = edge_begs + translation;
                edge_ends = edge_ends + translation;
            }
            //points = points.rowwise() + V_uv_2d.row(ordered_cut[semi_cut]);

            /*Eigen::MatrixXd points3d = Eigen::MatrixXd::Zero(points.rows(), 3);
            points3d.col(0) = points.col(0);
            points3d.col(1) = points.col(1);
            points3d = (rot * points3d.transpose()).transpose();
            points3d = points3d + translation;*/


            viewer.data().add_edges(edge_begs, edge_ends, Eigen::RowVector3d(0., 1., 1.));
            //viewer.data().add_points(points3d, Eigen::RowVector3d(1., 1., 1.));

            //Eigen::MatrixXd points = dart.getMiddlePoints(V);
            //points = (rot * points.transpose()).transpose();
            //points = points + translation;
            //viewer.data().add_points(points, Eigen::RowVector3d(1., 0., 0.));

            Eigen::MatrixXd sym_points_2d = dart.getSymmetricPoints(V_uv_2d, axis);
            Eigen::MatrixXd sym_points = Eigen::MatrixXd::Zero(sym_points_2d.rows(), 3);
            sym_points.col(0) = sym_points_2d.col(0);
            sym_points.col(1) = sym_points_2d.col(1);
            if (false && match_rotation) {
                sym_points = (rot * sym_points.transpose()).transpose();
                sym_points = sym_points + translation;
            }
            viewer.data().add_points(sym_points, Eigen::RowVector3d(1., 1., 0.));

        }
    };

    //helper function for menu
    auto make_checkbox = [&](const char *label, unsigned int &option) {
        return ImGui::Checkbox(
            label,
            [&]() { return viewer.core().is_set(option); },
            [&](bool value) { return viewer.core().set(option, value); });
    };

    bool morphing = false;
    float uv_factor = 1.0;
    float time_increment = 0.01;

    menu.callback_draw_custom_window = [&]() {
        bool show = true;
        ImGui::SetNextWindowPos(ImVec2(10, 10), ImGuiCond_FirstUseEver);
        ImGui::SetNextWindowSize(ImVec2(350, -1), ImGuiCond_FirstUseEver);
        if (ImGui::Begin("Flat")) {

            float w = ImGui::GetContentRegionAvailWidth();
            float p = ImGui::GetStyle().FramePadding.x;

            if (ImGui::CollapsingHeader("3D/2D", ImGuiTreeNodeFlags_DefaultOpen)){
                
                if (ImGui::Button("3D",  ImVec2((w - p) / 2.f, 0))){
                    morphing = false;
                    V = V_3D;
                    N = N_3D;
                    updateViz();
                    anim_time = 3.14/2.0;                
                }
                
                ImGui::SameLine(0, p);
                if (ImGui::Button("2D",  ImVec2((w - p) / 2.f, 0))){
                    morphing = false;
                    V = V_uv;
                    N = N_uv;
                    updateViz();
                    anim_time = 0.0;
                }
                
                ImGui::Checkbox("3D <-> 2D Morphing", &morphing);
                if (morphing) viewer.core().is_animating = true;
            }

            ImGui::Dummy(ImVec2(0.0f, 10.0f));
            if (ImGui::CollapsingHeader("Distortion", ImGuiTreeNodeFlags_DefaultOpen)){

                if (ImGui::Button("Symmetrize darts", ImVec2(-1, 0))){
                    for (auto dart: simple_darts){
                        Eigen::RowVector3d sym = dart.computeSymmetryAxis(V_uv_2d);
                        dart.snapSymmetric(V_uv_2d, sym);
                        V_uv.col(0) = V_uv_2d.col(0);
                        V_uv.col(1) = V_uv_2d.col(1);
                    }
                    V = V_uv;
                    updateViz();
                }

                if (ImGui::Button("Area", ImVec2((w - p) / 2.f, 0))){
                    dist_type = AREA_DIST;
                    viewer.data().show_texture = false;
                    updateViz();
                }
                ImGui::SameLine(0, p);
                if (ImGui::Button("Angle", ImVec2((w - p) / 2.f, 0))){
                    dist_type = ANGLE_DIST;
                    viewer.data().show_texture = false;
                    updateViz();
                }

                if (dist_type == AREA_DIST){
                    ImGui::Text("Area_UV / Area_3D");
                    ImGui::Text("Min-Max (Avg): %f - %f (%f)", area_dist.minCoeff(), area_dist.maxCoeff(), area_dist.mean());

                    ImGui::Text(std::to_string(min_area_ratio).substr(0,5).c_str());
                    std::string text = std::to_string(max_area_ratio).substr(0,5);
                    auto textWidth   = ImGui::CalcTextSize(text.c_str()).x;
                    ImGui::SameLine();
                    ImGui::SetCursorPosX((ImGui::GetWindowSize().x - 1) * 0.5f);
                    ImGui::Text("1");
                    ImGui::SameLine();
                    ImGui::SetCursorPosX(ImGui::GetWindowSize().x - 1.5*textWidth);
                    ImGui::Text(text.c_str());

                    drawColormap(igl::COLOR_MAP_TYPE_JET);
                                    
                }

                if (dist_type == ANGLE_DIST){
                    ImGui::Text("|(Area_UV - Area_3D)/ Area_3D|");
                    ImGui::Text("Min-Max (Avg): %f - %f (%f)", ang_dist.minCoeff(), ang_dist.maxCoeff(), ang_dist.mean());
                    drawColormap(igl::COLOR_MAP_TYPE_MAGMA);
                                    
                }

                std::vector<double> energies = flattening_energies;

                int arr_size = energies.size();
                double xs2[arr_size], ys2[arr_size];
                for (int i = 0; i < arr_size; ++i) {
                    xs2[i] = i * 1.0f;
                    ys2[i] = (float) energies[i];
                }


                #ifdef USE_IMPLOT
                ImPlot::CreateContext();

                //ImPlot::ShowDemoWindow();
                double max = *std::max_element(energies.begin(), energies.end());
                double min = *std::min_element(energies.begin(), energies.end());
                ImPlot::SetNextAxesLimits(0, arr_size, min, max);

                if (ImPlot::BeginPlot("Flattening energy", "Iteration", "Energy", ImVec2(-1, 200))) {
                    ImPlot::PlotLine("SCAF", xs2, ys2, arr_size);
                    ImPlot::EndPlot();
                }
                ImPlot::DestroyContext();
                #endif


            }

            ImGui::Dummy(ImVec2(0.0f, 10.0f));
            if (ImGui::CollapsingHeader("Display", ImGuiTreeNodeFlags_DefaultOpen)){
                make_checkbox("Show texture", viewer.data().show_texture);
                make_checkbox("Show lines", viewer.data().show_lines);
                make_checkbox("Show faces", viewer.data().show_faces);
                make_checkbox("Show v labels", viewer.data().show_vertex_labels);
                if (ImGui::Checkbox("Display X pattern", &display_x_pattern))
                    updateViz();
                ImGui::SliderFloat("Light factor", &viewer.core().lighting_factor, 0.0f, 5.0f, "%.3f");
                ImGui::SliderFloat("Time increment", &time_increment, 0.001f, 0.3f, "%.3f");
                if (ImGui::SliderFloat("UV factor", &uv_factor, 0.1f, 15.0f, "%.3f")){
                    viewer.data().set_uv(V_uv * uv_factor);
                    updateViz();
                }
            }
            
            //ImGui::ShowDemoWindow();

            ImGui::End();
        }
    };


    viewer.callback_post_draw = [&](igl::opengl::glfw::Viewer &v) {
        if (morphing)
        {
            anim_time += time_increment;
            double coeff = std::pow(std::sin(anim_time), 2);
            V = coeff * V_3D + (1-coeff) * V_uv;
            N = coeff * N_3D + (1-coeff) * N_uv;
            updateViz();
        }
        return true;
    };
    

    viewer.data().show_lines = 0u;
    viewer.data().line_width = 2;
    viewer.data().point_size = 6;
    viewer.core().orthographic = false;
    viewer.data().set_face_based(0u);
    viewer.data().show_texture = false;
    viewer.data().set_colors(Eigen::RowVector3d(1.0, 1.0, 1.0));
    viewer.core().lighting_factor = 0.7;
    //viewer.core().is_animating = true;

    updateViz();

    viewer.core().set_rotation_type(igl::opengl::ViewerCore::ROTATION_TYPE_TRACKBALL);
    viewer.core().background_color = Eigen::Vector4f(202.0/255.0, 190.0/255.0, 232.0/255.0, 1.0);
    viewer.launch();
}