/**
 * @author Corentin Dumery
 * @brief Finds the best reflection between two similar lines in a least squares sense,
 * and forces reflectability
 * @date 2022-02-04
 * 
 */
#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <cmath>

#include <param/param_utils.h>

void forceReflectability(Eigen::MatrixXd& line1,
                         Eigen::MatrixXd& line2){
    Eigen::MatrixXd R_est;
    Eigen::VectorXd T_est;
    procrustes(line1, line2, R_est, T_est);

    Eigen::MatrixXd line1t = line1.transpose();
    Eigen::MatrixXd line2t = line2.transpose();
    Eigen::MatrixXd line3;

    /* FROM 2 TO 1 */
    line3 = line2t.colwise() - T_est;
    line3 = (R_est.transpose() * line3);

    Eigen::MatrixXd diff = line1t - line3;
    line1t = line1t - diff / 2.0;
    line3 = line1t;

    /* FROM 1 TO 2 */
    line3 = (R_est * line3);
    line3 = line3.colwise() + T_est;

    line2t = line3;
    
    line1 = line1t.transpose();
    line2 = line2t.transpose();
}

int main(int argc, char *argv[]){

    Eigen::MatrixXd line1(17, 2), line2(17, 2);

    line1 << -0.9393371343612671, -0.30847856402397156,
            1.0, 0.17045824229717255,
            0.5, 0.06720391660928726,
            0.00527503015473485, -0.09328988939523697,
            -0.481537401676178, -0.08783920109272003,
            -0.7104372978210449, -0.062260814011096954,
            0.75, 0.09282717108726501,
            -0.21834981441497803, -0.2868318259716034,
            0.27901268005371094, -0.026927676051855087,
            -0.5854372978210449, -0.062260814011096954,
            0.875, 0.10548053681850433,
            -0.09334981441497803, -0.2598654329776764,
            0.393462598323822, 0.01758262887597084,
            -0.8143371343612671, -0.2282152622938156,
            0.6329125165939331, 0.028557751327753067,
            -0.356537401676178, -0.1447737216949463,
            0.125, -0.02570086158812046;

    line2 << -1.0, 0.30847856402397156,
            1.0, -0.23112109303474426,
            0.5, -0.1384168267250061,
            0.0, 0.06691473722457886,
            -0.5, 0.10893931984901428,
            -0.75, 0.03588566184043884,
            0.75, -0.13766492903232574,
            -0.25, 0.2551816403865814,
            0.25, 0.0031900405883789062,
            -0.625, 0.03588566184043884,
            0.875, -0.15031829476356506,
            -0.125, 0.5282152622938156,
            0.375, -0.11253317445516586,
            -0.875, 0.2282152622938156,
            0.625, -0.049657873809337616,
            -0.375, 0.16587384045124054,
            0.125, 0.02570086158812046;

    double theta = 0.75;
    Eigen::MatrixXd R(2,2);
    R << std::cos(theta), - std::sin(theta),
         std::sin(theta), std::cos(theta);

    Eigen::RowVector2d T(2.1, 1.3);

    line2 = (R * line2.transpose()).transpose();
    line2 = line2.rowwise() + T;
    
    Eigen::MatrixXd line3, line3t;
    auto compLine3 = [&](){
        Eigen::MatrixXd R_est;
        Eigen::VectorXd T_est;
        procrustes(line1, line2, R_est, T_est);
        
        Eigen::MatrixXd line3t = line1.transpose();
        line3t = (R_est * line3t);
        line3t = line3t.colwise() + T_est;
        line3 = line3t.transpose();
    };

    compLine3();

    // --- VISUALIZATION --- //

    Eigen::MatrixXi E(16, 2);
    E << 9, 4, 
        10, 1, 
        11, 3, 
        12, 2, 
        13, 5, 
        14, 6, 
        15, 7, 
        16, 8, 
         5, 9, 
        6, 10, 
        7, 11, 
        8, 12, 
        0, 13, 
        2, 14, 
        4, 15, 
        3, 16;

    auto lineFromPoints = [&](Eigen::MatrixXd points){
        Eigen::MatrixXd begs(E.rows(), 2), ends(E.rows(), 2);
        for (int i=0; i<E.rows(); i++){
            begs.row(i) = points.row(E(i,0));
            ends.row(i) = points.row(E(i,1));
        }
        return std::make_pair(begs, ends);
    };

    igl::opengl::glfw::Viewer viewer;

    igl::opengl::glfw::imgui::ImGuiMenu menu;
    menu.callback_draw_viewer_window = []() {};
    viewer.plugins.push_back(&menu);

    auto updateViz = [&](){
        viewer.data().clear_edges();
        viewer.data().clear_points();
        compLine3();
        
        viewer.data().add_points(line1, Eigen::RowVector3d(1.0, 0.0, 0.0));
        auto edges1 = lineFromPoints(line1);
        viewer.data().add_edges(edges1.first, edges1.second, Eigen::RowVector3d(1.0, 0.0, 0.0));
        
        viewer.data().add_points(line2, Eigen::RowVector3d(0.0, 0.0, 1.0));
        auto edges2 = lineFromPoints(line2);
        viewer.data().add_edges(edges2.first, edges2.second, Eigen::RowVector3d(0.0, 0.0, 1.0));

        viewer.data().add_points(line3, Eigen::RowVector3d(1.0, 0.0, 1.0));
        auto edges3 = lineFromPoints(line3);
        viewer.data().add_edges(edges3.first, edges3.second, Eigen::RowVector3d(1.0, 0.0, 1.0));
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
        if (ImGui::Begin("Procrustes")) {
            ImGui::Text("Distance Blue <-> Pink: %f", (line2 - line3).norm());
            if (ImGui::Button("Align lines", ImVec2(-1, 0))){
                forceReflectability(line1, line2);
                updateViz();
            }
            ImGui::End();
        }
    };

    updateViz();

    viewer.data().line_width = 5;
    viewer.data().point_size = 10;
    viewer.core().orthographic = true;
    viewer.core().set_rotation_type(igl::opengl::ViewerCore::ROTATION_TYPE_TRACKBALL);
    viewer.core().background_color = Eigen::Vector4f(202.0/255.0, 190.0/255.0, 232.0/255.0, 1.0);
    viewer.launch();
}