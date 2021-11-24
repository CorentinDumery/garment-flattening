

#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/readOBJ.h>

#include <igl/triangle_triangle_adjacency.h>
#include <igl/boundary_facets.h>


int main(int argc, char *argv[]){

    Eigen::MatrixXd V_3d, V;
    Eigen::MatrixXi F;

    igl::readOBJ("../data/flat_dress.obj", V_3d, F);

    V_3d = V_3d.rowwise() - V_3d.colwise().minCoeff();

    int n_fibers = 50;

    V_3d *= static_cast<double>(n_fibers) / (V_3d.maxCoeff() - V_3d.minCoeff());

    std::cout << "Mesh range:" << std::endl;
    std::cout << "X: " << V_3d.col(0).minCoeff() << " -> " << V_3d.col(0).maxCoeff()  << std::endl;
    std::cout << "Y: " << V_3d.col(1).minCoeff() << " -> " << V_3d.col(1).maxCoeff()  << std::endl;
    std::cout << "Z: " << V_3d.col(2).minCoeff() << " -> " << V_3d.col(2).maxCoeff()  << std::endl;

    V.resize(V_3d.rows(), 2);
    V.col(0) = V_3d.col(0);
    V.col(1) = V_3d.col(1);


    Eigen::MatrixXi TT, Eb;
    igl::triangle_triangle_adjacency(F, TT);

    Eigen::VectorXi J, K;
    igl::boundary_facets(F, Eb, J, K);

    Eigen::MatrixXd edge_begs(Eb.rows(), 2), edge_ends(Eb.rows(), 2);

    for (int i=0; i<Eb.rows(); i++){
        edge_begs.row(i) = V.row(Eb(i,0));
        edge_ends.row(i) = V.row(Eb(i,1));
    }

    auto interpolateToVal = [&V](int target, const Eigen::RowVector2d& v0, const Eigen::RowVector2d& v1, int axis){
        double alpha0 = (static_cast<double>(target) - v0(axis)) / (v1(axis) - v0(axis)); 
        return (1-alpha0) * v0 + (alpha0) * v1;
    };

    std::vector<Eigen::MatrixXd> fiber_begs_list;
    std::vector<Eigen::MatrixXd> fiber_ends_list;

    for (int axis=0; axis<2; axis++){ // For U and V

        int min_ax = std::floor(V.col(axis).minCoeff()) + 1;
        int max_ax = std::floor(V.col(axis).maxCoeff()) - 1;

        std::map<int, std::vector<int>> fiber_axis_intersec;

        for (int i=min_ax; i<=max_ax; i++){
            fiber_axis_intersec[i] = {};
        }

        for (int i=0; i<Eb.rows(); i++){
            Eigen::RowVector2d v0 = V.row(Eb(i,0));
            Eigen::RowVector2d v1 = V.row(Eb(i,1));
            int e_min_ax = std::floor(std::min(v0(axis), v1(axis))) + 1;
            int e_max_ax = std::floor(std::max(v0(axis), v1(axis))) + 0;
            for (int j = e_min_ax; j<= e_max_ax; j++){
                fiber_axis_intersec[j].push_back(i);
            }
        }

        int e_count = 0;

        for (auto const& x : fiber_axis_intersec){
            std::vector<int> vals = x.second; 
            for (int j: vals) {
                e_count ++;
            }
        }

        if (e_count % 2 != 0) std::cout << "ERROR: odd number of edge intersections?" << std::endl;

        e_count /= 2;

        // Sort following other axis
        for (auto & x : fiber_axis_intersec){
            std::vector<int> vals = x.second; 
            for (int j=0; j<vals.size(); j ++) {
                for (int k=j+1; k<vals.size(); k ++) {
                    double jv = ((V.row(Eb(vals[j], 0)) + V.row(Eb(vals[j], 1)))/2.0)((axis + 1) % 2); // not worth calling interpolateToVal here?
                    double kv = ((V.row(Eb(vals[k], 0)) + V.row(Eb(vals[k], 1)))/2.0)((axis + 1) % 2);
                    if (jv > kv){
                        int temp = vals[k];
                        vals[k] = vals[j];
                        vals[j] = temp;
                    }
                }
            }
            x.second = vals;
        }

        // Visualize fiber edges
        Eigen::MatrixXd fiber_begs(e_count, 2), fiber_ends(e_count, 2);
        int curr_fib_id = 0;
        for (auto const& x : fiber_axis_intersec){
            std::vector<int> vals = x.second; 
            for (int j=0; j<vals.size(); j += 2) {
                fiber_begs.row(curr_fib_id) = interpolateToVal(x.first, V.row(Eb(vals[j], 0)), V.row(Eb(vals[j], 1)), axis);
                fiber_ends.row(curr_fib_id) = interpolateToVal(x.first, V.row(Eb(vals[j+1], 0)), V.row(Eb(vals[j+1], 1)), axis);
                curr_fib_id ++;
            }
        }

        fiber_begs_list.push_back(fiber_begs);
        fiber_ends_list.push_back(fiber_ends);
        
    }





    // --- VISUALIZATION ---

   
    igl::opengl::glfw::Viewer viewer;

    igl::opengl::glfw::imgui::ImGuiMenu menu;
    menu.callback_draw_viewer_window = []() {};
    viewer.plugins.push_back(&menu);

    viewer.data().add_edges(edge_begs, edge_ends, Eigen::RowVector3d(1.0, 0.5, 0.6));
    viewer.data().add_edges(fiber_begs_list[0], fiber_ends_list[0], Eigen::RowVector3d(0.5, 1.0, 0.6));
    viewer.data().add_edges(fiber_begs_list[1], fiber_ends_list[1], Eigen::RowVector3d(0.5, 0.5, 1.0));
    viewer.data().set_mesh(V, F);

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
        if (ImGui::Begin("Net")) {
            
            make_checkbox("Show mesh", viewer.data().show_lines);
            ImGui::End();
        }
    
    };

    //updateViz();


    viewer.data().line_width = 5;

    viewer.data().show_lines = false;

    viewer.data().point_size = 10;
    viewer.core().orthographic = true;

    viewer.core().set_rotation_type(igl::opengl::ViewerCore::ROTATION_TYPE_NO_ROTATION);
    viewer.core().background_color = Eigen::Vector4f(202.0/255.0, 190.0/255.0, 232.0/255.0, 1.0);
    viewer.launch();
}