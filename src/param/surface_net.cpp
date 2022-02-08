#include "surface_net.h"

#include <igl/triangle_triangle_adjacency.h>
#include <igl/boundary_facets.h>
#include <igl/barycentric_coordinates.h>
#include <iostream>
#include <vector>
#include <map>

SurfaceNet::SurfaceNet(const Eigen::MatrixXi& F,
                   const Eigen::MatrixXf& V_3d,
                   const Eigen::MatrixXf& V_2d)
                   : F_(F), V_3d_(V_3d), V_2d_(V_2d){
    n_tris_ = F_.rows();

    igl::triangle_triangle_adjacency(F_, TT_);

    Eigen::VectorXi J, K;
    igl::boundary_facets(F_, Eb_, J, K);

    if (V_2d_.cols() > 2){
        Eigen::MatrixXf temp = V_2d_;
        V_2d_.resize(V_2d_.rows(), 2);
        V_2d_.col(0) = temp.col(0);
        V_2d_.col(1) = temp.col(1);
    }
}
    
void SurfaceNet::vizBoundaryEdges(Eigen::MatrixXd& edge_begs, Eigen::MatrixXd& edge_ends) const {
    edge_begs.resize(Eb_.rows(), 2);
    edge_ends.resize(Eb_.rows(), 2);

    for (int i=0; i<Eb_.rows(); i++){
        edge_begs.row(i) = V_2d_.row(Eb_(i,0)).cast<double>();
        edge_ends.row(i) = V_2d_.row(Eb_(i,1)).cast<double>();
    }

    //f romRenderToInitCoords(edge_begs);
    //f romRenderToInitCoords(edge_ends);
}

void SurfaceNet::computeNet(){
    
    auto interpolateToVal = [](int target, const Eigen::RowVector2f& v0, const Eigen::RowVector2f& v1, int axis){
        double alpha0 = (static_cast<double>(target) - v0(axis)) / (v1(axis) - v0(axis)); 
        return (1-alpha0) * v0 + (alpha0) * v1;
    };

    // Eigen::MatrixXf V_2d_ = V_2d_;
    // f romRenderToInitCoords(V_2d_);


    for (int axis=0; axis<2; axis++){ // For U and V
        #ifdef NET_PARAM_DEBUG
        std::cout << "Compute fibers, axis " << axis << std::endl;
        #endif

        int min_ax = std::floor(V_2d_.col(axis).minCoeff()) + 1;
        int max_ax = std::floor(V_2d_.col(axis).maxCoeff()) - 1;

        std::map<int, std::vector<int>> fiber_axis_intersec;

        for (int i=min_ax; i<=max_ax; i++){
            fiber_axis_intersec[i] = {};
        }

        for (int i=0; i<Eb_.rows(); i++){
            Eigen::RowVector2f v0 = V_2d_.row(Eb_(i,0));
            Eigen::RowVector2f v1 = V_2d_.row(Eb_(i,1));
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

        #ifdef NET_PARAM_DEBUG
        std::cout << "Compute fibers, e_count " << e_count << std::endl;
        if (e_count <= 0) std::cout << "ERROR: not enough edge intersections?" << std::endl;
        #endif

        // Sort following other axis
        for (auto & x : fiber_axis_intersec){
            std::vector<int> vals = x.second; 
            for (int j=0; j<vals.size(); j ++) {
                for (int k=j+1; k<vals.size(); k ++) {
                    double jv = ((V_2d_.row(Eb_(vals[j], 0)) + V_2d_.row(Eb_(vals[j], 1)))/2.0)((axis + 1) % 2); // not worth calling interpolateToVal here?
                    double kv = ((V_2d_.row(Eb_(vals[k], 0)) + V_2d_.row(Eb_(vals[k], 1)))/2.0)((axis + 1) % 2);
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
                fiber_begs.row(curr_fib_id) = interpolateToVal(x.first, V_2d_.row(Eb_(vals[j], 0)), V_2d_.row(Eb_(vals[j], 1)), axis).cast<double>();
                fiber_ends.row(curr_fib_id) = interpolateToVal(x.first, V_2d_.row(Eb_(vals[j+1], 0)), V_2d_.row(Eb_(vals[j+1], 1)), axis).cast<double>();
                curr_fib_id ++;
            }
        }

        fiber_begs_list.push_back(fiber_begs);
        fiber_ends_list.push_back(fiber_ends);
        
    }

    #ifdef NET_PARAM_DEBUG
    std::cout << "Compute fibers done" << std::endl;
    #endif
}

std::vector<std::vector<int>> SurfaceNet::nearestFibers(const Eigen::RowVector2f& vertex) const {

    // Vertex in init coordinates

    // TODO OPTIMIZE !!
    // with integer lines and/or proper data structure, shouldn't even need to search

    // TODO consider different lines with = constant

    std::vector<std::vector<int>> selected_fibers;
    for (int axis=0; axis<2; axis++){
        Eigen::VectorXd fibers = fiber_begs_list[axis].col(axis);
        double t = vertex(axis);

        int fib1 = -1; // left fiber 
        int fib2 = 0; // right
        
        int pos = 0;
        //std::cout << "Fibrows " << fibers.rows() << std::endl;
        //std::cout << t << " <? " << fibers(pos + 1) << std::endl;
        while (pos < fibers.rows() && t >= fibers(pos)){
            pos ++;
            fib1 ++;
            fib2 ++;
        }

        if (pos >= fibers.rows()) fib2 = -1;

        // there can be several fibers with that same value
        // pick the one on which you can project
        
        // the following assumes fibers are ordered by main axis, and fib1 fib2 are the last (resp first)
        // before (resp after) vertex's axis value. But since there can be many with that value,
        // try to fib1 -- or fib2 ++ until you find one you can project on
        // if you can't project on such a line with = axis value to fib1, then return -1 
        auto canProjOnFiber = [](float v_other_axis, double fib_beg, double fib_end){
            return v_other_axis > fib_beg && v_other_axis < fib_end;
        };

        float init_dist1 = vertex(axis) - fiber_begs_list[axis](fib1, axis);
        while (true && fib1 >= 0){
            if (!canProjOnFiber(vertex((axis+1) % 2), fiber_begs_list[axis](fib1, (axis+1) % 2), fiber_ends_list[axis](fib1, (axis+1) % 2))){
                fib1 --;
            }
            else break;
            if (vertex(axis) - fiber_begs_list[axis](fib1, axis) > init_dist1 * 1.1){
                fib1 = -1;
                break;
            }
        }

        float init_dist2 = fiber_begs_list[axis](fib2, axis) - vertex(axis);
        while (true && fib2 < fiber_begs_list[axis].rows() && fib2 != -1){
            if (!canProjOnFiber(vertex((axis+1) % 2), fiber_begs_list[axis](fib2, (axis+1) % 2), fiber_ends_list[axis](fib2, (axis+1) % 2))){
                fib2 ++;
            }
            else break;
            if (fiber_begs_list[axis](fib2, axis) - vertex(axis) > init_dist2 * 1.1 || fib2 == fiber_begs_list[axis].rows()){
                fib2 = -1;
                break;
            }
        }

        selected_fibers.push_back({fib1, fib2});
    }


    std::cout << "SELECTED FIBERS:" << std::endl;
    std::cout << "U:" << selected_fibers[0][0] << " " <<  selected_fibers[0][1] << std::endl;
    std::cout << "V:" << selected_fibers[1][0] << " " <<  selected_fibers[1][1] << std::endl;

    return selected_fibers;
}
