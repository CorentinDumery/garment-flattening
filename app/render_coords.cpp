// see https://openglbook.com/chapter-3-index-buffer-objects-and-primitive-types.html

#include <glad/glad.h>
#include <GLFW/glfw3.h>

#include <igl/readOBJ.h>

#include <iostream>
#include <chrono>

#include "net_param.h"

#define OFFSCREEN_RENDERING_COORDS false // TODO remove, also in net_param.cpp

using std::chrono::steady_clock;
using std::chrono::duration_cast;
using std::chrono::microseconds; 
 
int main(){

    auto time_init = steady_clock::now();

    Eigen::MatrixXd V_2d_d, V_3d_d;
    Eigen::MatrixXi F;
    igl::readOBJ("../data/dress_front_cut.obj", V_3d_d, F);
    igl::readOBJ("../data/flat_dress.obj", V_2d_d, F); // note: could also be encoded as UV in a single OBJ

    int n_vertices = V_2d_d.rows();
    Eigen::MatrixXf V_2d(n_vertices, 2), V_3d(n_vertices, 3);
    V_2d.col(0) = V_2d_d.col(0).cast<float>();
    V_2d.col(1) = V_2d_d.col(1).cast<float>();
    V_3d = V_3d_d.cast<float>();

    if (V_2d.rows() != V_3d.rows() || V_2d.rows() < 1){
        std::cout << "ERROR, mesh mismatch:" << std::endl;
        std::cout <<  V_2d.rows() << " vs " << V_3d.rows() << std::endl;
        return 0;
    }

    NetParam net_param(F, V_3d, V_2d);


    net_param.initializeRendering();
 
    auto time_pre_render = steady_clock::now();
    auto time_post_render = steady_clock::now(); // note: value will be overwritten, but couldn't initialize otherwise
    auto time_post_read = steady_clock::now();
    // render loop

    net_param.render();

    std::cout << "Rendered" << std::endl;
    
    Eigen::RowVector2d start, end;
    start << 0.4, 0.4;
    end << 0.5, 0.4;
    double length = net_param.measureFiber(start, end);

    std::cout << "Length: " << length << std::endl;


    net_param.freeRenderingBuffers();
 


    auto time_final = steady_clock::now();

    if (OFFSCREEN_RENDERING_COORDS){
        std::cout << "Total time (μs): "<< duration_cast<microseconds>(time_final - time_init).count() << std::endl;
        std::cout << "\tInit  : "<< duration_cast<microseconds>(time_pre_render - time_init).count() << std::endl;
        std::cout << "\tRender: "<< duration_cast<microseconds>(time_post_render - time_pre_render).count() << std::endl;
        std::cout << "\tRead  : "<< duration_cast<microseconds>(time_post_read - time_post_render).count() << std::endl;
    }
 
    return 0;
}