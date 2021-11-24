#include "net_param.h"

#include <igl/triangle_triangle_adjacency.h>
#include <igl/boundary_facets.h>

#include <iostream>
#include <vector>
#include <map>

#define OFFSCREEN_RENDERING_COORDS false

// TODO put these functions somewhere
// glfw: when window size changed this callback function executes // glfwSetFramebufferSizeCallback
void framebuffer_size_callback(GLFWwindow* window, int width, int height) {
    // glfw: the viewport matches the new window dimensions
    glViewport(0, 0, width, height);
}
 
// glfw: process keys
void processInput(GLFWwindow* window) {
    if (glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS)
        glfwSetWindowShouldClose(window, true);
}


NetParam::NetParam(const Eigen::MatrixXi& F,
                   const Eigen::MatrixXf& V_3d,
                   const Eigen::MatrixXf& V_2d)
                   : F_(F), V_3d_(V_3d), V_2d_(V_2d){
    n_tris_ = F_.rows();

    igl::triangle_triangle_adjacency(F_, TT_);

    Eigen::VectorXi J, K;
    igl::boundary_facets(F_, Eb_, J, K);
}

void NetParam::initializeRendering(){


    std::vector<float> V_2d_vec;
    // TODO REPLACE THIS
    for (int i=0; i<V_2d_.rows(); i++){
        for (int j=0; j<V_2d_.cols(); j++){
            V_2d_vec.push_back(V_2d_(i,j));
        }    
    }

    std::vector<float> V_3d_vec;
    // TODO REPLACE THIS
    for (int i=0; i<V_3d_.rows(); i++){
        for (int j=0; j<V_3d_.cols(); j++){
            V_3d_vec.push_back(V_3d_(i,j));
        }    
    }

    std::vector<unsigned int> indices_vec; 
    // TODO REPLACE THIS
    for (int i=0; i<F_.rows(); i++){
        for (int j=0; j<F_.cols(); j++){
            indices_vec.push_back(static_cast<unsigned int>(F_(i,j)));
        }    
    }


    //*
    float vertices_2d[V_2d_vec.size()];
    std::copy(V_2d_vec.begin(), V_2d_vec.end(), vertices_2d);

    float vertices_3d[V_3d_vec.size()];
    std::copy(V_3d_vec.begin(), V_3d_vec.end(), vertices_3d);

    unsigned int indices[indices_vec.size()];
    std::copy(indices_vec.begin(), indices_vec.end(), indices);

    //int n_tris = indices_vec.size() / 3;
    //*/



    // glfw: initialize
    glfwInit();
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 5);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

    // OFFSCREEN RENDERING https://github.com/glfw/glfw/blob/master/examples/offscreen.c
    if (OFFSCREEN_RENDERING_COORDS){
        glfwWindowHint(GLFW_VISIBLE, GLFW_FALSE);
        glfwInitHint(GLFW_COCOA_MENUBAR, GLFW_FALSE); // no menubar for mac users
    }

    // glfw: create window
    window = glfwCreateWindow(RENDER_WIDTH, RENDER_HEIGHT, "Window Title", NULL, NULL);
    if (window == NULL) {
        std::cout << "Failed to create GLFW window" << "\n";
        glfwTerminate();
    }
    glfwMakeContextCurrent(window);
    glfwSetFramebufferSizeCallback(window, framebuffer_size_callback);
 
    // glad: load all OpenGL function pointers
    if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress)) {
        std::cout << "Failed to initialize GLAD" << "\n";
    }
 
    // build and compile shader program
    int success;
    char errorInfo[512] = "";

    const char* vertexShaderSource = 
    "#version 450 core\n"
    "layout (location = 0) in vec2 aPos;\n"
    "layout (location = 1) in vec3 pos3d;\n"
    "out vec3 interpColor;\n"
    "void main()\n"
    "{\n"
    "   gl_Position = vec4(aPos.x, aPos.y, 0.0, 1.0);\n"
    "   interpColor= pos3d;\n"
    "}\n\0";
    
    const char* fragmentShaderSource = 
    "#version 450 core\n"
    "in vec3 interpColor;\n"
    "out vec4 FragmentColor;\n"
    "void main()\n"
    "{\n"
    "   FragmentColor = vec4(interpColor, 1.0);\n"
    "}\n\0";
 
    // vertex shader
    unsigned int vertexShader = glCreateShader(GL_VERTEX_SHADER);
    glShaderSource(vertexShader, 1, &vertexShaderSource, NULL);
    glCompileShader(vertexShader);
 
    glGetShaderiv(vertexShader, GL_COMPILE_STATUS, &success);
    if (!success) {
        glGetShaderInfoLog(vertexShader, 512, NULL, errorInfo);
        std::cout << "ERROR::VERTEX::SHADER::COMPILATION_FAILED\n" << errorInfo << "\n";
    }
 
    // fragment shader
    unsigned int fragmentShader = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(fragmentShader, 1, &fragmentShaderSource, NULL);
    glCompileShader(fragmentShader);
 
    glGetShaderiv(fragmentShader, GL_COMPILE_STATUS, &success);
    if (!success) {
        glGetShaderInfoLog(fragmentShader, 512, NULL, errorInfo);
        std::cout << "ERROR::FRAGMENT::SHADER::COMPILATION_FAILED\n" << errorInfo << "\n";
    }
 
    // link shaders
    shaderProgram = glCreateProgram();
    glAttachShader(shaderProgram, vertexShader);
    glAttachShader(shaderProgram, fragmentShader);
    glLinkProgram(shaderProgram);
 
    glGetProgramiv(shaderProgram, GL_LINK_STATUS, &success);
    if (!success) {
        glGetProgramInfoLog(shaderProgram, 512, NULL, errorInfo);
        std::cout << "ERROR::PROGRAM::LINKING_FAILED\n" << errorInfo << "\n";
    }
 
    // delete shaders after linking
    glDeleteShader(vertexShader);
    glDeleteShader(fragmentShader);
 
    // set vertex data 
    /*
    std::vector<float> vec2d = {
       0.5f,  0.5f,  // top right // 0
       0.5f, -0.5f,  // bottom right // 1
      -0.5f, -0.5f,  // bottom left // 2
      -0.5f,  0.5f,   // top left  // 3
        0.0,  1.0 // mid top
    };

    float vertices_2d[vec2d.size()];
    std::copy(vec2d.begin(), vec2d.end(), vertices_2d);

    float vertices_3d[] = {
       0.0f,  0.0f, 0.0,  // top right // 0
       0.0f, -0.0f, 0.25,  // bottom right // 1
      -0.0f, -0.0f, 0.5,  // bottom left // 2
      -0.0f,  0.0f,  0.75,  // top left  // 3
      -0.0f,  1.0f,  0.75,  // mid top
    };
 
    // index buffer // Element Buffer Objects (EBO)
    unsigned int indices[] = {  
        0, 3, 1,   // first triangle
        1, 2, 3,    // second triangle
        0, 3, 4
    };
    int n_tris = 3;
    //*/
 
    // set vertex buffer object anb vertex array object and element buffer objects 
    
    glGenVertexArrays(1, &VAO);
    glGenBuffers(1, &VBO);
    glGenBuffers(1, &EBO);
     
    // bind vertex array object
    glBindVertexArray(VAO);
 
    // bind vertex buffer object
    glBindBuffer(GL_ARRAY_BUFFER, VBO);
    glBufferData(GL_ARRAY_BUFFER, sizeof(vertices_2d), vertices_2d, GL_STATIC_DRAW);
 
    

    // bind element buffer objects
    // EBO is stored in the VAO
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(indices), indices, GL_STATIC_DRAW);
 
    // registered VBO as the vertex attributes
    glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 2 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);
 
    // --- Adding data to vertices, see https://learnopengl.com/code_viewer_gh.php?code=src/2.lighting/1.colors/colors.cpp

    //*
    unsigned int VBO_3D;
    glGenBuffers(1, &VBO_3D);

    glBindBuffer(GL_ARRAY_BUFFER, VBO_3D);
    glBufferData(GL_ARRAY_BUFFER, sizeof(vertices_3d), vertices_3d, GL_STATIC_DRAW);

    GLint pos3d_attrib = glGetAttribLocation(shaderProgram, "pos3d");
    glEnableVertexAttribArray(pos3d_attrib);
    glVertexAttribPointer(pos3d_attrib, 3, GL_FLOAT, GL_FALSE, 0, 0);
    //*/

    // --- 

    // unbind the VAO
    glBindBuffer(GL_ARRAY_BUFFER, 0);
    glBindVertexArray(0);
}

void NetParam::freeRenderingBuffers(){
 // de-allocate all resources
    glDeleteVertexArrays(1, &VAO);
    glDeleteBuffers(1, &VBO);
    glDeleteBuffers(1, &EBO);
    glDeleteProgram(shaderProgram);
 
    // glfw: terminate and clear all previously GLFW allocated resources
    glfwTerminate();

}


void NetParam::render(){
    while (!glfwWindowShouldClose(window)) {
        if (! OFFSCREEN_RENDERING_COORDS) 
            processInput(window);
 
        // render
        glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT);
 
        // draw triangle
        glUseProgram(shaderProgram);
        glBindVertexArray(VAO); 
        glDrawElements(GL_TRIANGLES, 3 * n_tris_, GL_UNSIGNED_INT, 0);
 
        //time_post_render = steady_clock::now();

        int mx = 30;
        float pixels[ 1 * RENDER_HEIGHT * 3 ];
        glReadPixels(mx, 0, 1, RENDER_HEIGHT, GL_RGB, GL_FLOAT, pixels);

        double dist = 0;
        for (int i=0; i<RENDER_HEIGHT-1; i+=3){
            if (pixels[3*i] == 0 && pixels[3*i+1] == 0 && pixels[3*i+2] == 0) continue; // note: we might miss a pixel at the origin in 3D
            if (pixels[3*(i+1)] == 0 && pixels[3*(i+1)+1] == 0 && pixels[3*(i+1)+2] == 0) continue;
            //std::cout << "(" << pixels[3*i] << " " << pixels[3*i+1] << " " << pixels[3*i+2] << ")" << std::endl;
            dist += std::sqrt(std::pow(pixels[3*(i+1)  ] - pixels[3*i  ], 2) 
                            + std::pow(pixels[3*(i+1)+1] - pixels[3*i+1], 2) 
                            + std::pow(pixels[3*(i+1)+2] - pixels[3*i+2], 2));
        }

        std::cout << "dist = " << dist << std::endl; 

        //time_post_read = steady_clock::now();

        if (OFFSCREEN_RENDERING_COORDS) break;

        // glfw: swap buffers
        glfwSwapBuffers(window);
 
        // glfw: poll IO events (keys & mouse) 
        // (including X close window button)
        glfwPollEvents();
 
    }
}

void NetParam::vizBoundaryEdges(Eigen::MatrixXd& edge_begs, Eigen::MatrixXd& edge_ends){
    edge_begs.resize(Eb_.rows(), 2);
    edge_ends.resize(Eb_.rows(), 2);

    for (int i=0; i<Eb_.rows(); i++){
        edge_begs.row(i) = V_2d_.row(Eb_(i,0)).cast<double>();
        edge_ends.row(i) = V_2d_.row(Eb_(i,1)).cast<double>();
    }
}

void NetParam::computeFibers(){
    
    auto interpolateToVal = [](int target, const Eigen::RowVector2f& v0, const Eigen::RowVector2f& v1, int axis){
        double alpha0 = (static_cast<double>(target) - v0(axis)) / (v1(axis) - v0(axis)); 
        return (1-alpha0) * v0 + (alpha0) * v1;
    };


    for (int axis=0; axis<2; axis++){ // For U and V

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
}