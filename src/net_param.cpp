#include "net_param.h"

#include <igl/triangle_triangle_adjacency.h>
#include <igl/boundary_facets.h>
#include <igl/barycentric_coordinates.h>
#include <igl/jet.h>

#include <igl/opengl/glfw/background_window.h>
#include <igl/opengl/init_render_to_texture.h>

#include <iostream>
#include <vector>
#include <map>

//#define NET_PARAM_DEBUG
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

void NetParam::fromInitToRenderCoords(Eigen::MatrixXf& points) const {
    points.col(0) = points.col(0).array() - off_u_;
    points.col(1) = points.col(1).array() - off_v_;
    points.col(0) /= scale_u_; // (actually make it [-0.95,-0.95] so that UV can change a bit without needing to re-scale)
    points.col(1) /= scale_v_;
    points = points.array() - 0.95;
}

void NetParam::fromRenderToInitCoords(Eigen::MatrixXf& points) const{
    points = points.rowwise() + Eigen::RowVector2f(0.95, 0.95);

    points.col(0) *= scale_u_;
    points.col(1) *= scale_v_;

    points = points.rowwise() + Eigen::RowVector2f(off_u_, off_v_);
}

void NetParam::fromRenderToInitCoords(Eigen::MatrixXd& points) const{
    Eigen::MatrixXf bis = points.cast<float>();
    fromRenderToInitCoords(bis); // TODO template function instead
    points = bis.cast<double>();
}


NetParam::NetParam(const Eigen::MatrixXi& F,
                   const Eigen::MatrixXf& V_3d,
                   const Eigen::MatrixXf& V_2d)
                   : F_(F), V_3d_(V_3d), V_2d_(V_2d){
    n_tris_ = F_.rows();

    igl::triangle_triangle_adjacency(F_, TT_);

    Eigen::VectorXi J, K;
    igl::boundary_facets(F_, Eb_, J, K);

    if (V_2d_.cols() > 2){
        std::cout << "WARNING: received 2D mesh in a 3D matrix (at least), discarding other columns" << std::endl;
        Eigen::MatrixXf temp = V_2d_;
        V_2d_.resize(V_2d_.rows(), 2);
        V_2d_.col(0) = temp.col(0);
        V_2d_.col(1) = temp.col(1);
    }

    //std::cout << Eb_ << std::endl;

    // Scaling mesh:    
    off_u_ = V_2d_.col(0).minCoeff();
    off_v_ = V_2d_.col(1).minCoeff();
    //V_2d_.col(0) = V_2d_.col(0).array() - off_u_;
    //V_2d_.col(1) = V_2d_.col(1).array() - off_v_;
    //V_2d_.col(0) /= V_2d_.col(0).maxCoeff()/2.0; // directly scale 2D mesh to [-1,1]x[-1,1] for rendering
    //V_2d_.col(1) /= V_2d_.col(1).maxCoeff()/2.0;
    // V_2d_ = V_2d_.array() - 1.0;

    scale_u_ = (V_2d_.col(0).array() - off_u_).maxCoeff()/1.9;
    scale_v_ = (V_2d_.col(1).array() - off_v_).maxCoeff()/1.9; 
    //V_2d_.col(0) /= scale_u_; // (actually make it [-0.95,-0.95] so that UV can change a bit without needing to re-scale)
    //V_2d_.col(1) /= scale_v_;
    //V_2d_ = V_2d_.array() - 0.95;

    //V_3d_ = V_3d_.array() - V_3d_.minCoeff();
    //V_3d_ /= V_3d_.maxCoeff();

    RENDER_HEIGHT = static_cast<int>(static_cast<float>(RENDER_WIDTH) * scale_v_ / scale_u_);
}

void NetParam::prepareShaders(){
     
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
    "layout(location = 0) out vec3 FragmentColor;\n"
    "void main()\n"
    "{\n"
    "   FragmentColor = interpColor;\n"
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
}

void NetParam::prepareBuffers(){


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

void NetParam::initializeRendering(){

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

    prepareShaders();
 
    prepareBuffers();

    std::cout << "Window initialization done" << std::endl;
}

void NetParam::freeRenderingBuffers(){
    // de-allocate all resources
    glBindVertexArray(0);
    glDeleteVertexArrays(1, &VAO);
    glDeleteBuffers(1, &VBO);
    glDeleteBuffers(1, &EBO);
    glUseProgram(0);
    glDeleteProgram(shaderProgram);
 
    glfwDestroyWindow(window);
    // glfw: terminate and clear all previously GLFW allocated resources
    // glfwTerminate();

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

        if (OFFSCREEN_RENDERING_COORDS) break;

        // glfw: swap buffers
        glfwSwapBuffers(window);
 
        // glfw: poll IO events (keys & mouse) 
        // (including X close window button)
        glfwPollEvents();
    }
}

double NetParam::simpleMeasureFiber(const Eigen::RowVector2f& start, const Eigen::RowVector2f& end) const {
    /*std::cout << "start/end" << std::endl;
    std::cout << start << std::endl;
    std::cout << end << std::endl;*/
    int resolution = 50;

    //Eigen::MatrixXf V_2d_ter = V_2d_;
    //f romRenderToInitCoords(V_2d_ter);

    auto findPoint = [&](const Eigen::RowVector2f& point,
                         int& id,
                         Eigen::MatrixXd& baryp){

        for (int j=0; j<F_.rows(); j++){
            Eigen::MatrixXd L;
            Eigen::MatrixXd V_2d_bis = V_2d_.cast<double>();
            Eigen::MatrixXd point_bis = point.cast<double>();
            Eigen::MatrixXd V0 =  V_2d_bis.row(F_(j, 0));
            Eigen::MatrixXd V1 =  V_2d_bis.row(F_(j, 1));
            Eigen::MatrixXd V2 =  V_2d_bis.row(F_(j, 2));
            igl::barycentric_coordinates(point_bis, V0,
                                                V1,
                                                V2,
                                                L);
            if (L.maxCoeff() > 1 || L.minCoeff() < 0){
                continue;
            }
            baryp = L;
            id = j;
            return;
        }
        
        // out of mesh
        id = -1;
        baryp = Eigen::MatrixXd(0,0);
    };

    Eigen::MatrixXf interpolated(resolution + 1, 2);
    Eigen::MatrixXd bary(resolution + 1, 3);
    Eigen::VectorXi tri_ids(resolution + 1);
    Eigen::MatrixXf interpolated_3d(resolution + 1, 2);

    for (int i=0; i<=resolution; i++){
        float t = static_cast<float>(i)/static_cast<float>(resolution);
        Eigen::RowVector2f interpoint1 = (1-t) * start + t * end;
        interpolated.row(i) = interpoint1;
        int id = -1;
        Eigen::MatrixXd baryp; 
        findPoint(interpoint1, id, baryp);
        
        tri_ids(i) = id;

        if (id != -1){
            bary.row(i) = baryp.row(0);
            interpolated_3d.row(i) = bary(i,0) * V_3d_.row(F_(tri_ids(i), 0)) 
                                + bary(i,1) * V_3d_.row(F_(tri_ids(i), 1)) 
                                + bary(i,2) * V_3d_.row(F_(tri_ids(i), 2));
        }
        else {
            interpolated_3d.row(i) = Eigen::RowVector2f(0,0);
        }
    }

    double dist = 0;
    for (int i=0; i<resolution; i++){
        if (tri_ids(i) == -1 || tri_ids(i+1) == -1) continue;
        dist += (interpolated_3d.row(i) - interpolated_3d.row(i+1)).norm();
    }

    return dist;
}

double NetParam::measureFiber(const Eigen::RowVector2d& start, const Eigen::RowVector2d& end) const {

    //glfwMakeContextCurrent(window);
    double du = start(0) - end(0);
    double dv = start(1) - end(1);

    int start_x = static_cast<int>(RENDER_WIDTH * (start(0)+1.0)/2.0);
    int start_y = static_cast<int>(RENDER_HEIGHT * (start(1)+1.0)/2.0);
    int end_x = static_cast<int>(RENDER_WIDTH * (end(0)+1.0)/2.0);
    int end_y = static_cast<int>(RENDER_HEIGHT * (end(1)+1.0)/2.0);

    std::cout << "Pixel range:" << std::endl;
    std::cout << "start: " << start(0) << " " << start(1) << std::endl;
    std::cout << "end: " << end(0) << " " << end(1) << std::endl;
    std::cout << "X: " << start_x << " " << end_x << std::endl;
    std::cout << "Y: " << start_y << " " << end_y << std::endl;

    if (start_x > end_x || start_y > end_y){
        std::cout << "ERROR, fiber given in wrong order: " << start << " " << end << std::endl;
        throw "Invalid fiber";
    }

    int axis, width_read, height_read;
    if (start_x == end_x){
        axis = 0;
        width_read = 1;
        height_read = end_y - start_y;
    }
    else if (start_y == end_y){
        axis = 1;
        height_read = 1;
        width_read = end_x - start_x;
    }
    else {
        std::cout << "ERROR, measuring fiber with invalid values: " << du << " " << dv << std::endl;
        throw "Invalid fiber";
    }

    float pixels[width_read * height_read * 3 ];
    glReadPixels(start_x, start_y, width_read, height_read, GL_RGB, GL_FLOAT, pixels);

    double dist = 0;
    for (int i=0; i<width_read * height_read-1; i+=3){
        if (pixels[3*i] == 0 && pixels[3*i+1] == 0 && pixels[3*i+2] == 0){
            std::cout << "shouldn't happen " << std::endl; // NOTE: I could also find the origin of the fiber that way... but then I have to read all black pixels, is it worth it?
            continue; // note: we might miss a pixel at the origin in 3D
        } 
        if (pixels[3*(i+1)] == 0 && pixels[3*(i+1)+1] == 0 && pixels[3*(i+1)+2] == 0) {
            std::cout << "also shouldn't happen " << std::endl;
            continue;
        }
        std::cout << "(" << pixels[3*i] << " " << pixels[3*i+1] << " " << pixels[3*i+2] << ")" << std::endl;
        dist += std::sqrt(std::pow(pixels[3*(i+1)  ] - pixels[3*i  ], 2) 
                        + std::pow(pixels[3*(i+1)+1] - pixels[3*i+1], 2) 
                        + std::pow(pixels[3*(i+1)+2] - pixels[3*i+2], 2));
    }

    std::cout << "dist = " << dist << std::endl;
    return dist;
}


void NetParam::alternativeRenderingAttempt(float start_u, float start_v, float measure_length){


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

    //GLFWwindow * window;
    if(!igl::opengl::glfw::background_window(window)){
        std::cout << "Could not initialize glfw window" << std::endl;
        glfwTerminate();
        return;
    }

    prepareShaders();

    prepareBuffers();

    /*glUniform1i(glGetUniformLocation(shaderProgram, "tex"),0);
    // Prepare texture
    GLuint in_tex;
    GLenum format;
    {
        int nc = 3;
        int w = RENDER_WIDTH;
        int h = RENDER_HEIGHT;
        format = nc==1 ? GL_RED : (nc==3 ? GL_RGB : (nc == 4 ? GL_RGBA : GL_FALSE));
        glGenTextures(1, &in_tex);
        glBindTexture(GL_TEXTURE_2D, in_tex);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
        glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
        glTexImage2D(GL_TEXTURE_2D, 0,GL_RGB, w, h, 0,format, GL_UNSIGNED_BYTE, in_data);
    }*/




    // Render to FBO http://www.opengl-tutorial.org/fr/intermediate-tutorials/tutorial-14-render-to-texture/

    // The framebuffer, which regroups 0, 1, or more textures, and 0 or 1 depth buffer.
    GLuint FramebufferName = 0;
    glGenFramebuffers(1, &FramebufferName);
    glBindFramebuffer(GL_FRAMEBUFFER, FramebufferName);

    // The texture we're going to render to
    GLuint renderedTexture;
    glGenTextures(1, &renderedTexture);

    // "Bind" the newly created texture : all future texture functions will modify this texture
    glBindTexture(GL_TEXTURE_2D, renderedTexture);

    // Give an empty image to OpenGL ( the last "0" )
    glTexImage2D(GL_TEXTURE_2D, 0,GL_RGB, RENDER_WIDTH, RENDER_HEIGHT, 0,GL_RGB, GL_UNSIGNED_BYTE, 0);

    // Poor filtering. Needed !
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);


    // Set "renderedTexture" as our colour attachement #0
    glFramebufferTexture(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, renderedTexture, 0);

    // Set the list of draw buffers.
    GLenum DrawBuffers[1] = {GL_COLOR_ATTACHMENT0};
    glDrawBuffers(1, DrawBuffers); // "1" is the size of DrawBuffers

    // Always check that our framebuffer is ok
    if(glCheckFramebufferStatus(GL_FRAMEBUFFER) != GL_FRAMEBUFFER_COMPLETE)
        std::cout << "Framebuffer ERROR" << std::endl;

    glBindFramebuffer(GL_FRAMEBUFFER, FramebufferName);
    glViewport(0,0,RENDER_WIDTH, RENDER_HEIGHT);


    glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT);

    // draw triangle
    glUseProgram(shaderProgram);
    glBindVertexArray(VAO); 

    glDrawElements(GL_TRIANGLES, 3 * n_tris_, GL_UNSIGNED_INT, 0);
    
    std::cout << RENDER_WIDTH << " x " << RENDER_HEIGHT << std::endl;
    //float start_u = 0.4, start_v=0.4, measure_length = 3.0;
    int measure_axis = 0;
    double measurement = -1;
    bool first_menu_call = true;
    Eigen::MatrixXf start_mat, end_mat;
    Eigen::RowVector2f start_point(start_u, start_v);
    Eigen::RowVector2f end_point;

    start_mat = start_point;
    //f romInitToRenderCoords(start_mat); // TODO get rid of this ?

    end_point = start_point;
    end_point(measure_axis) += measure_length;
    end_mat = end_point;
    //f romInitToRenderCoords(end_mat); // TODO get rid of this ?
    Eigen::RowVector2d temp_start, temp_end;
    temp_start.row(0) = start_mat.row(0).cast<double>();
    temp_end.row(0) = end_mat.row(0).cast<double>();
    measureFiber(temp_start, temp_end);


    glDeleteFramebuffers(1, &FramebufferName);
    glDeleteTextures(1, &renderedTexture);
    glBindFramebuffer(GL_FRAMEBUFFER, 0);
    glBindVertexArray(0);
    freeRenderingBuffers();
}



void NetParam::otherRenderingAttempt(float start_u, float start_v, float measure_length){
    GLuint tex_id, fbo_id, d_id;
    bool depth_texture = true;
    igl::opengl::init_render_to_texture(RENDER_WIDTH, RENDER_HEIGHT, depth_texture, tex_id, fbo_id, d_id);
    glBindFramebuffer(GL_FRAMEBUFFER, fbo_id);
    if(!depth_texture)
    {
      glBindRenderbuffer(GL_RENDERBUFFER, d_id);
    }
    
    //draw scene ...


    /// ---------------------

    prepareBuffers();


    glClearColor(0.1f, 0.0f, 0.0f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT);

    glActiveTexture(GL_TEXTURE0+0);
    glBindTexture(GL_TEXTURE_2D,tex_id);

    // draw triangle
    glUseProgram(shaderProgram);
    glBindVertexArray(VAO); 

    glDrawElements(GL_TRIANGLES, 3 * n_tris_, GL_UNSIGNED_INT, 0);


    std::cout << RENDER_WIDTH << " x " << RENDER_HEIGHT << std::endl;
    //float start_u = 0.4, start_v=0.4, measure_length = 3.0;
    int measure_axis = 0;
    double measurement = -1;
    bool first_menu_call = true;
    Eigen::MatrixXf start_mat, end_mat;
    Eigen::RowVector2f start_point(start_u, start_v);
    Eigen::RowVector2f end_point;

    start_mat = start_point;
    //f romInitToRenderCoords(start_mat); // TODO get rid of this ?

    end_point = start_point;
    end_point(measure_axis) += measure_length;
    end_mat = end_point;
    // f romInitToRenderCoords(end_mat); // TODO get rid of this ?
    Eigen::RowVector2d temp_start, temp_end;
    temp_start.row(0) = start_mat.row(0).cast<double>();
    temp_end.row(0) = end_mat.row(0).cast<double>();
    measureFiber(temp_start, temp_end);


    /// ---------------------
    
    //clean up
    glBindFramebuffer(GL_FRAMEBUFFER,0);
    if(!depth_texture)
    {
      glBindRenderbuffer(GL_RENDERBUFFER, 0);
    }


    // Later ...
    /*
    glActiveTexture(GL_TEXTURE0+0);
    glBindTexture(GL_TEXTURE_2D,tex_id);
    if(depth_texture)
    {
      glActiveTexture(GL_TEXTURE0+1);
      glBindTexture(GL_TEXTURE_2D,d_id);
    }*/

    

    glUseProgram(0);
}


    

void NetParam::vizBoundaryEdges(Eigen::MatrixXd& edge_begs, Eigen::MatrixXd& edge_ends) const {
    edge_begs.resize(Eb_.rows(), 2);
    edge_ends.resize(Eb_.rows(), 2);

    for (int i=0; i<Eb_.rows(); i++){
        edge_begs.row(i) = V_2d_.row(Eb_(i,0)).cast<double>();
        edge_ends.row(i) = V_2d_.row(Eb_(i,1)).cast<double>();
    }

    //f romRenderToInitCoords(edge_begs);
    //f romRenderToInitCoords(edge_ends);
}

void NetParam::computeFibers(){
    
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

std::vector<Eigen::VectorXd> NetParam::computeFibersDistortion() const {

    // In initial coordinates

    if (fiber_begs_list.size() < 2){
        std::cout << "Error, insufficient fibers, cannot color" << std::endl;
        return {};
    }

    // Eigen::MatrixXf V_2d_bis = V_2d_;
    // f romRenderToInitCoords(V_2d_bis);

    std::vector<Eigen::VectorXd> fiber_dist_per_axis;
    for (int axis=0; axis<2; axis++){
        int n_fibs = fiber_begs_list[axis].rows();
        Eigen::VectorXd fiber_dist(n_fibs);

        for (int fib=0; fib<fiber_begs_list[axis].rows(); fib++){
            Eigen::MatrixXf v1, v2;
            v1 = fiber_begs_list[axis].row(fib).cast<float>();
            v2 = fiber_ends_list[axis].row(fib).cast<float>();
            // v1r = v1;
            // v2r = v2;
            // f romInitToRenderCoords(v1r);
            // f romInitToRenderCoords(v2r);

            //double len = measureFiber(v1b, v2b);
            double len = simpleMeasureFiber(v1, v2);
            double len2 = (v1 - v2).norm();

            if (axis == 0) len *= scale_u_ / scale_v_; // TODO investigate ???
            //if (axis == 1) len *= scale_v_;

            fiber_dist(fib) = len / len2;
        }
        fiber_dist_per_axis.push_back(fiber_dist); 
    }

    return fiber_dist_per_axis;
}

std::vector<Eigen::MatrixXd> NetParam::computeFibersColor() const {
    std::vector<Eigen::VectorXd> distos = computeFibersDistortion();
    std::vector<Eigen::MatrixXd> dist_colors;
    for (int i=0; i<distos.size(); i++){
        Eigen::MatrixXd colors;
        igl::jet(distos[i], 0.75, 1.25, colors);
        dist_colors.push_back(colors);
    }

    return dist_colors;    
}

std::vector<std::vector<int>> NetParam::nearestFibers(const Eigen::RowVector2f& vertex) const {

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


    bool show_selected_id = false;
    if (show_selected_id){
        std::cout << "SELECTED FIBERS:" << std::endl;
        std::cout << "U:" << selected_fibers[0][0] << " " <<  selected_fibers[0][1] << std::endl;
        std::cout << "V:" << selected_fibers[1][0] << " " <<  selected_fibers[1][1] << std::endl;
    }

    return selected_fibers;
}

void NetParam::adjustVertices(){
    //computeFibers();
    Eigen::MatrixXf V_2d_new = V_2d_;
    // f romRenderToInitCoords(V_2d_new);
    fiber_dist = computeFibersDistortion(); // TODO do not recompute if not needed?
    //fiber_dist[1] = Eigen::VectorXd::Ones(fiber_dist[1].rows());

    for (int v_id=0; v_id<V_2d_.rows(); v_id++){
        // Find 4 closest fibers
        
        std::vector<std::vector<int>> fibs = nearestFibers(V_2d_.row(v_id));

        std::cout << fibs[0][0] << " " << fibs[0][1] << " " << fibs[1][0] << " " << fibs[1][1] << " " << std::endl;

        for (int axis=0; axis<2; axis++){

            //if (axis == 1) continue; // TEST 

            float off_axis = 0;
            int valid_fibs = 0;

            for (int i: fibs[axis]){
                if (i < 0) continue;
                valid_fibs ++;
                float t1 = V_2d_(v_id, (axis+1) % 2);
                float f1 = fiber_begs_list[axis](i, (axis+1) % 2);
                float f2 = fiber_ends_list[axis](i, (axis+1) % 2);
                float m = (f1 + f2) / 2.0;
                float d = t1 - m;

                // less weight when far on axis. Assumes we're at fiber scale
                float coeff = 1.0 - std::fabs(V_2d_(v_id, axis) - fiber_begs_list[axis](i, axis));
                
                if (t1 < f1 || t1 > f2) std::cout << "ERROR: vertex out of fiber range" << std::endl;
                if (coeff < 0 || coeff > 1) std::cout << "ERROR: invalid coeff, are we at fiber scale?" << std::endl;
                std::cout << "fiber range " << f1 << " -> " << t1 << " -> " << f2 << std::endl;
                std::cout << "fiber_dist[axis](i) " << fiber_dist[axis](i) << std::endl;
                off_axis = coeff * d * (static_cast<float>(fiber_dist[axis](i)) - 1.0);
                
            }

            std::cout << "off " << off_axis << std::endl;
            std::cout << "valid_fibs " << valid_fibs << std::endl;

            if (valid_fibs == 0) continue;
            off_axis /= static_cast<float>(valid_fibs);

            V_2d_new(v_id, (axis+1)%2) += off_axis;
        }
        

    }

    // f romInitToRenderCoords(V_2d_new);
    V_2d_ = V_2d_new;

    // overwrite data that is now invalid
    fiber_begs_list = {}; 
    fiber_ends_list = {}; 
    fiber_dist = {};
}

/*void NetParam::adjustVertices(){

    Eigen::MatrixXf V_2d_new = V_2d_;
    // f romRenderToInitCoords(V_2d_new);
    fiber_dist = computeFibersDistortion(); // TODO do not recompute if not needed?
    fiber_dist[1] = Eigen::VectorXd::Ones(fiber_dist[1].rows());

    for (int v_id=0; v_id<V_2d_.rows(); v_id++){
        // Find 4 closest fibers
        
        std::vector<std::vector<int>> fibs = nearestFibers(V_2d_.row(v_id));

        for (int axis=0; axis<2; axis++){
            float off_axis = 0;
            int valid_fibs = 0;

            for (int i: fibs[axis]){
                if (i < 0) continue;
                valid_fibs ++;
                float t1 = V_2d_((axis+1) % 2);
                float f1 = fiber_begs_list[axis](i, (axis+1) % 2);
                float f2 = fiber_ends_list[axis](i, (axis+1) % 2);
                float m = (f1 + f2) / 2.0;
                float d = t1 - m;
                std::cout << "fiber_dist[axis](i) " << fiber_dist[axis](i) << std::endl;
                off_axis = d * (static_cast<float>(fiber_dist[axis](i)) - 1.0);
                
            }

            std::cout << "off " << off_axis << std::endl;
            std::cout << "valid_fibs " << valid_fibs << std::endl;

            if (valid_fibs == 0) continue;
            off_axis /= static_cast<float>(valid_fibs);

            V_2d_new(v_id, (axis+1)%2) += off_axis;
        }
        

    }

    // f romInitToRenderCoords(V_2d_new);
    V_2d_ = V_2d_new;

    // overwrite data that is now invalid
    fiber_begs_list = {}; 
    fiber_ends_list = {}; 
    fiber_dist = {};
}*/