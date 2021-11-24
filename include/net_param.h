
#include <Eigen/Core>

#include <glad/glad.h>
#include <GLFW/glfw3.h>

#include <vector>

class NetParam {

public:
    NetParam(const Eigen::MatrixXi& F,
             const Eigen::MatrixXf& V_3d,
             const Eigen::MatrixXf& V_2d);

    void initializeRendering();
    void render();
    void freeRenderingBuffers();

    void vizBoundaryEdges(Eigen::MatrixXd& edge_begs, Eigen::MatrixXd& edge_ends);

    void computeFibers();
    std::vector<std::vector<Eigen::MatrixXd>> vizFibers(){return {fiber_begs_list, fiber_ends_list};}

private:
    
    
    double measureFiber();
    void modifyParam();

    // CONFIG
    const unsigned int RENDER_WIDTH = 800;
    const unsigned int RENDER_HEIGHT = 600;

    // Set by constructor
    const Eigen::MatrixXi F_;
    const Eigen::MatrixXf V_3d_;
    Eigen::MatrixXf V_2d_; // Note: ideally we'd get rid of that and directly work on something closer to the rendering buffers
    Eigen::MatrixXi Eb_; // boundary edges
    Eigen::MatrixXi TT_; // triangle-triangle adjacency
    int n_tris_;

    // Rendering buffers/variables
    unsigned int VBO, VAO, EBO;
    unsigned int shaderProgram; // TODO naming
    GLFWwindow* window;

    // Net fiber variables
    std::vector<Eigen::MatrixXd> fiber_begs_list; // Fiber viz
    std::vector<Eigen::MatrixXd> fiber_ends_list; 
};