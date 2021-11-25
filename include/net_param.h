
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

    void computeFibers();

    void vizBoundaryEdges(Eigen::MatrixXd& edge_begs, Eigen::MatrixXd& edge_ends) const;
    std::vector<std::vector<Eigen::MatrixXd>> vizFibers() const {return {fiber_begs_list, fiber_ends_list};}

    double measureFiber(const Eigen::RowVector2d& start, const Eigen::RowVector2d& end) const;

    // 2D mesh is in [-1, 1] for rendering, this function converts back to initial coords
    // Note: we could keep normal coords and just send min/max to shader instead
    void fromInitToRenderCoords(Eigen::MatrixXf& points) const;
    void fromRenderToInitCoords(Eigen::MatrixXf& points) const;
    void fromRenderToInitCoords(Eigen::MatrixXd& points) const;

private:
    
    void modifyParam();

    

    // CONFIG
    const unsigned int RENDER_WIDTH = 300;

    // Set by constructor
    const Eigen::MatrixXi F_;
    Eigen::MatrixXf V_3d_;
    Eigen::MatrixXf V_2d_; // Note: ideally we'd get rid of that and directly work on something closer to the rendering buffers
    Eigen::MatrixXi Eb_; // boundary edges
    Eigen::MatrixXi TT_; // triangle-triangle adjacency
    int n_tris_;
    float scale_u_;
    float scale_v_;
    float off_u_;
    float off_v_;
    unsigned int RENDER_HEIGHT; // based on RENDER_WIDTH and initial X/Y ratio in V_2d_ 

    // Rendering buffers/variables
    unsigned int VBO, VAO, EBO;
    unsigned int shaderProgram; // TODO naming
    GLFWwindow* window;

    // Net fiber variables
    std::vector<Eigen::MatrixXd> fiber_begs_list; // Fiber viz
    std::vector<Eigen::MatrixXd> fiber_ends_list; 
};