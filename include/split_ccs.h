#include <Eigen/Core>
#include <igl/vertex_components.h>

void splitMeshIntoCCs(const Eigen::MatrixXd& V,
                      const Eigen::MatrixXi& F,
                      std::vector<Eigen::MatrixXd>& V_comps,
                      std::vector<Eigen::MatrixXi>& F_comps,
                      Eigen::VectorXi& vertex_components,
                      Eigen::VectorXi& face_components,
                      std::vector<Eigen::VectorXi>& v_maps,
                      std::vector<Eigen::VectorXi>& f_maps){
    
    // for a vertex i in V.row(i),
    // it is moved to
    // V_comps[vertex_components(i)].row(v_maps[vertex_components(i)](i))
    
    V_comps = {};
    F_comps = {};
    v_maps = {};
    f_maps = {};

    // Compute vertex and face components
    igl::vertex_components(F, vertex_components);
    face_components.resize(F.rows());
    for (int j = 0; j < F.rows(); j++) {
        face_components(j) = vertex_components(F(j,0));
    }

    // Isolate mesh for each component
    int num_components = vertex_components.maxCoeff() + 1;
    for (int i = 0; i < num_components; i++) {
        Eigen::VectorXi component_vertices = (vertex_components.array() == i).select(Eigen::VectorXi::Ones(V.rows()), Eigen::VectorXi::Zero(V.rows()));
        Eigen::VectorXi component_faces = (face_components.array() == i).select(Eigen::VectorXi::Ones(F.rows()), Eigen::VectorXi::Zero(F.rows()));
        
        std::cout << "Component " << i << " vertices: " << component_vertices.sum() << std::endl;
        std::cout << "Component " << i << " faces: " << component_faces.sum() << std::endl;

        // Compute vertex component
        Eigen::MatrixXd V_component(component_vertices.sum(), V.cols());
        Eigen::VectorXi v_component_map(V.rows()); // TODO initialize to -1?
        int count = 0;
        for (int j = 0; j < V.rows(); j++) {
            if (component_vertices(j)) {
                V_component.row(count) = V.row(j);
                v_component_map(j) = count;
                count ++;
            }
        }
        if (count != V_component.rows()){
            std::cout << "Error non matching v sizes:" << count << ", " << V_component.rows() << std::endl;
        }
        
        // Compute face component
        Eigen::MatrixXi F_component(component_faces.sum(), F.cols());
        Eigen::VectorXi f_component_map(F.rows());
        count = 0;
        for (int j = 0; j < F.rows(); j++) {
            if (component_faces(j)) {   
                F_component.row(count) = F.row(j);

                for (int k = 0; k < F.cols(); k++) {
                    F_component(count, k) = v_component_map(F(j, k));
                }
                f_component_map(j) = count;
                count ++;
            }
        }
        if (count != F_component.rows()){
            std::cout << "Error non matching f sizes:" << count << ", " << F_component.rows() << std::endl;
        }

        V_comps.push_back(V_component);
        F_comps.push_back(F_component);
        v_maps.push_back(v_component_map);
        f_maps.push_back(f_component_map);
    }
}


void mergeCCsBack(const std::vector<Eigen::MatrixXd>& V_comps, 
                  const Eigen::VectorXi& vertex_components,
                  const std::vector<Eigen::VectorXi>& v_maps,
                  int n_vs,
                  Eigen::MatrixXd& V_merged){
    
    V_merged.resize(n_vs, 3);
    for (int i=0; i<n_vs; i++){
        V_merged.row(i) = V_comps[vertex_components(i)].row(v_maps[vertex_components(i)](i));
    }
}