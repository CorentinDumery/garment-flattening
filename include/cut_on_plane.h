/**
 * @author Corentin Dumery
 * @brief Cut a mesh on a plane. Remeshes around the plane intersection
 * so that the cut matches perfectly.
 * @date 2022-07-04
 * 
 */

#pragma once
#include <Eigen/Core>
#include <vector>
#include <iostream>

#include <igl/remove_duplicate_vertices.h>
#include <igl/facet_components.h>
#include <igl/cut_mesh.h>
#include <igl/remove_unreferenced.h>

bool isPointOnLine(const Eigen::Vector3d& point1, const Eigen::Vector3d& point2, const Eigen::Vector3d& point3){
    Eigen::Vector3d lineDirection = point2 - point1;
    Eigen::Vector3d point1ToPoint3 = point3 - point1;

    // Check if the vectors are collinear
    double crossProductNorm = lineDirection.cross(point1ToPoint3).norm();
    double epsilon = 1e-6; // A small tolerance for numerical stability

    return crossProductNorm < epsilon;
}

// Function to cut a triangle mesh along a plane with constant Z value
void precutMeshOnPlane(const std::vector<Eigen::Vector3d>& vertices, 
                    const std::vector<Eigen::Vector3i>& triangles,
                    double z_cut,
                    std::vector<Eigen::Vector3d>& cutVertices,
                    std::vector<Eigen::Vector3d>& V_cut,
                    std::vector<Eigen::Vector3i>& F_cut,
                    std::vector<Eigen::Vector3i>& cuts,
                    int axis = 0){
                        
    //std::vector<Eigen::Vector3d> cutVertices;
    std::vector<Eigen::Vector3i> cutTriangles;

    Eigen::VectorXi tri_intersects = Eigen::VectorXi::Zero(triangles.size());
    std::vector<std::vector<Eigen::Vector3d>> inter_points_per_tri;
    std::vector<std::vector<int>> lone_vertex;

    // Step 2: Determine intersecting triangles
    int tri_id = -1;
    for (const auto& triangle : triangles){
        tri_id ++;
        bool isIntersecting = false;
        std::vector<int> edgesOnOppositeSides;

        std::vector<Eigen::Vector3d> inter_points;

        for (int i = 0; i < 3; i++)
        {
            double z1 = vertices[triangle[i]](axis);
            double z2 = vertices[triangle[(i + 1) % 3]](axis);

            // Check if vertex is on opposite sides of the cutting plane
            if ((z1 < z_cut && z2 > z_cut) ||
                (z1 > z_cut && z2 < z_cut))
            {
                isIntersecting = true;
                tri_intersects(tri_id) = 1;
                edgesOnOppositeSides.push_back(i);
            }
        }

        // Step 3: Compute intersection points
        if (isIntersecting){
            for (int i=0; i<edgesOnOppositeSides.size(); i++){
                int id1 = triangle[edgesOnOppositeSides[i]];
                int id2 = triangle[(edgesOnOppositeSides[i] + 1) % 3];
                Eigen::Vector3d v1 = vertices[id1];
                Eigen::Vector3d v2 = vertices[id2];
                double t = (z_cut - v1(axis)) / (v2(axis) - v1(axis));
                Eigen::Vector3d intersectionPoint = v1 + t * (v2 - v1);
                cutVertices.push_back(intersectionPoint);
                inter_points.push_back(intersectionPoint);
            }
        }

        // store lone vertex

        if (isIntersecting){
            std::vector<int> count = {0, 0, 0};
            count[edgesOnOppositeSides[0]] ++;
            count[(edgesOnOppositeSides[0] + 1) % 3] ++;
            count[edgesOnOppositeSides[1]] ++;
            count[(edgesOnOppositeSides[1] + 1) % 3] ++;
            if (count[0] == 2) lone_vertex.push_back({0});
            else if (count[1] == 2) lone_vertex.push_back({1});
            else if (count[2] == 2) lone_vertex.push_back({2});
            else std::cout << "Algo error" << std::endl;
        } 
        else lone_vertex.push_back({});

        // Add intersecting triangles to the new mesh
        if (isIntersecting){
            cutTriangles.push_back(triangle);
        }

        inter_points_per_tri.push_back(inter_points);
    }

    

    
    V_cut = vertices;

    /* For a triangle 012 intersecting the cutting plane,
       consider vertices 3 and 4 on the edges [01] and [02] 
       respectively, and their middlepoint 5.
       
          0
          /\
         /  \
        /    \
      3/__5___\4
      /        \
     /          \
    /____________\
    1             2

    */

    for (int tri=0; tri<triangles.size(); tri++){
        if (!tri_intersects(tri)){
            F_cut.push_back(triangles[tri]);
            cuts.push_back(Eigen::Vector3i::Zero());
        }
        else {
            // get the two vertices of interest
            int loner_id = lone_vertex[tri][0];
            int v0_id = triangles[tri](loner_id);
            int v1_id = triangles[tri]((loner_id + 1) % 3);
            int v2_id = triangles[tri]((loner_id + 2) % 3);
            int v3_id = V_cut.size(); // these will be added on top
            int v4_id = V_cut.size() + 1;
            int v5_id = V_cut.size() + 2;
            Eigen::Vector3d v3 = inter_points_per_tri[tri][0];
            Eigen::Vector3d v4 = inter_points_per_tri[tri][1];
            if (isPointOnLine(V_cut[v0_id], V_cut[v1_id], v4)){ // make sure we match the scheme, otherwise, swap v3 and v4
                Eigen::Vector3d temp = v4; 
                v4 = v3;
                v3 = temp;
            }
            Eigen::Vector3d v5 = (v3 + v4) / 2;
            V_cut.push_back(v3);
            V_cut.push_back(v4);
            V_cut.push_back(v5);

            F_cut.push_back(Eigen::Vector3i(v0_id, v3_id, v5_id));
            cuts.push_back(Eigen::Vector3i(0, 1, 0));
            
            F_cut.push_back(Eigen::Vector3i(v0_id, v5_id, v4_id));
            cuts.push_back(Eigen::Vector3i(0, 1, 0));

            F_cut.push_back(Eigen::Vector3i(v1_id, v5_id, v3_id));
            cuts.push_back(Eigen::Vector3i(0, 1, 0));

            F_cut.push_back(Eigen::Vector3i(v2_id, v4_id, v5_id));
            cuts.push_back(Eigen::Vector3i(0, 1, 0));

            F_cut.push_back(Eigen::Vector3i(v1_id, v2_id, v5_id));
            cuts.push_back(Eigen::Vector3i(0, 0, 0));

        }
    }
}

void cutMeshOnPlane(const Eigen::MatrixXd& V, 
                    const Eigen::MatrixXi& F,
                    std::vector<Eigen::MatrixXd>& V_list,
                    std::vector<Eigen::MatrixXi>& F_list, 
                    Eigen::VectorXi& cut0, 
                    Eigen::VectorXi& cut1, 
                    int axis = 2,
                    double z_cut = 0){
    std::vector<Eigen::Vector3d> cutVertices;
    std::vector<Eigen::Vector3d> V_cut;
    std::vector<Eigen::Vector3i> F_cut;
    std::vector<Eigen::Vector3i> cuts;

    std::vector<Eigen::Vector3d> V_vec;
    std::vector<Eigen::Vector3i> F_vec;
    for (int i=0; i<V.rows(); i++) V_vec.push_back(V.row(i));
    for (int i=0; i<F.rows(); i++) F_vec.push_back(F.row(i));
    precutMeshOnPlane(V_vec, F_vec, z_cut, cutVertices, V_cut, F_cut, cuts, axis);
    
    Eigen::MatrixXd points_cut = Eigen::MatrixXd(cutVertices.size(), 3);
    for (int i=0; i<cutVertices.size(); i++) points_cut.row(i) = cutVertices[i];

    Eigen::MatrixXd V_cut2 = Eigen::MatrixXd(V_cut.size(), 3);
    for (int i=0; i<V_cut.size(); i++) V_cut2.row(i) = V_cut[i];

    Eigen::MatrixXi F_cut2 = Eigen::MatrixXi(F_cut.size(), 3);
    for (int i=0; i<F_cut.size(); i++) F_cut2.row(i) = F_cut[i];

    Eigen::MatrixXi cuts2 = Eigen::MatrixXi(cuts.size(), 3);
    for (int i=0; i<cuts.size(); i++) cuts2.row(i) = cuts[i];

    Eigen::MatrixXd SV;
    Eigen::VectorXi SVI, SVJ;
    igl::remove_duplicate_vertices(V_cut2, 1e-7, SV, SVI, SVJ);
    // remap faces
    Eigen::MatrixXi SF(F_cut2.rows(), F_cut2.cols());
    for (int i=0; i<F_cut2.rows(); i++){
        SF(i,0) = SVJ(F_cut2(i,0));
        SF(i,1) = SVJ(F_cut2(i,1));
        SF(i,2) = SVJ(F_cut2(i,2));
    } 

    Eigen::MatrixXd V_cut3;
    Eigen::MatrixXi F_cut3;
    Eigen::VectorXi corres;
    igl::cut_mesh<Eigen::Matrix<double, -1, -1, 0, -1, -1>, // ambiguous if template not specified
                  Eigen::Matrix<int, -1, -1, 0, -1, -1>, 
                  Eigen::Matrix<int, -1, -1, 0, -1, -1>, 
                  Eigen::Matrix<int, -1, 1, 0, -1, 1> >(SV, SF, cuts2, V_cut3, F_cut3, corres);

    // Find pairs of new vertices that map to the same old vertex 
    Eigen::VectorXi pointed_at_by = Eigen::VectorXi::Constant(V_cut3.rows(), -1); 
    std::vector<std::pair<int,int>> pairs;
    for (int i=0; i<corres.rows(); i++){
        if (pointed_at_by(corres[i]) != -1){
            // vertex already assigned, so it has been duplicated
            pairs.push_back(std::make_pair(corres[i], i));
        }
        else {
            pointed_at_by(corres[i]) = i;
        }
    }

    Eigen::VectorXi C;
    igl::facet_components(F_cut3, C);

    if (C.maxCoeff() > 1) {
        std::cout << "ERROR while cutting mesh: more than two components, case not handled." << std::endl; 
    }

    // Extend per-face component id vector C to per vertex component C_v 
    Eigen::VectorXi C_v = Eigen::VectorXi::Constant(V_cut3.rows(), -1);
    for (int f_id=0; f_id<C.rows(); f_id++){
        C_v(F_cut3(f_id, 0)) = C(f_id);
        C_v(F_cut3(f_id, 1)) = C(f_id);
        C_v(F_cut3(f_id, 2)) = C(f_id);
    }


    cut0 = Eigen::VectorXi::Constant(pairs.size(), -1);
    cut1 = Eigen::VectorXi::Constant(pairs.size(), -1);
    for (int i=0; i<pairs.size(); i++){
        if (C_v(pairs[i].first) == 0){
            cut0(i) = pairs[i].first;
            cut1(i) = pairs[i].second;
        }
        else {
            cut0(i) = pairs[i].second;
            cut1(i) = pairs[i].first;    
        }
    }

    int size_comp0 = C.unaryExpr([](double val) { return val == 0 ? 1 : 0; }).sum();
    int size_comp1 = C.unaryExpr([](double val) { return val == 1 ? 1 : 0; }).sum();
    
    Eigen::MatrixXi F_comp0(size_comp0, 3); 
    int curr = 0;
    for (int i=0; i<F_cut3.rows(); i++){
        if (C(i) == 0) {
            F_comp0.row(curr) = F_cut3.row(i);
            curr ++;
        }
    }

    Eigen::MatrixXi F_comp1(size_comp1, 3); 
    curr = 0;
    for (int i=0; i<F_cut3.rows(); i++){
        if (C(i) == 1) {
            F_comp1.row(curr) = F_cut3.row(i);
            curr ++;
        }
    }

    Eigen::MatrixXd V0, V1;
    Eigen::MatrixXi F0, F1;
    Eigen::VectorXi I0, I1;
    igl::remove_unreferenced(V_cut3, F_comp0, V0, F0, I0);
    igl::remove_unreferenced(V_cut3, F_comp1, V1, F1, I1);

    V_list = {V0, V1};
    F_list = {F0, F1}; // TODO potentially more components here

    for (int i=0; i<cut0.rows(); i++){ // TODO make more compact?
        cut0(i) = I0(cut0(i));
        cut1(i) = I1(cut1(i));
    }
}