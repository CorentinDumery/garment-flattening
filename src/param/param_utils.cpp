#include "param/param_utils.h"

#include <Eigen/Geometry>
#include <Eigen/SVD>
#include <Eigen/LU>
#include <iostream>
#include <cmath>
#include <igl/doublearea.h>

#include <igl/lscm.h>

#include <igl/remove_unreferenced.h>
#include <igl/remove_duplicate_vertices.h>

// Required for ARAP
#include <igl/arap.h>
#include <igl/boundary_loop.h>
#include <igl/harmonic.h>
#include <igl/map_vertices_to_circle.h>

// Required for SCAF
#ifdef USE_SCAF_PARAM
#include <igl/triangle/scaf.h>
#include <igl/arap.h>
#include <igl/boundary_loop.h>
#include <igl/harmonic.h>
#include <igl/map_vertices_to_circle.h>
#include <igl/MappingEnergyType.h>
#include <igl/doublearea.h>
#include <igl/PI.h>
#include <igl/flipped_triangles.h>
#include <igl/topological_hole_fill.h>
#endif

void meshCleanup(Eigen::MatrixXd& V_3d, Eigen::MatrixXi& F){
    Eigen::MatrixXd V_3db(V_3d.rows(), 3);
    V_3db.col(0) = V_3d.col(0);
    V_3db.col(1) = V_3d.col(1);
    V_3db.col(2) = V_3d.col(2);
    V_3d = V_3db;

    Eigen::MatrixXd NV;
    Eigen::MatrixXi NF;
    Eigen::VectorXi I, J;
    igl::remove_unreferenced(V_3d, F, NV, NF, I, J);
    std::cout << "Removing unref vertices: " << V_3d.rows() << " down to " << NV.rows() << std::endl;
    V_3d = NV;
    F = NF;

    // DISABLED BECAUSE THE FOLLOWING OPERATIONS MODIFY F'S SIZE
    /*
    Eigen::MatrixXd SV;
    Eigen::VectorXi SVI, SVJ;
    igl::remove_duplicate_vertices(V_3d, 1e-5, SV, SVI, SVJ);
    // remap faces
    Eigen::MatrixXi SF(F.rows(), F.cols());
    for (int i=0; i<F.rows(); i++){
        SF(i,0) = SVJ(F(i,0));
        SF(i,1) = SVJ(F(i,1));
        SF(i,2) = SVJ(F(i,2));
    } 

    V_3d = SV;
    F = SF;

    
    std::vector<int> degenerateRows;
    degenerateRows.reserve(F.rows());

    // Identify degenerate triangles (rows with two or more equal indices)
    for (int i = 0; i < F.rows(); ++i) {
        if (F(i, 0) == F(i, 1) || F(i, 0) == F(i, 2) || F(i, 1) == F(i, 2)) {
            degenerateRows.push_back(i);
        }
    }

    // Create the new matrix F2 without degenerate triangles
    int rowsToRemove = degenerateRows.size();
    int rowsOriginal = F.rows();
    int rowsNew = rowsOriginal - rowsToRemove;
    Eigen::MatrixXi F2(rowsNew, 3);

    int srcRow = 0;
    for (int dstRow = 0; dstRow < rowsNew; ++dstRow) {
        while (srcRow < rowsOriginal && std::binary_search(degenerateRows.begin(), degenerateRows.end(), srcRow)) {
            // Skip degenerate rows
            ++srcRow;
        }
        if (srcRow > rowsOriginal) break;
        F2.row(dstRow) = F.row(srcRow);
        ++srcRow;
    }

    F = F2;*/
}


void procrustes(const Eigen::MatrixXd& points1, // to
               const Eigen::MatrixXd& points2, // from
               Eigen::MatrixXd& R_est,
               Eigen::VectorXd& T_est){

    // https://igl.ethz.ch/projects/ARAP/svd_rot.pdf
    // https://math.stackexchange.com/questions/849217/estimate-rotation-and-translation-from-two-sets-of-points-in-different-coordinat
    // https://en.wikipedia.org/wiki/Procrustes_analysis

    Eigen::MatrixXd points1t = points1.transpose();
    Eigen::MatrixXd points2t = points2.transpose();

    Eigen::VectorXd pb = points1t.rowwise().mean();
    Eigen::VectorXd qb = points2t.rowwise().mean();

    Eigen::MatrixXd X = (points1t.colwise() - pb);
    Eigen::MatrixXd Y = (points2t.colwise() - qb);
    Eigen::MatrixXd S = X * Y.transpose();
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(S, Eigen::ComputeThinU | Eigen::ComputeThinV);
    Eigen::MatrixXd sigma = Eigen::MatrixXd::Identity(svd.matrixU().cols(), svd.matrixV().cols());
    sigma(sigma.rows() - 1, sigma.cols() - 1) = -(svd.matrixV() * svd.matrixU().transpose()).determinant();
    R_est = svd.matrixV() * sigma * svd.matrixU().transpose();
    T_est = qb - R_est * pb;
}

// Move triangle from 3D (Z != 0) to 2D (Z = 0), with
// same lengths but arbitrary orientation
Eigen::MatrixXd move3Dto2D(const Eigen::MatrixXd& V_tri){
    Eigen::MatrixXd V_2d(3,3); // TODO Matrix3d
    // First, move 3D triangle to 2D plane
    // V_2d: put A in (0,0), B in (0, |AB|), and find C 
    double r0 = (V_tri.row(1) - V_tri.row(0)).norm();
    double r1 = (V_tri.row(2) - V_tri.row(0)).norm();
    double r2 = (V_tri.row(2) - V_tri.row(1)).norm();
    V_2d.row(0) = Eigen::RowVector3d(0, 0, 0);
    V_2d.row(1) = Eigen::RowVector3d(r0, 0, 0);
    double CAB_angle = std::acos((r0*r0 + r1*r1 - r2*r2)/(2*r0*r1));
    double l1 = r1 * std::cos(CAB_angle);
    double h = l1 * std::tan(CAB_angle);
    V_2d.row(2) = Eigen::RowVector3d(l1, h, 0);

    r0 = (V_2d.row(1) - V_2d.row(0)).norm();
    r1 = (V_2d.row(2) - V_2d.row(0)).norm();
    r2 = (V_2d.row(2) - V_2d.row(1)).norm();
    if ( std::fabs(r0 - (V_tri.row(1) - V_tri.row(0)).norm()) // just checking...
        +std::fabs(r1 - (V_tri.row(2) - V_tri.row(0)).norm())
        +std::fabs(r2 - (V_tri.row(2) - V_tri.row(1)).norm()) > 0.0001){
        std::cout << "ERROR, flat triangle is different:" << std::endl;
        std::cout << r0 << " vs " << (V_tri.row(1) - V_tri.row(0)).norm() << std::endl;
    }

    return V_2d;
}

// TODO Matrix3d ?
Eigen::MatrixXd makeTriPoints(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, int f_id){
    Eigen::MatrixXd p(3,3);
    p.row(0) = V.row(F(f_id,0));
    p.row(1) = V.row(F(f_id,1));
    p.row(2) = V.row(F(f_id,2));
    return p;
};

void makeTriPoints(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, int f_id, Eigen::MatrixXd& out){
    out.row(0) = V.row(F(f_id,0));
    out.row(1) = V.row(F(f_id,1));
    out.row(2) = V.row(F(f_id,2));
}

Eigen::VectorXd vertices2dToVector(const Eigen::MatrixXd& V){
    Eigen::VectorXd res(2 * V.rows()); // TODO faster?
    for (int i=0; i<V.rows(); i++){
        res(2 * i) = V(i,0);
        res(2 * i + 1) = V(i,1);
    }
    return res;
}

Eigen::Vector3d barycentricCoords(const Eigen::RowVector3d& p, const Eigen::RowVector3d& a, 
                                         const Eigen::RowVector3d& b, const Eigen::RowVector3d& c){
    Eigen::RowVector3d v0 = b - a;
    Eigen::RowVector3d v1 = c - a;
    Eigen::RowVector3d v2 = p - a;
    float d00 = v0.dot(v0);
    float d01 = v0.dot(v1);
    float d11 = v1.dot(v1);
    float d20 = v2.dot(v0);
    float d21 = v2.dot(v1);
    float denom = d00 * d11 - d01 * d01;
    double v = (d11 * d20 - d01 * d21) / denom;
    double w = (d00 * d21 - d01 * d20) / denom;
    double u = 1.0f - v - w;
    return Eigen::Vector3d(u, v, w);
}

Eigen::Matrix3d computeRotation(const Eigen::RowVector3d& from,
                                const Eigen::RowVector3d& to){
    // There might already be something like this in Eigen? Couldn't find it
    // https://math.stackexchange.com/questions/180418/calculate-rotation-matrix-to-align-vector-a-to-vector-b-in-3d
    Eigen::RowVector3d a = from.normalized();
    Eigen::RowVector3d b = to.normalized();
    if (a==b || a == -b) std::cout << "ERROR: case not handled in computeRotation" << std::endl;
    Eigen::RowVector3d v = a.cross(b);
    double s = v.norm();
    double c = a.dot(b);
    Eigen::Matrix3d vs;
    vs <<     0, -v[2],  v[1], 
           v[2],     0, -v[0], 
          -v[1],  v[0],     0;

    return Eigen::Matrix3d::Identity() + vs + vs * vs * 1.0 / (1.0+c);
}

Eigen::MatrixXd paramARAP(const Eigen::MatrixXd& V_3d, const Eigen::MatrixXi& F,
                          const Eigen::VectorXi& bnd){

    Eigen::MatrixXd V_2d, V_3db;
    V_3db = V_3d;
    double scale = V_3d.maxCoeff();
    V_3db = V_3db.array() / scale; 

    // Compute the initial solution for ARAP (harmonic parametrization)
    Eigen::MatrixXd bnd_uv;
    igl::map_vertices_to_circle(V_3db, bnd, bnd_uv);

    Eigen::MatrixXd initial_guess;
    igl::harmonic(V_3db, F, bnd, bnd_uv, 1, initial_guess);

    // Add dynamic regularization to avoid to specify boundary conditions
    igl::ARAPData arap_data;
    arap_data.with_dynamics = true;
    Eigen::VectorXi b  = Eigen::VectorXi::Zero(0);
    Eigen::MatrixXd bc = Eigen::MatrixXd::Zero(0,0);

    // Initialize ARAP
    arap_data.max_iter = 100;
    // 2 means that we're going to *solve* in 2d
    arap_precomputation(V_3db, F, 2, b, arap_data);

    // Solve arap using the harmonic map as initial guess
    V_2d = initial_guess;
    arap_solve(bc, arap_data, V_2d);


    Eigen::MatrixXd V_2db = Eigen::MatrixXd::Zero(V_2d.rows(), 3);
    V_2db.col(0) = V_2d.col(0);
    V_2db.col(1) = V_2d.col(1);
    
    return V_2db.array() * scale;
}

Eigen::MatrixXd paramSCAF(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F,
                          const Eigen::VectorXi& bnd){
    
    std::cout << "ERROR: SCAF parameterization disabled" << std::endl;
    std::cout << "(Remember to add igl::triangle to cmake project if you want to enable it)" << std::endl;

    // NOT WORKING: runs into Numerical issue?

    // https://github.com/libigl/libigl/blob/main/tutorial/710_SCAF/main.cpp
    /*igl::triangle::SCAFData scaf_data;
    Eigen::MatrixXd bnd_uv, uv_init;
    Eigen::VectorXd M;
    igl::doublearea(V, F, M);

    float uv_scale = 1.0f;
    igl::map_vertices_to_circle(V, bnd, bnd_uv);
    bnd_uv *= sqrt(M.sum() / (2 * igl::PI));
    
    if (bnd.rows() == V.rows()) { // case: all vertex on boundary
        uv_init.resize(V.rows(), 2);
        for (int i = 0; i < bnd.rows(); i++)
        uv_init.row(bnd(i)) = bnd_uv.row(i);
    }
    else {
        igl::harmonic(V, F, bnd, bnd_uv, 1, uv_init);
        if (igl::flipped_triangles(uv_init, F).size() != 0)
        igl::harmonic(F, bnd, bnd_uv, 1, uv_init); // fallback uniform laplacian
    }
    
    
    Eigen::VectorXi b; Eigen::MatrixXd bc;
    igl::triangle::scaf_precompute(V, F, uv_init, scaf_data, igl::MappingEnergyType::SYMMETRIC_DIRICHLET, b, bc, 0);
    igl::triangle::scaf_solve(scaf_data, 1);
    Eigen::MatrixXd V_uv = uv_scale * scaf_data.w_uv.topRows(V.rows());
    return V_uv;*/
}


Eigen::MatrixXd paramLSCM(const Eigen::MatrixXd& V_3d, const Eigen::MatrixXi& F,
                          const Eigen::VectorXi& bnd){
    Eigen::VectorXi b(2,1);
    b(0) = bnd(0);
    b(1) = bnd(bnd.size()/2);
    Eigen::MatrixXd bc(2, 2);

    double dist = (V_3d.row(b(0)) - V_3d.row(b(1))).norm();

    bc<<0,0,dist,0;

    Eigen::MatrixXd V_2d;
    igl::lscm(V_3d,F,b,bc,V_2d);
    Eigen::MatrixXd V_2db = Eigen::MatrixXd::Zero(V_2d.rows(), 3);
    V_2db.col(0) = V_2d.col(0);
    V_2db.col(1) = V_2d.col(1);
    return V_2db;
}

Eigen::MatrixXd paramLSCMwithConstraint(const Eigen::MatrixXd& V_3d, const Eigen::MatrixXi& F,
                                        int v1_id, double v1_u, double v1_v,
                                        int v2_id, double v2_u, double v2_v){
    Eigen::VectorXi b(2,1);
    b(0) = v1_id;
    b(1) = v2_id;
    Eigen::MatrixXd bc(2, 2);
    std::cout << b << std::endl;
    //bc<< v1_u, v2_u, v1_v,  v2_v;
    bc<< v1_u, v1_v, v2_u, v2_v;
    
    //std::cout << bc << std::endl;
    /*
    bc.row(0) = Eigen::RowVector2d(v1_u, v1_v);
    bc.row(1) = Eigen::RowVector2d(v2_u, v2_v);


    bc.row(0) = Eigen::RowVector2d(v1_u, v2_u);
    bc.row(1) = Eigen::RowVector2d(v1_v, v2_v);
    bc.row(2) = Eigen::RowVector2d(0, 0);*/

    std::cout << bc << std::endl;

    Eigen::MatrixXd V_2d;
    igl::lscm(V_3d, F, b, bc, V_2d);
    Eigen::MatrixXd V_2db = Eigen::MatrixXd::Zero(V_2d.rows(), 3);
    V_2db.col(0) = V_2d.col(0);
    V_2db.col(1) = V_2d.col(1);
    return V_2db;
}

Eigen::Matrix3d rotationVote(const Eigen::MatrixXd& V_3d,
                             const Eigen::MatrixXd& V_2d,
                             const Eigen::MatrixXi& F,
                             const Eigen::RowVector3d& target_3d,
                             const Eigen::RowVector3d& target_2d){
                                 
    if (target_2d(2) != 0) std::cout << "Wrong usage of rotationVote" << std::endl;
    if (V_2d.cols() != 3) std::cout << "Wrong usage of rotationVote, V_2d should have 3 cols for convenience" << std::endl;

    // For each triangle
    //     proj axis on 3D plane
    //     transport to 2D
    //     add to average
    // 

    // ABC triangle in 2d
    // ApBpCp triangle in 3d
    // Ep = Ap + vector
    // then project Ep in Ep_proj, and convert to 2D


    // For the case of 90 degrees angle with the desired grain direction we have no good solution, 
    // it’s a singularity (branching point). So we need to isolate the singularity. 
    // Meaning, we compute all the individual angles in the range [-90, 90], 
    // and then we discard the angles who are in the range [90-epsilon, 90] and [-90, -90+epsilon]. 
    // We simply don’t include those in the averaging. We average the remaining ones just normally, 
    // arithmetic average, and get the rotation matrix from that.
    // When you compute the individual rotation angle, you take the smallest absolute value angle between the options alpha, 180-alpha
    // and then you take the equivalent of that angle in the -90,90 quadrant
    
    Eigen::VectorXd area;
    igl::doublearea(V_3d, F, area);
    
    double total_weights = 0.0;
    Eigen::RowVector3d average_proj = Eigen::RowVector3d::Zero();
    for (int f_id=0; f_id<F.rows(); f_id++){
        Eigen::RowVector3d A = V_2d.row(F(f_id, 0));
        Eigen::RowVector3d B = V_2d.row(F(f_id, 1));
        Eigen::RowVector3d C = V_2d.row(F(f_id, 2));
        Eigen::RowVector3d Ap = V_3d.row(F(f_id, 0));
        Eigen::RowVector3d Bp = V_3d.row(F(f_id, 1));
        Eigen::RowVector3d Cp = V_3d.row(F(f_id, 2));
        Eigen::RowVector3d Ep = Ap + target_3d;
 
        Eigen::Vector3d BpAp = (Bp - Ap).transpose();
        Eigen::Vector3d CpAp = (Cp - Ap).transpose();
        Eigen::RowVector3d n = (BpAp.cross(CpAp)).transpose();
        n = n.normalized();
        Eigen::RowVector3d v = Ep - Ap;
        double dist = v(0) * n(0) + v(1) * n(1) + v(2) * n(2);
        Eigen::RowVector3d Ep_proj = Ep - dist * n;

        Eigen::Vector3d bary_E_proj = barycentricCoords(Ep_proj, Ap, Bp, Cp);
        Eigen::RowVector3d E = bary_E_proj(0) * A
                             + bary_E_proj(1) * B
                             + bary_E_proj(2) * C;
        Eigen::RowVector3d AE = E - A;

        double weight = area(f_id);
        if (std::isnan(weight) || weight <= 10e-6){
            //std::cout << "0 area face found during rotationVote" << std::endl;
            continue;
        }
        if (AE.hasNaN()){
            std::cout << "Failed to find orientation for a face, ignoring it." << std::endl;
            continue;
        }

        total_weights += weight;
        average_proj += AE * weight;
    }

    // average isn't really what we need
    average_proj /= total_weights;
    return computeRotation(average_proj, target_2d);
}