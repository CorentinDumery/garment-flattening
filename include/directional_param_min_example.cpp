#include <Eigen/Core>
#include <Eigen/IterativeLinearSolvers> 

Eigen::MatrixXd localGlobal(const Eigen::MatrixXd& V_2d, const Eigen::MatrixXd& V_3d, 
                            const Eigen::MatrixXi& F){

    Eigen::SparseMatrix<double> A;
    Eigen::DiagonalMatrix<double, Eigen::Dynamic> W, Wt;
    Eigen::VectorXd b, x;

    makeSparseMatrix(V_2d, V_3d, F, A, b, W); // Fills A, b, and W. Examples shown below

    Wt = W; // weights are diagonal

    Eigen::SparseMatrix<double> At = A;
    At = At.transpose(); // I get a bit paranoid with Eigen's transpose being in place or not
    
    Eigen::SparseMatrix<double> Ap = At * Wt * W * A;
    Eigen::VectorXd bp = At * Wt * W * b;

    //Ap = W * A; // simpler and faster formulation but doesn't work
    //bp = W * b;

    x = vertices2dToVector(V_2d); // Initial solution
    
    Eigen::LeastSquaresConjugateGradient<Eigen::SparseMatrix<double>> solver(Ap); 
    x = solver.solve(bp);

    Eigen::MatrixXd res = vectorToVertices2d(x);

    return res;
}

// NOTES
// 1. During "solver.solve(bp)", we often get: solver.info() == Eigen::NoConvergence
// 2. The simpler formulation "Ap = W * A" isn't working and I'm not sure why
// 3. The current weighted version is about 10x slower to solve than the initial Ax=b


// EXAMPLE makeSparseMatrix
// Given the following input consisting of two triangles, where one is distorted:
/*
V_2d.resize(4, 3);
V_3d.resize(4, 3);
F.resize(2, 3);

V_2d <<  0,   0, 0,
        1.0,   0, 0,
            0, 1.0, 0,
        1.0, 1.5, 0;

V_3d <<  0,   0, 0,
        1.0,   0, 0,
            0, 1.0, 0,
        1.0, 1.0, 0;

F << 0, 1, 2,
        1, 3, 2;

F0 = F;

Then makeSparseMatrix produces the following matrices:

W's diagonal coefficients: 
  1 // stretch u
  1 // stretch v
0.8 // angle u
0.8 // angle v
  1 // repeat for second triangle
  1
0.8
0.8

A
-1 0 1 0 0 0 0 0 // at most three non-zero values per row (barycentric coordinates)
0 -1 0 0 0 1 0 0 
-2 0 1 0 1 0 0 0 
0 -2 0 1 0 1 0 0 
0 0 0.333334 0 -1 0 0.666667 0 
0 0 0 -0.666666 0 0 0 0.666667 
0 0 -0.333333 0 -1 0 1.33333 0 
0 0 0 -0.333333 0 -1 0 1.33333 

b
       1
       1
       1
       1
 1.05409
0.666667
0.843275
0.333335
*/