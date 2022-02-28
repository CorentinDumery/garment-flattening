#include "param/dart.h"

void UnorderedDart::print() const {
    std::cout << "Dart:" << std::endl;
    std::cout << "\t tip: " << tip_ << std::endl;
    std::cout << "\t pairs: " << std::endl;
    for (int i=0; i<pairs_.size(); i++){
        std::cout << "\t\t" << pairs_[i].first << " " << pairs_[i].second << std::endl;
    }
}

Eigen::RowVector2d UnorderedDart::computeSymmetryAxis(const Eigen::MatrixXd& V_2d) const {
    if (pairs_.size() < 1) std::cout << "ERROR: dart too short to compute sym axis" << std::endl;
    if (V_2d.cols() != 2) std::cout << "ERROR: expected 2D mesh in computeSymmetryAxis" << std::endl;
    Eigen::RowVector2d semi_point = V_2d.row(tip_);

    // Simple linear regression without the intercept term:
    // find optimal dart symmetry axis, using regression on middle points
    // between duplicated vertices, but constrained to go through
    // dart end
    // https://en.wikipedia.org/wiki/Simple_linear_regression
    Eigen::MatrixXd points(pairs_.size(), 2);
    for (int i=0; i<pairs_.size(); i++){
        points.row(i) = (V_2d.row(pairs_[i].first) + V_2d.row(pairs_[i].second))/2.0;
    }

    points = points.rowwise() - V_2d.row(tip_);

    Eigen::RowVectorXd product_xy(points.rows());
    Eigen::RowVectorXd product_xx(points.rows());
    for (int i=0; i<points.rows(); i++) {
        product_xy(i) = points(i,0) * points(i,1); 
        product_xx(i) = points(i,0) * points(i,0); 
    }

    double beta = product_xy.mean() / product_xx.mean();
    Eigen::RowVector2d sym_2d = Eigen::RowVector2d(-1, -beta);
    Eigen::RowVector2d symmetry_axis;
    symmetry_axis(0) = sym_2d(0);
    symmetry_axis(1) = sym_2d(1);

    // orient towards dart
    Eigen::RowVector2d symm_end = (V_2d.row(pairs_[0].first) + V_2d.row(pairs_[0].second))/2.0;
    Eigen::RowVector2d out_axis = symm_end - semi_point;

    if (sym_2d.dot(out_axis) < 0) symmetry_axis = -symmetry_axis;

    symmetry_axis = symmetry_axis * out_axis.norm() / sym_2d.norm();
    return symmetry_axis;
}

std::vector<std::pair<Eigen::RowVector2d, Eigen::RowVector2d>> UnorderedDart::getSymmetricPoints(const Eigen::MatrixXd& V_2d,
                                                                                    const Eigen::RowVector2d& sym_axis) const {
    std::vector<std::pair<Eigen::RowVector2d, Eigen::RowVector2d>> res;
    res.resize(pairs_.size());
    for (int i=0; i<pairs_.size(); i++){
        Eigen::RowVector2d A, B, C, AD, CE;
        A = V_2d.row(tip_);
        B = A + sym_axis;
        C = V_2d.row(pairs_[i].second);
        AD = (B - A) * (C - A).dot(B - A) / std::pow((B - A).norm(), 2);
        CE = 2 * (AD - (C - A));
        res[i].first = C + CE;

        A = V_2d.row(tip_);
        B = A + sym_axis;
        C = V_2d.row(pairs_[i].first);
        AD = (B - A) * (C - A).dot(B - A) / std::pow((B - A).norm(), 2);
        CE = 2 * (AD - (C - A));
        res[i].second = C + CE;
    }

    return res;
};

std::vector<std::pair<Eigen::RowVector2d, Eigen::RowVector2d>> UnorderedDart::computeSymmetryTargets(const Eigen::MatrixXd& V_2d) const {
    Eigen::RowVector2d sym_axis = computeSymmetryAxis(V_2d);
    std::vector<std::pair<Eigen::RowVector2d, Eigen::RowVector2d>> sym_p = getSymmetricPoints(V_2d, sym_axis);
    std::vector<std::pair<Eigen::RowVector2d, Eigen::RowVector2d>> new_pos;
    new_pos.resize(pairs_.size());
    for (int i=0; i<pairs_.size(); i++){
        new_pos[i].first = (V_2d.row(pairs_[i].first) + sym_p[i].first)/2.0;
        new_pos[i].second = (V_2d.row(pairs_[i].second) + sym_p[i].second)/2.0;
    }
    return new_pos;
}

void UnorderedDart::snapSymmetric(Eigen::MatrixXd& V_2d,
                                  const Eigen::RowVector2d& sym_axis) const {

    if (V_2d.cols() > 2) std::cout << "Warning, replacing 3D rows with 2D" << std::endl;

    std::vector<std::pair<Eigen::RowVector2d, Eigen::RowVector2d>> new_pos = computeSymmetryTargets(V_2d);
    for (int i=0; i<pairs_.size(); i++){
        V_2d.row(pairs_[i].first) = new_pos[i].first;
        V_2d.row(pairs_[i].second) = new_pos[i].second;
    }
}
