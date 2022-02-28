#pragma once
/**
 * @author Corentin Dumery 
 * @brief Dart class(es) to represent darts and compute symmetry targets.
 * @date 2022-02-28
 * 
 */

#include <Eigen/Core>
#include <vector>
#include <iostream>

/**
 * @brief Triangle dart, represented as a tip + pairs of matching vertices (unordered)
 * 
 */
class UnorderedDart {
    
public:
    const std::vector<std::pair<int, int>> pairs_;
    int tip_;

    UnorderedDart() {};
    UnorderedDart(const std::vector<std::pair<int, int>>& pairs, int tip)
        : pairs_(pairs), tip_(tip){};

    void print() const;

    /**
     * @brief For a triangle dart, compute ideal symmetry axis (best fitting line
     * to matching pairs centerpoints, constrained in tip vertex) 
     */
    Eigen::RowVector2d computeSymmetryAxis(const Eigen::MatrixXd& V_2d) const;

    /**
     * @brief Compute vertex targets positions as defined in our associated paper
     * 
     * @return vector of pairs where vec[i].first is the target of vertex pairs_[i].first
     */
    std::vector<std::pair<Eigen::RowVector2d, Eigen::RowVector2d>> computeSymmetryTargets(const Eigen::MatrixXd& V_2d) const;

private:
    /**
     * @brief Compute symmetric points
     * 
     * @param V_2d 
     * @param sym_axis see computeSymmetryAxis
     * @return vector of pairs where vec[i].first is the symmetric point of pairs_[i].first
     */
    std::vector<std::pair<Eigen::RowVector2d, Eigen::RowVector2d>> getSymmetricPoints(const Eigen::MatrixXd& V_2d,
                                                                                      const Eigen::RowVector2d& sym_axis) const;
    
    /**
     * @brief Move vertices to desired symmetry locations where
     * they are defined
     */
    void snapSymmetric(Eigen::MatrixXd& V_2d,
                       const Eigen::RowVector2d& sym_axis) const;
};