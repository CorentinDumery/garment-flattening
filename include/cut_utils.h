/**
 * @author Corentin Dumery
 * @brief Some functions to cut a mesh given a input set of vertices.
 * @date 2022-02-16
 * 
 */

#pragma once
#include <Eigen/Core>
#include <vector>

// Example usage:
/*
Eigen::MatrixXi cuts = cutsMatrix(selection, F);
Eigen::VectorXi I; // vertex correspondence pre/post cut
igl::cut_mesh(V, F, cuts, I);
std::vector<int> ordered_cut = identifyCut(I, F, selection);
*/

Eigen::MatrixXi cutsMatrix(const std::vector<int>& selected,
                           const Eigen::MatrixXi& F){
    Eigen::MatrixXi cuts = Eigen::MatrixXi::Zero(F.rows(), 3);

    for (int i=0; i<F.rows(); i++){
        for (int j=0; j<selected.size(); j++){
            for (int k=0; k<selected.size(); k++){
                if (k==j) continue;
                for (int l=0; l<3; l++){
                    if (selected[j] == F(i,l)){
                        int m = (l+1)%3;
                        if (selected[k] == F(i,m)){
                            cuts(i, l) = 1;
                        } 
                    }
                }
            }
        }
    }
    return cuts;
}

std::vector<std::vector<int>> computeConnectivity(const std::vector<int>& selected,
                                                  const Eigen::MatrixXi& F){
    std::vector<std::vector<int>> neighbors;
    for (int i=0; i<selected.size(); i++){
        neighbors.push_back({});
    }
    for (int i=0; i<F.rows(); i++){
        for (int j=0; j<3; j++){
            for (int k=0; k<selected.size(); k++){
                if (F(i,j) == selected[k]){
                    for (int l=0; l<selected.size(); l++){
                        if (l==k) continue;
                        if (selected[l] == F(i,(j+1)%3)) neighbors[k].push_back(l);
                        if (selected[l] == F(i,(j+2)%3)) neighbors[k].push_back(l);
                    }
                }
            }    
        }
    }

    for (int i=0; i<neighbors.size(); i++){
        std::sort(neighbors[i].begin(), neighbors[i].end());
        neighbors[i].erase(std::unique(neighbors[i].begin(), neighbors[i].end()), neighbors[i].end());

    }

    return neighbors;
}

std::vector<int> restoreCutOrder(const std::vector<int>& selected,
                                 const Eigen::MatrixXi& F){

    std::vector<std::vector<int>> neigbhors = computeConnectivity(selected, F);

    int curr_v = -1;
    for (int i=0; i<neigbhors.size(); i++){
        if (neigbhors[i].size()==1){
            curr_v = i;
            break;
        }
    }

    if (curr_v < 0){
        std::cout << "ERROR: cut starting point not found" << std::endl;
        return {};
    }
    
    std::vector<int> new_selected = {selected[curr_v]};
    int next_v = neigbhors[curr_v][0];
    while (neigbhors[next_v].size() > 1){
        if (neigbhors[next_v][0] == curr_v){
            curr_v = next_v;
            next_v = neigbhors[next_v][1];
        }
        else if (neigbhors[next_v][1] == curr_v){
            curr_v = next_v;
            next_v = neigbhors[next_v][0];
        }
        else {
            std::cout << "Failed to order cutting path." << std::endl;
            break;
        }
        new_selected.push_back(selected[curr_v]);
    }
    new_selected.push_back(selected[next_v]);

    return new_selected;
}

std::vector<int> identifyCut(const Eigen::VectorXi& corres, 
                             const Eigen::MatrixXi& F,
                             const std::vector<int>& selected){ // assumes pre-cut ids are still used in cut (on either side)
    std::vector<int> selected_dupl = selected;
    for (int i=0; i<corres.size(); i++){
        if (corres[i] != i){
            for (int j=0; j<selected.size(); j++){
                if (corres[i] == selected[j]){
                    selected_dupl.push_back(i);
                }
            }
        }
    }
    std::vector<int> cut = restoreCutOrder(selected_dupl, F);
    return cut;
}