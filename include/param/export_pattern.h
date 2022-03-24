#pragma once

#include <clipper.hpp>
#include <Eigen/Core>
#include <iostream>
#include <string>
#include <fstream>

const std::string SVG_INTRO = "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"no\"?> \
 \n\
<svg  \n\
   xmlns:dc=\"http://purl.org/dc/elements/1.1/\"  \n\
   xmlns:cc=\"http://creativecommons.org/ns#\"  \n\
   xmlns:rdf=\"http://www.w3.org/1999/02/22-rdf-syntax-ns#\"  \n\
   xmlns:svg=\"http://www.w3.org/2000/svg\"  \n\
   xmlns=\"http://www.w3.org/2000/svg\"  \n\
   xmlns:sodipodi=\"http://sodipodi.sourceforge.net/DTD/sodipodi-0.dtd\"  \n\
   xmlns:inkscape=\"http://www.inkscape.org/namespaces/inkscape\"  \n\
   id=\"svg8\"  \n\
   sodipodi:docname=\"drawing.svg\" \n";

const std::string SVG_LAYER =  "<g  \n\
     inkscape:label=\"Layer 1\" \n\
     inkscape:groupmode=\"layer\" \n\
     id=\"layer1\"> \n\
";

const std::string SVG_OUTRO = "</g>\n</svg>";

Eigen::MatrixXd createPatternLine(const Eigen::MatrixXd& V_2d, const Eigen::VectorXi& bnd, double allowance_size = 0.2){
    ClipperLib::ClipperOffset co;

    double scaling = 1000.0;
    Eigen::MatrixXd V_2db = V_2d * scaling;

    // Single patch
    ClipperLib::Path clipOutline;
    for (size_t j=0; j<bnd.rows(); j++)
        clipOutline << ClipperLib::IntPoint(V_2db(bnd(j), 0), V_2db(bnd(j), 1));
    co.AddPath(clipOutline, ClipperLib::jtRound, ClipperLib::etClosedPolygon);
    

    ClipperLib::Paths solution;
    co.Execute(solution, allowance_size * scaling);


    int n_points = 0;
    for (size_t i=0; i<solution.size(); i++){
        n_points +=solution[i].size();
    }


    Eigen::MatrixXd allowance(n_points, 2);
    int next = 0;
    for (size_t i=0; i<solution.size(); i++){
        for (size_t j=0;j<solution[i].size();j++){
            ClipperLib::IntPoint clipPos=solution[i][j];
            allowance(next, 0) =  clipPos.X;
            allowance(next, 1) =  clipPos.Y;
            next ++;
        }
    }

    allowance /= scaling;

    return allowance;
}

void createPattern(const Eigen::MatrixXd& V_2d, const Eigen::VectorXi& bnd, double allowance_size,
                   Eigen::MatrixXd& edge_begs, Eigen::MatrixXd& edge_ends){
    auto lineToEdges = [](const Eigen::MatrixXd& line, 
                          Eigen::MatrixXd& edge_begs, 
                          Eigen::MatrixXd& edge_ends){
        int n = line.rows();
        edge_begs = line;
        edge_ends.resize(n, line.cols());
        edge_ends.topRows(n - 1) = edge_begs.bottomRows(n - 1);
        edge_ends.row(n - 1) = edge_begs.row(0);
    };
    Eigen::MatrixXd outline = createPatternLine(V_2d, bnd, 0.0);
    Eigen::MatrixXd allowance = createPatternLine(V_2d, bnd, allowance_size);
    Eigen::MatrixXd edge_begs0, edge_ends0, edge_begs1, edge_ends1;
    lineToEdges(outline, edge_begs0, edge_ends0);
    lineToEdges(allowance, edge_begs1, edge_ends1);

    edge_begs.resize(edge_begs0.rows()+edge_begs1.rows(), edge_begs0.cols());
    edge_begs << edge_begs0, edge_begs1;
    edge_ends.resize(edge_ends0.rows()+edge_ends1.rows(), edge_ends0.cols());
    edge_ends << edge_ends0, edge_ends1;
}

void savePattern(const Eigen::MatrixXd& V_2d, const Eigen::VectorXi& bnd, double allowance_size){

    Eigen::MatrixXd outline = createPatternLine(-V_2d, bnd, 0.0);
    Eigen::MatrixXd allowance = createPatternLine(-V_2d, bnd, allowance_size);

    auto pointToSvgLines = [](const Eigen::MatrixXd& points){
        std::string res = "<path \n style=\"fill:none;stroke:#000000;stroke-width:1px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1\" d=\"M ";
        for (int i=0; i<points.rows(); i++){
            res += std::to_string(100.0 * points(i,0)) + "," + std::to_string(100.0 * points(i,1)) + " "; 
        }

        res += "Z\" id=\"path3717\" inkscape:connector-curvature=\"0\" />\n";
        return res;
    };

    std::string line;
    std::ofstream myfile("../example.svg");
    if (myfile.is_open()){
        myfile << SVG_INTRO;

        int mm_scale0 = static_cast<int>(100.0 *(allowance.col(0).maxCoeff() - allowance.col(0).minCoeff()));
        int mm_scale1 = static_cast<int>(100.0 *(allowance.col(1).maxCoeff() - allowance.col(1).minCoeff()));

        int start0 = -20*0 + static_cast<int>(100.0 * allowance.col(0).minCoeff());
        int start1 = -20*0 + static_cast<int>(100.0 * allowance.col(1).minCoeff());
        int end0 = 20*0 + static_cast<int>(100.0 * allowance.col(0).maxCoeff());
        int end1 = 20*0 + static_cast<int>(100.0 * allowance.col(1).maxCoeff());

        myfile << "width=\"" + std::to_string(mm_scale0) + "mm\"  \n\
                   height=\"" + std::to_string(mm_scale1) + "mm\"  \n\
                   viewBox=\""+ std::to_string(start0) +" "+ std::to_string(start1) +" "+ std::to_string(mm_scale0) +" "+ std::to_string(mm_scale1) +"\">   \n";


        myfile << SVG_LAYER;
        myfile << pointToSvgLines(outline);
        myfile << pointToSvgLines(allowance);
        myfile << SVG_OUTRO;
        myfile.close();
    }
    else std::cout << "Unable to open file" << std::endl; 
}
