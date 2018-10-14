#include "solver.h"
#include <string>
#include <fstream>
#include <iostream>
#include <cassert>
#include <random>

using namespace std;

namespace pcenter_solver {

void Solver::LoadGraph(const char * path) {
    string strpath(path);
    if (strpath.rfind(".txt") != string::npos) {
        instance_type_ = InstanceType::TXT;
        //first line:<node_number> <edge_number> <P>; other lines:<node1> <node2> <distance>
        ifstream ifs(path);
        if (!ifs.is_open()) {
            cout << "ERROR:Can't open instance file!" << endl;
            return;
        } else {
            double dist = 0.0, max_dist = 0;
            int edge_num = 0, node1 = 0, node2 = 0;
            cout << "Load instance: " << strpath.substr(strpath.rfind("instance/") + 9) << endl;
            ifs >> vertex_num_ >> edge_num >> facility_num_;
            //resize graph and set it's default value as 0
            graph_matrix_.resize(vertex_num_);
            for (int i = 0; i < vertex_num_; ++i)
                graph_matrix_[i].resize(vertex_num_);
            for (int i = 0; i < vertex_num_; ++i)
                for (int j = 0; j < vertex_num_; ++j)
                    graph_matrix_[i][j] = kMaxWeight;
            for (int i = 0; i < vertex_num_; ++i)
                graph_matrix_[i][i] = 0;
            //save value to graph
            for (int i = 0; i < edge_num; ++i) {
                ifs >> node1 >> node2 >> dist;
                graph_matrix_[node1 - 1][node2 - 1] = dist;  //node id in the instance start from 1
                graph_matrix_[node2 - 1][node1 - 1] = dist;
                if (max_dist < dist)max_dist = dist;
            }
            assert(kMaxWeight > max_dist);
            ifs.close();
        }
    } else if (strpath.rfind(".tsp") != string::npos) {
        instance_type_ = InstanceType::TSP;
        //first line:<node_number>; other lines:<node> <position_x> <position_y>
        ifstream ifs(path);
        if (!ifs.is_open()) {
            cout << "ERROR:Can't open instance file!" << endl;
            return;
        } else {
            int node = 0;
            double pos_x = 0.0, pos_y = 0.0;
            cout << "Load instance: " << strpath.substr(strpath.rfind("instance/") + 9) << endl;
            ifs >> vertex_num_;
            graph_posmap_.resize(vertex_num_);
            for (int i = 0; i < vertex_num_; ++i) {
                ifs >> node >> pos_x >> pos_y;
                graph_posmap_[node - 1].push_back(pos_x);
                graph_posmap_[node - 1].push_back(pos_y);
            }
            ifs.close();
        }
    }
}

int Solver::Solve() {
    return 0;
}

std::vector<int> Solver::GetResult() {
    return std::vector<int>();
}

//Initialization shortest_dist_,shortest_path_,tabu_table,fd_table,choosed_facility
void Solver::Init() {
    //Init shortest_dist_, shortest_path_ and degrees_
    shortest_dist_.resize(vertex_num_);
    shortest_path_.resize(vertex_num_);
    degrees_.resize(vertex_num_);
    for (int i = 0; i < vertex_num_; ++i) {
        shortest_dist_[i].resize(vertex_num_);
        shortest_path_[i].resize(vertex_num_);
        degrees_[i] = 0;
    }
    if (instance_type_ == InstanceType::TXT) {
        shortest_dist_ = graph_matrix_;
        for (int i = 0; i < vertex_num_; ++i)
            for (int j = 0; j < vertex_num_; ++j)
                if (graph_matrix_[i][j] == kMaxWeight) {
                    shortest_path_[i][j] = PathType::NoPath;
                }
                else {
                    shortest_path_[i][i] = PathType::DirectTo;
                    degrees_[i]++;
                }
    } else if (instance_type_ == InstanceType::TSP) {  //get distance by position
        for (int i = 0; i < vertex_num_; ++i) {
            for (int j = i; j < vertex_num_; ++j) {
                shortest_dist_[i][j] = shortest_dist_[j][i] = sqrt(
                    pow(graph_posmap_[i][0] - graph_posmap_[j][0], 2) + pow(graph_posmap_[i][1] - graph_posmap_[j][1], 2));
            }
        }
        for (int i = 0; i < vertex_num_; ++i)
            for (int j = 0; j < vertex_num_; ++j) {
                shortest_path_[i][j] = PathType::DirectTo;
                degrees_[i]++;
            }
    }
    //floyd-warshall
    for (int k = 0; k < vertex_num_; ++k)
        for (int i = 0; i < vertex_num_; ++i)
            for (int j = 0; j < vertex_num_; ++j)
                if (shortest_dist_[i][j] < shortest_dist_[i][k] + shortest_dist_[k][j]) {
                    shortest_dist_[i][j] = shortest_dist_[i][k] + shortest_dist_[k][j];
                    shortest_path_[i][j] = k;  //record the nearest node to j in the shorest path from i to j
                }
    //Init tabu table
    tabu_table_.resize(vertex_num_);
    for (int i = 0; i < vertex_num_; ++i)
        tabu_table_[i].resize(vertex_num_);
    for (int i = 0; i < vertex_num_; ++i)
        for (int j = 0; j < vertex_num_; ++j)
            tabu_table_[i][j] = 0;
    //Init FD_table, choosed_facility
    fd_table_.resize(vertex_num_);
    choosed_facility_.resize(facility_num_);
}

void Solver::GenInitSolution() {
    int choosed_facility_num = 0;
    //choose the nodes who has no neighbor as a center
    for (int n = 0; n < vertex_num_; ++n)
        if (degrees_[n] == 0)
            choosed_facility_[choosed_facility_num++] = n;
    //choose the first center who has max degrees_
    if (choosed_facility_num < facility_num_)
        choosed_facility_[choosed_facility_num++] = 
        *max_element(degrees_.begin(), degrees_.end());
    while (choosed_facility_num < facility_num_) {
        
    }
}

bool Solver::Check() {
    return false;
}

}