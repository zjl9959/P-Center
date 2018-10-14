#include "solver.h"
#include <string>
#include <fstream>
#include <iostream>
#include <cassert>

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
                    graph_matrix_[i][j] = kINF;
            for (int i = 0; i < vertex_num_; ++i)
                graph_matrix_[i][i] = 0;
            //save value to graph
            for (int i = 0; i < edge_num; ++i) {
                ifs >> node1 >> node2 >> dist;
                graph_matrix_[node1 - 1][node2 - 1] = dist;  //node id in the instance start from 1
                graph_matrix_[node2 - 1][node1 - 1] = dist;
                if (max_dist < dist)max_dist = dist;
            }
            assert(kINF > 10 * max_dist);
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
            vector<vector<double>> graph_posmap;
            graph_posmap.resize(vertex_num_);
            for (int i = 0; i < vertex_num_; ++i) {
                ifs >> node >> pos_x >> pos_y;
                graph_posmap[node - 1].push_back(pos_x);
                graph_posmap[node - 1].push_back(pos_y);
            }
            for (int i = 0; i < vertex_num_; ++i) {
                for (int j = i; j < vertex_num_; ++j) {
                    graph_matrix_[i][j] = graph_matrix_[j][i] = sqrt(
                        pow(graph_posmap[i][0] - graph_posmap[j][0], 2) + pow(graph_posmap[i][1] - graph_posmap[j][1], 2));
                }
            }
            ifs.close();
        }
    }
}

int Solver::Solve() {
    return 0;
}

std::vector<int> Solver::GetResult() {
    return choosed_facility_;
}

void Solver::GenInitSolution() {
    int choosed_facility_num = 0;
    fd_table_.resize(vertex_num_);
    choosed_facility_.resize(facility_num_);
    //choose the first center randomly
    choosed_facility_[choosed_facility_num++] = RandomVertex();
    while (choosed_facility_num < facility_num_) {
        //[TODO]choose other centers
    }
    //Init tabu table
    tabu_table_.resize(vertex_num_);
    for (int i = 0; i < vertex_num_; ++i)
        tabu_table_[i].resize(vertex_num_);
    for (int i = 0; i < vertex_num_; ++i)
        for (int j = 0; j < vertex_num_; ++j)
            tabu_table_[i][j] = 0;
}

bool Solver::Check() {
    return false;
}

}