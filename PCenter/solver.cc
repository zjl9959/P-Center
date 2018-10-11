#include "solver.h"
#include <string>
#include <fstream>
#include <iostream>

using namespace std;

namespace pcenter_solver {

void Solver::LoadGraph(char * path) {
    string strpath(path);
    if (strpath.rfind(".txt") != string::npos) {
        instance_type_ = InstanceType::TXT;
        //first line:<node_number> <edge_number>; other lines:<node1> <node2> <distance>
        ifstream ifs(path);
        if (!ifs.is_open()) {
            cout << "ERROR:Can't open instance file!" << endl;
            return;
        } else {
            double dist = 0.0;
            int edge_num = 0, node1 = 0, node2 = 0;
            ifs >> vertex_num_ >> edge_num;
            //resize graph and set it's default value as 0
            graph_.resize(vertex_num_);
            for (int i = 0; i < vertex_num_; ++i)
                graph_[i].resize(vertex_num_);
            for (int i = 0; i < vertex_num_; ++i)
                for (int j = 0; j < vertex_num_; ++j)
                    graph_[i][j] = 0;
            //save value to graph
            for (int i = 0; i < edge_num; ++i) {
                ifs >> node1 >> node2 >> dist;
                graph_[node1 - 1][node2 - 1] = dist;  //node id in the instance start from 1
                graph_[node2 - 1][node1 - 1] = dist;
            }
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

void Solver::Init() {}

bool Solver::Check() {
    return false;
}

}