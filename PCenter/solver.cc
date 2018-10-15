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
    return facility_nodes_;
}

void Solver::GenInitSolution() {
    //sort the neighbors by distance for every node
    dist_table_.resize(vertex_num_);
    vector<Edge> neighbors;
    for (int i = 0; i < vertex_num_; ++i) {
        neighbors.clear();
        for (int j = 0; j < vertex_num_; ++j) {
            if (i != j && graph_matrix_[i][j] != kINF)
                neighbors.push_back(Edge(j, graph_matrix_[i][j]));
        }
        sort_heap(neighbors.begin(),neighbors.end(), [](Edge &lhs, Edge &rhs) {
            return lhs.dist < rhs.dist; });
        dist_table_.push_back(neighbors);
    }
    int choosed_facility_num = 0;
    facility_nodes_.resize(facility_num_);
    isfacility_.resize(vertex_num_);
    //choose the first center randomly
    facility_nodes_[choosed_facility_num++] = RandomVertex();
    isfacility_[facility_nodes_[0]] = true;
    while (choosed_facility_num < facility_num_) {
        //get the longest service edge
        double longest_service_edge = -kINF;
        int longest_service_node = NONE;
        for (int u = 0; u < vertex_num_; ++u) {  //u duputy user node identity
            double service_edge = kINF;
            if (!isfacility_[u]) {
                for (int f : facility_nodes_)  //for every user node, choose a facility
                    if (graph_matrix_[f][u] < service_edge)
                        service_edge = graph_matrix_[f][u];
                if (service_edge != kINF && longest_service_edge < service_edge) {
                    longest_service_edge = service_edge;
                    longest_service_node = u;
                }
            }
        }
        //random select a node form the kth nearest neighbor of longest_user_node
        vector<int> kth_neighbor;
        for (Edge n : dist_table_[longest_service_node])
            if (!isfacility_[n.node]) kth_neighbor.push_back(n.node);
        facility_nodes_[choosed_facility_num++] = kth_neighbor[rand() % kConsiderRange];
        isfacility_[facility_nodes_[choosed_facility_num - 1]] = true;
    }
    //Init tabu table
    tabu_table_.resize(vertex_num_);
    for (int i = 0; i < vertex_num_; ++i)
        tabu_table_[i].resize(vertex_num_);
    for (int i = 0; i < vertex_num_; ++i)
        for (int j = 0; j < vertex_num_; ++j)
            tabu_table_[i][j] = 0;
    //Init FDTable
    FDtable_.resize(vertex_num_);
    vector<Edge> edge_pair(2);
    for (int i = 0; i < vertex_num_; ++i) {
        edge_pair.clear();  //must be clear
        GetKthNeighbors(i, 2, edge_pair);
        if (isfacility_[i])
            FDtable_[i] = FDPair(i, 0, edge_pair[0].node, edge_pair[0].dist);
        else
            FDtable_[i] = FDPair(edge_pair[0].node, edge_pair[0].dist, edge_pair[1].node, edge_pair[1].dist);
    }
}

//GetKthNeighbor who is not facility, res should be empty
void Solver::GetKthNeighbors(int node, int k, vector<Edge>& res) {
    int index = 0;
    while (k--) {
        while (isfacility_[dist_table_[node][index++].node]) {
            if (index == vertex_num_)return;
        }
        res.push_back(dist_table_[node][index - 1]);
    }
}

bool Solver::Check() {
    return false;
}

}