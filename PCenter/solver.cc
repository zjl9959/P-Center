#include "solver.h"
#include <string>
#include <fstream>
#include <iostream>
#include <cassert>

using namespace std;

namespace pcenter_solver {

void Solver::LoadGraph(const char * path) {
    string strpath(path);
    instance_name_ = strpath.substr(strpath.rfind("instance/") + 9);
    if (strpath.rfind(".txt") != string::npos) {
        //first line:<node_number> <edge_number> <P>; other lines:<node1> <node2> <distance>
        ifstream ifs(path);
        if (!ifs.is_open()) {
            cout << "ERROR:Can't open instance file!" << endl;
            return;
        } else {
            double dist = 0.0, max_dist = 0;
            int edge_num = 0, node1 = 0, node2 = 0;
            cout << "Load instance: " << instance_name_ << endl;
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
            //floyd-warshall, construct a complete graph
            for (int k = 0; k < vertex_num_; ++k)
                for (int i = 0; i < vertex_num_; ++i)
                    for (int j = 0; j < vertex_num_; ++j)
                        if (graph_matrix_[i][j] < graph_matrix_[i][k] + graph_matrix_[k][j])
                            graph_matrix_[i][j] = graph_matrix_[i][k] + graph_matrix_[k][j];
            assert(kINF > 10 * max_dist);
            ifs.close();
        }
    } else if (strpath.rfind(".tsp") != string::npos) {
        //first line:<node_number>; other lines:<node> <position_x> <position_y>
        ifstream ifs(path);
        if (!ifs.is_open()) {
            cout << "ERROR:Can't open instance file!" << endl;
            return;
        } else {
            int node = 0;
            double pos_x = 0.0, pos_y = 0.0;
            cout << "Load instance: " << instance_name_ << endl;
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

void Solver::Solve() {
    TabuSearch();
    if (Check()) {
        cout << "--------------------------------------" << endl;
        cout << "\t instance \t iterations \t obj \t facilitys" << endl;
        cout << "\t" << instance_name_ << iterator_num_ << best_objval_;
        for (int i = 0; i < best_facility_nodes_.size(); ++i) {
            cout << best_facility_nodes_[i] << ",";
        }
        cout << "\b\a";
    }
}

std::vector<int> Solver::GetResult() {
    return best_facility_nodes_;
}

void Solver::GenInitSolution() {
    //sort the neighbors by distance for every node
    sorted_neighbors_.resize(vertex_num_);
    vector<Edge> neighbors;
    for (int i = 0; i < vertex_num_; ++i) {
        neighbors.clear();
        for (int j = 0; j < vertex_num_; ++j) {
            if (graph_matrix_[i][j] != kINF)
                neighbors.push_back(Edge(j, graph_matrix_[i][j]));
        }
        sort_heap(neighbors.begin(),neighbors.end(), [](Edge &lhs, Edge &rhs) {
            return lhs.dist < rhs.dist; });
        sorted_neighbors_.push_back(neighbors);
    }
    int choosed_facility_num = 0;
    facility_nodes_.resize(facility_num_);
    isfacility_.resize(vertex_num_);
    //choose the first center randomly
    facility_nodes_[choosed_facility_num++] = RandomVertex();
    isfacility_[facility_nodes_[0]] = true;
    while (choosed_facility_num < facility_num_) {
        //get the longest service edge
        double longest_service_dist = -kINF;
        int longest_service_node = NONE;
        for (int u = 0; u < vertex_num_; ++u) {  //u duputy user node identity
            double service_edge = kINF;
            if (!isfacility_[u]) {
                for (int f : facility_nodes_)  //for every user node, choose a facility
                    if (graph_matrix_[f][u] < service_edge)
                        service_edge = graph_matrix_[f][u];
                if (service_edge != kINF && longest_service_dist < service_edge) {
                    longest_service_dist = service_edge;
                    longest_service_node = u;
                }
            }
        }
        //random select a node form the kth nearest neighbor of longest_user_node
        vector<Edge> neighborks;
        GetKthNeighbors(longest_service_node, kConsiderRange, neighborks);
        assert(neighborks.size());  //neighborks size should larger than 0
        facility_nodes_[choosed_facility_num++] = neighborks[rand() % neighborks.size()].node;
        isfacility_[facility_nodes_[choosed_facility_num - 1]] = true;
    }
    //Init tabu table
    user_tabu_table_.resize(vertex_num_);
    facility_tabu_table_.resize(vertex_num_);
    for (int i = 0; i < vertex_num_; ++i)
        user_tabu_table_[i] = facility_tabu_table_[i] = 0;
    //Init FDTable and object value
    best_objval_ = -kINF;
    FDtable_.resize(vertex_num_);
    for (int i = 0; i < vertex_num_; ++i) {
        Edge first_facility, second_facility;
        for (int f : facility_nodes_)
            if (graph_matrix_[i][f] < first_facility.dist) {
                first_facility.node = f;
                first_facility.dist = graph_matrix_[i][f];
            } else if (graph_matrix_[i][f] < second_facility.dist) {
                second_facility.node = f;
                second_facility.node = graph_matrix_[i][f];
            }
        if (current_objval_ < first_facility.dist)
            current_objval_ = first_facility.dist;
    }
    //record the best solution
    best_objval_ = current_objval_;
    best_facility_nodes_ = facility_nodes_;
    base_user_tabu_steps_ = facility_num_ / 10;
    base_facility_tabu_steps_ = (vertex_num_ - facility_num_) / 10;
}

//GetKthNeighbor who is not facility, res should be empty
void Solver::GetKthNeighbors(int node, int k, vector<Edge>& res) {
    int index = 0;
    while (k--) {
        while (isfacility_[sorted_neighbors_[node][index++].node]) {
            if (index == vertex_num_)return;
        }
        res.push_back(sorted_neighbors_[node][index - 1]);
    }
}

void Solver::FindMove(int & choosed_user, int & choosed_facility) {
    int hail_user;  //the user node on the longest service edge
    double hail_dist = -kINF;
    int tabu_best_user = NONE, tabu_best_facility = NONE, notabu_best_user = NONE, 
        notabu_best_facility = NONE, tabu_same_count = 0, notabu_same_count = 0;
    double tabu_best_obj = kINF, notabu_best_obj = kINF;
    vector<Edge> neighborks;  //the k's nearest neighbors of the hail user
    vector<FDPair> temp_FDtable;
    //get the longest service node and the distance
    for (int i = 0; i < vertex_num_; ++i) {
        if (FDtable_[i].nearest.dist > hail_dist) {
            hail_user = i;
            hail_dist = FDtable_[i].nearest.dist;
        }
    }
    GetKthNeighbors(hail_user, kConsiderRange, neighborks);
    for (Edge e:neighborks) {
        if (e.dist < hail_dist) {  //consider the neighbor who can improve the object value
            temp_FDtable = FDtable_;  //copy the oldversion FDtable
            //update the FDtable after add a facility
            for (int i = 0; i < vertex_num_; ++i) {
                if (temp_FDtable[i].nearest.dist > graph_matrix_[e.node][i]) {
                    temp_FDtable[i].nearest = e;
                    temp_FDtable[i].second_nearest = temp_FDtable[i].nearest;
                }
                else if (temp_FDtable[i].second_nearest.dist > graph_matrix_[e.node][i]) {
                    temp_FDtable[i].second_nearest = e;
                }
            }
            //remove a facility and evluate it
            for (int f : facility_nodes_) {  //f is the facility to remove
                double longest_service_dist = -kINF;
                for (int i = 0; i < vertex_num_; ++i)
                    if (f = temp_FDtable[i].nearest.node && longest_service_dist < temp_FDtable[i].second_nearest.dist)
                        longest_service_dist = temp_FDtable[i].second_nearest.dist;
                //record best move
                if (facility_tabu_table_[e.node] < iterator_num_ || user_tabu_table_[f] < iterator_num_) {  //no tabu
                    if (longest_service_dist < tabu_best_obj) {
                        tabu_best_obj = longest_service_dist;
                        tabu_best_user = e.node;
                        tabu_best_facility = f;
                    } else if (longest_service_dist == tabu_best_obj && (rand() % tabu_same_count) == 0) {
                        //pool sampling for same movment
                        tabu_best_obj = longest_service_dist;
                        tabu_best_user = e.node;
                        tabu_best_facility = f;
                        tabu_same_count++;
                    }
                } else if (longest_service_dist < notabu_best_obj) {  //tabu
                    notabu_best_obj = longest_service_dist;
                    notabu_best_user = e.node;
                    notabu_best_facility = f;
                } else if (longest_service_dist == notabu_best_obj && (rand() % notabu_same_count) == 0) {
                    notabu_best_obj = longest_service_dist;
                    notabu_best_user = e.node;
                    notabu_best_facility = f;
                    notabu_same_count++ ;
                }
            }
        }
    }
    if (tabu_best_obj > best_objval_&&tabu_best_obj > notabu_best_obj) {  //release tabu
        current_objval_ = tabu_best_obj;
        choosed_user = tabu_best_user;
        choosed_facility = tabu_best_facility;
    } else {
        current_objval_ = notabu_best_obj;
        choosed_user = notabu_best_user;
        choosed_facility = notabu_best_facility;
    }
}

void Solver::MakeMove(int choosed_user, int choosed_facility) {
    //update facility_nodes_, is_facility_
    for(int i=0;i<facility_nodes_.size();++i)
        if (facility_nodes_[i] = choosed_facility) {
            facility_nodes_[i] = choosed_user;
        }
    isfacility_[choosed_user] = true;
    isfacility_[choosed_facility] = false;
    //update FDtable
    for (int i = 0; i < vertex_num_; ++i) {  //add choosed user as facility
        if (FDtable_[i].nearest.dist > graph_matrix_[choosed_user][i]) {
            FDtable_[i].second_nearest = FDtable_[i].nearest;
            FDtable_[i].nearest.node = choosed_user;
            FDtable_[i].nearest.dist = graph_matrix_[choosed_user][i];
        } else if (FDtable_[i].second_nearest.dist > graph_matrix_[choosed_user][i]) {
            FDtable_[i].second_nearest.node = choosed_user;
            FDtable_[i].second_nearest.dist = graph_matrix_[choosed_user][i];
        }
    }
    for (int i = 0; i < vertex_num_; ++i) {  //remove choosed facility
        if (FDtable_[i].nearest.node == choosed_facility) {
            FDtable_[i].nearest = FDtable_[i].second_nearest;
            FDtable_[i].second_nearest.dist = kINF;
            for (int f : facility_nodes_) {  //update second nearest facility
                if (graph_matrix_[f][i] >= FDtable_[i].nearest.dist&&
                    graph_matrix_[f][i] < FDtable_[i].second_nearest.dist && f != FDtable_[i].nearest.node) {
                    FDtable_[i].second_nearest.dist = graph_matrix_[f][i];
                    FDtable_[i].second_nearest.node = f;
                }
            }
        }
    }
    //update tabu tables
    user_tabu_table_[choosed_user] = iterator_num_ + base_user_tabu_steps_ + rand() % 10;
    facility_tabu_table_[choosed_facility] = iterator_num_ + base_facility_tabu_steps_ + rand() % 10;
    //update best objective
    if (current_objval_ < best_objval_) {
        best_objval_ = current_objval_;
        best_facility_nodes_ = facility_nodes_;
        cout << "A better obj: " << best_objval_ << " on iteration " << iterator_num_ << endl;
    }
}

void Solver::TabuSearch() {
    int choose_user, choose_facility;
    while (!IsTimeOut() && iterator_num_<kMaxTabuSteps) {
        FindMove(choose_user,choose_facility);
        MakeMove(choose_user,choose_facility);
        iterator_num_++;
    }
}

bool Solver::Check() {
    return false;
}

}