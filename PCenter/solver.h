#pragma once
#ifndef PCENTER_SOLVER_H
#define PCENTER_SOLVER_H

#include <vector>
#include <ctime>
#include <random>

namespace pcenter_solver {

#define NONE -1

class Solver {
public:
    enum InstanceType { TSP, TXT };
    const double kINF = INT32_MAX;
    const int kConsiderRange = 3;
    const unsigned int kRandomSeed;
public:
    explicit Solver(int P = 0, unsigned int seed = 99995011 / (time(NULL) % 255) + 1) :
        facility_num_(P), vertex_num_(0), kRandomSeed(seed) { srand(kRandomSeed); }
    ~Solver() {}
    void LoadGraph(const char *path);
    int Solve();
    std::vector<int> GetResult();
protected:
    struct Edge {
        int node;
        double dist;
        Edge(int Node = NONE, double distance = INT32_MAX) :node(Node), dist(distance) {}
    };
    struct FDPair {  //nearest/second_nearest facility/distance table
        Edge nearest;
        Edge second_nearest;
        FDPair(int nearest_facility = NONE, double nearest_dist = INT32_MAX, int second_nearest_facility = NONE, double second_nearest_dist = INT32_MAX) :
            nearest(nearest_facility, nearest_dist), second_nearest(second_nearest_facility, second_nearest_dist) {}
    };
protected:
    void GenInitSolution();
    int RandomVertex() { return rand() % vertex_num_; }
    void GetKthNeighbors(int node, int k, std::vector<Edge> &res);
    bool Check();  //Check the validity of the solution
private:
    int facility_num_;
    int vertex_num_;
    int instance_type_;
    double best_objval_;  //best objective value
    std::vector<std::vector<double>> graph_matrix_;  //the origin graph form instance
    std::vector<std::vector<Edge>> dist_table_;
    std::vector<std::vector<int>> tabu_table_;
    std::vector<FDPair> FDtable_;  //facility and distance table for every node
    std::vector<int> facility_nodes_;
    std::vector<bool> isfacility_;
};

}

#endif // !PCENTER_SOLVER_H
