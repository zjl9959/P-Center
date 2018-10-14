#pragma once
#ifndef PCENTER_SOLVER_H
#define PCENTER_SOLVER_H

#include <vector>
#include <map>

namespace pcenter_solver {

class Solver {
public:
    enum InstanceType { TSP, TXT };
    enum PathType {DirectTo = -1, NoPath = -2};
    const double kMaxWeight = 1000000.0;
public:
    explicit Solver(int P = 0) :facility_num_(P), vertex_num_(0) {}
    ~Solver() {}
    void LoadGraph(const char *path);
    int Solve();
    std::vector<int> GetResult();
protected:
    struct FD {
        int f1;  //Nearest facility identity
        double d1;  //Nearest facility distance
        int f2;  //Second nearest facility identity
        double d2;  //Second nearest facility distance
        FD(int facility1, double distance1, int facility2, double distance2) :
            f1(facility1), d1(distance1), f2(facility2), d2(distance2) {}
    };
protected:
    void Init();  //Initialization dist_table_,tabu_table,fd_table,choosed_facility
    void GenInitSolution();
    bool Check();  //Check the validity of the solution
private:
    int facility_num_;
    int vertex_num_;
    int instance_type_;
    double best_objval_;  //best objective value
    std::vector<std::vector<double>> graph_matrix_;  //the origin graph form instance
    std::vector<std::vector<double>> graph_posmap_;  //graph expressed as position map
    std::vector<std::vector<double>> shortest_dist_;  //the shortest distance table for all the nodes
    std::vector<std::vector<int>> shortest_path_;  //the shortest path table for all the nodes, i->...->shortest_path_[i][j]->j
    std::vector<std::vector<int>> tabu_table_;
    std::vector<int> degrees_;  //the neighbors of a node
    std::vector<FD> fd_table_;  //facility and distance table for every node
    std::vector<int> choosed_facility_;
};

}

#endif // !PCENTER_SOLVER_H
