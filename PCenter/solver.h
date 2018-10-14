#pragma once
#ifndef PCENTER_SOLVER_H
#define PCENTER_SOLVER_H

#include <vector>
#include <ctime>
#include <random>

namespace pcenter_solver {

class Solver {
public:
    enum InstanceType { TSP, TXT };
    const double kINF = 1000000.0;
public:
    explicit Solver(int P = 0, unsigned int seed = 0) :facility_num_(P), vertex_num_(0) {
        if (seed == 0) {  //set random seed
            unsigned int time_now = time(NULL);
            random_seed_ = time_now / ((time_now << 24) >> 24);
        } else {
            random_seed_ = seed;
        }
        srand(random_seed_);
    }
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
    void GenInitSolution();
    int RandomVertex() { return rand() % vertex_num_; }
    bool Check();  //Check the validity of the solution
private:
    int facility_num_;
    int vertex_num_;
    int instance_type_;
    unsigned int random_seed_;
    double best_objval_;  //best objective value
    std::vector<std::vector<double>> graph_matrix_;  //the origin graph form instance
    std::vector<std::vector<int>> tabu_table_;
    std::vector<FD> fd_table_;  //facility and distance table for every node
    std::vector<int> choosed_facility_;
};

}

#endif // !PCENTER_SOLVER_H
