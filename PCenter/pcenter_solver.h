#pragma once
#ifndef PCENTER_SOLVER_H
#define PCENTER_SOLVER_H

#include <vector>
#include <random>
#include <string>
#include <sstream>
#include <ctime>

namespace pcenter_solver {

#define NONE -1

class PCenterSolver {
public:
    const double kINF = INT32_MAX;
    const clock_t kStartClock = clock();
    const time_t kStartTime = time(NULL);
    const std::string kAlgorithm = "TabuSearch";
    const std::string kInstanceName;

    struct Configure {
        int consider_range = 50;
        int timeout_seconds = 10;
        int max_tabu_steps = INT32_MAX;
        unsigned int random_seed = 99995011 / (time(NULL) % 255) + 1;
        std::string ToString() {
            std::ostringstream res;
            res << "(" << consider_range << "/"
                << timeout_seconds << "/"
                << max_tabu_steps << "/"
                << random_seed << ")";
            return res.str();
        }
    };
public:
    explicit PCenterSolver(const char* inst_name, const int facility_num = 10):kInstanceName(inst_name){
        facility_num_ = facility_num;
        srand(cfg_.random_seed);
    }
    explicit PCenterSolver(const char* inst_name, Configure &cfg, const int facility_num = 10) : cfg_(cfg), kInstanceName(inst_name) {
        facility_num_ = facility_num;
        srand(cfg_.random_seed);
    }
    ~PCenterSolver() {}

    void Solve();
    void Record(const char *log_path = "log.csv");
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

    struct Infomation {
        clock_t time_stamp = 0;
        int iteration_stamp = 0;
        int generation_stamp = 0;
    } info;
protected:
    void LoadGraph(const std::string path);
    void GenInitSolution();
    void GetKthNeighbors(int node, int k, std::vector<int> &res);
    void FindMove(int k, int &choosed_user, int &choosed_facility);
    void MakeMove(int choosed_user, int choosed_facility);
    void TabuSearch();
    const bool Check() const;  //Check the validity of the solution

    const int RandomVertex() const { return rand() % vertex_num_; }
    const bool IsTimeOut() const { return (clock() - kStartClock) / CLOCKS_PER_SEC > cfg_.timeout_seconds; }
    const std::string FormatTime(const time_t time, const char *format = "%Y-%m-%d %H:%M:%S") {
        char buf[64];
        tm t;
        localtime_s(&t, &time);
        strftime(buf, sizeof buf, format, &t);
        return std::string(buf);
    }
    const std::string GetTimeStamp() {
        std::ostringstream res;
        if (info.time_stamp - kStartClock > CLOCKS_PER_SEC)
            res << static_cast<double>(info.time_stamp - kStartClock) / CLOCKS_PER_SEC << "s";
        else
            res << info.time_stamp - kStartClock << "ms";
        return res.str();
    }
private:
    int facility_num_;
    int vertex_num_;
    double best_objval_;  //best objective value
    double current_objval_;
    int iterator_num_;
    int base_user_tabu_steps_;
    int base_facility_tabu_steps_;
    Configure cfg_;

    std::vector<std::vector<double>> graph_matrix_;  //the origin graph form instance
    std::vector<std::vector<int>> sorted_neighbors_;
    std::vector<int> user_tabu_table_;
    std::vector<int> facility_tabu_table_;
    std::vector<FDPair> FDtable_;  //facility and distance table for every node
    std::vector<int> facility_nodes_;
    std::vector<int> best_facility_nodes_;
    std::vector<bool> isfacility_;
};

}

#endif // !PCENTER_SOLVER_H
