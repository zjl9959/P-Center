#pragma once
#ifndef PCENTER_GRAPH_H
#define PCENTER_GRAPH_H

#include <vector>
#include <sstream>
#include <map>
#include <assert.h>

class Graph {

enum GraphType {
    UndirectedUnweighted = 0,
    UndirectedWeighted,
    DirectedUnweighted,
    DirectedWeighted,
};

template <typename T = int>
class AdjMatrixImpl {
public:
    AdjMatrixImpl(){}
    ~AdjMatrixImpl(){}
    /*istringstream format:
      [H]"\d \d\n"  --vertex_nums, edge_nums
      [E]"\d \d \d?\n"  --node1, node2, weight
    *///[H]:one line,[E]:one or more lines 
    void LoadGraph(std::istringstream iss, int graph_type = GraphType::UndirectedUnweighted) {
        iss >> vertex_num_ >> edge_num_;
        //resize adjacency matrix and init it
        adjmatrix_.resize(vertex_num_);
        for (int i = 0; i < vertex_num_; ++i)
            adjmatrix_[i].resize(vertex_num_);
        for (int i = 0; i < vertex_num_; ++i)
            for (int j = 0; j < vertex_num_; ++j)
                adjmatrix_[i][j] = 0;  //default weight is zero
        //load edges from stringstream
        int node1, node2; T weight;
        int vertex_count = 0;  //used for vertex_map_
        for (int i = 0; i < edge_num_; ++i) {
            iss >> node1 >> node2;
            //check has/insert mapping to vertex_map_
            auto map_node1 = vertex_map_.find(node1), map_node2 = vertex_map_.find(node2);
            if (map_node1 == vertex_map_.end())
                vertex_map_.insert(node1, vertex_count++);
            if (map_node2 == vertex_map_.end())
                vertex_map_.insert(node2, vertex_count++);
            //add node to adj_matrix_
            if (graph_type == GraphType::DirectedWeighted || graph_type == GraphType::UndirectedWeighted) {
                iss >> weight;
                adjmatrix_[vertex_map_[node1]][vertex_map_[node2]] = weight;
                //undirected graph need to expand on its edge
                if (graph_type == GraphType::UndirectedWeighted)
                    adjmatrix_[vertex_map_[node2]][vertex_map_[node1]] = weight;
            } else {
                adjmatrix_[vertex_map_[node1]][vertex_map_[node2]] = 1;
                if (graph_type == GraphType::UndirectedUnweighted)
                    adjmatrix_[vertex_map_[node2]][vertex_map_[node1]] = 1;
            }
        }
        assert(vertex_count == vertex_num_);
    }
    //call this function when save results
    int mapping_vertex(int vertex) { return vertex_map_[vertex]; }
    int vertex_num() { return vertex_num_ };
    int edge_num() { return edge_num_ };

protected:

private:
    std::vector<std::vector<T>> adjmatrix_;
    std::map<int, int> vertex_map_;
    int vertex_num_;
    int edge_num_;
};

template <typename T = int>
class AdjListImpl {
public:
    AdjListImpl(){}
    ~AdjListImpl(){}
protected:

private:

};

};

#endif // !PCENTER_GRAPH_H

