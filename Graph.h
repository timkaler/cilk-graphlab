#include <iostream>
#include<fstream>
#include <stdio.h>
#include <cstdlib> 
#include <vector> 
#include <list>
#include <map>
#include <string>
#include <set> 
#include <cmath>
#include <algorithm>
#include <cilk/cilk.h> 
#include <cilk/reducer_list.h>
#include <cilk/reducer_min.h>
#include <cilk/reducer_max.h>
#include <cilk/holder.h>

#ifndef GRAPH_H
#define GRAPH_H

//using namespace std; 

struct edge_info {
  int edge_id;
  int neighbor_id;
} edge_info;

template<typename VertexType, typename EdgeType>
class Graph {

  private:
    int vertexCount;
    int nextEdgeId;
    std::map<int, VertexType> vertexData;
    std::map<int, EdgeType> edgeData;
    std::map<int, int> outDegree;
    std::map<int, int> inDegree;
    
    // maps vertexId to start of out edges.
    int* out_edge_index;
    int* in_edge_index;

    struct edge_info* out_edges;
    struct edge_info* in_edges;

    //std::vector<int, edge_info> out_edges;
    //std::vector<int, edge_info> in_edges;

    std::map<int, std::vector<struct edge_info> > temp_out_edges;
    std::map<int, std::vector<struct edge_info> > temp_in_edges;

  public:
    Graph();
    int* vertexColors;
    void addEdge(int vid1, int vid2, EdgeType edgeInfo);
    void addVertex(int vid, VertexType vdata);
    void finalize();
    void resize(int size);
    int compute_coloring();
    int getVertexColor(int vid);
    int getOutDegree(int vid);
    int getInDegree(int vid);
    int num_vertices();
    int num_edges();
    struct edge_info* getOutEdges(int vid);
    struct edge_info* getInEdges(int vid);
    VertexType* getVertexData(int vid);
    EdgeType* getEdgeData(int eid);
};

#endif