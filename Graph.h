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

struct edge_info {
  int edge_id;
  int out_vertex;
  int in_vertex;
  int next;
} edge_info;

template<typename VertexType, typename EdgeType>
class Graph {

  private:
    int vertexCount;
    int nextEdgeId;
    std::vector<struct edge_info> temp_edges;
    
    VertexType* vertexData;
    EdgeType* edgeData;
    int* outDegree;
    int* inDegree;
 
    // maps vertexId to start of out edges.
    int* out_edge_index;
    int* in_edge_index;

    struct edge_info* out_edges;
    struct edge_info* in_edges;

  public:
    Graph();
    int* vertexColors;
    int compute_trivial_coloring();
    void addEdge(int vid1, int vid2, EdgeType edgeInfo);
    void addVertex(int vid, VertexType vdata);
    void partition(int v, int* order, int* partitionIndexIn, int* partitionIndexOut);
    bool updateIndices(int r, int v, int* order,
        int* partitionIndexIn, int* partitionIndexOut, int* currentIndexIn,
        int* currentIndexOut, int* currentIndexInDynamic, int* currentIndexOutDynamic);
    void colorVertex(int v);
    void asyncColor(int v, int* order, int* counters, cilk::holder< std::set<int> >* neighbor_set_holder);
    void finalize();
    void resize(int size);
    void prefetch_vertex(int vid);
    int compute_coloring();
    int compute_coloring_rootset();
    int compute_coloring_atomiccounter();
    void validate_coloring();
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
