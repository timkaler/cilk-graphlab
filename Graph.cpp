// Copyright Tim Kaler 2013

#include <utility>
#include <algorithm>
#include <vector>
#include <set>
#include "./Graph.h"

template<typename VertexType, typename EdgeType>
Graph<VertexType, EdgeType>::Graph<VertexType, EdgeType>() {
  nextEdgeId = 0;
}

template<typename VertexType, typename EdgeType>
void Graph<VertexType, EdgeType>::addEdge(
    int vid1, int vid2, EdgeType edgeInfo) {

  // Grow the edgeData array. This should probably be replaced with a vector.
  if (nextEdgeId%10000 == 0) {
    edgeData = reinterpret_cast<EdgeType*>(
        realloc(edgeData, (nextEdgeId+10000)*sizeof(EdgeType)));
  }

  edgeData[nextEdgeId] = edgeInfo;

  struct edge_info edge;
  edge.edge_id = nextEdgeId;
  edge.out_vertex = vid1;
  edge.in_vertex = vid2;
  temp_edges.push_back(edge);
  nextEdgeId++;
}

template<typename VertexType, typename EdgeType>
void Graph<VertexType, EdgeType>::addVertex(int vid, VertexType vdata) {
  assert(vertexCount > vid);
  vertexData[vid] = vdata;
}

template<typename VertexType, typename EdgeType>
void Graph<VertexType, EdgeType>::resize(int size) {
  vertexCount = size;
  vertexData = reinterpret_cast<VertexType*>(
      realloc(vertexData, vertexCount * sizeof(VertexType)));
}

template<typename VertexType, typename EdgeType>
void Graph<VertexType, EdgeType>::finalize() {
  out_edges = reinterpret_cast<struct edge_info*>(
      calloc(nextEdgeId, sizeof(struct edge_info)));
  in_edges = reinterpret_cast<struct edge_info*>(
      calloc(nextEdgeId, sizeof(struct edge_info)));

  out_edge_index = reinterpret_cast<int*>(calloc(nextEdgeId, sizeof(int)));
  in_edge_index = reinterpret_cast<int*>(calloc(nextEdgeId, sizeof(int)));

  outDegree = reinterpret_cast<int*>(calloc(vertexCount, sizeof(int)));
  inDegree = reinterpret_cast<int*>(calloc(vertexCount, sizeof(int)));

  // compute the in degree and out degrees.
  for (int i = 0; i < temp_edges.size(); i++) {
    outDegree[temp_edges[i].out_vertex] += 1;
    inDegree[temp_edges[i].in_vertex] += 1;
  }

  // lets compute the out_edge and in_edge indices.
  out_edge_index[0] = 0;
  in_edge_index[0] = 0;
  for (int i = 1; i < vertexCount; i++) {
    out_edge_index[i] = out_edge_index[i-1] + outDegree[i-1];
    in_edge_index[i] = in_edge_index[i-1] + inDegree[i-1];
  }

  // now put all the edges into position.
  for (int i = 0; i < temp_edges.size(); i++) {
    int out_index = out_edge_index[temp_edges[i].out_vertex]++;
    int in_index = in_edge_index[temp_edges[i].in_vertex]++;
    out_edges[out_index] = temp_edges[i];
    in_edges[in_index] = temp_edges[i];
  }

  // now we need to subtract the degrees from the indicies.
  for (int i = 0; i < vertexCount; i++) {
    out_edge_index[i] = out_edge_index[i] - outDegree[i];
    in_edge_index[i] = in_edge_index[i] - inDegree[i];
  }
}

template<typename VertexType, typename EdgeType>
void Graph<VertexType, EdgeType>::prefetch_vertex(int vid) {
  _mm_prefetch(reinterpret_cast<char*>(&vertexData[vid]), 3);
  _mm_prefetch(reinterpret_cast<char*>(&inDegree[vid]), 3);
  _mm_prefetch(reinterpret_cast<char*>(&outDegree[vid]), 3);

  _mm_prefetch(reinterpret_cast<char*>(&out_edges[vid]), 3);
  _mm_prefetch(reinterpret_cast<char*>(&in_edges[vid]), 3);
}

template<typename VertexType, typename EdgeType>
struct edge_info* Graph<VertexType, EdgeType>::getOutEdges(int vid) {
  return &out_edges[out_edge_index[vid]];
}

template<typename VertexType, typename EdgeType>
struct edge_info* Graph<VertexType, EdgeType>::getInEdges(int vid) {
  return &in_edges[in_edge_index[vid]];
}

template<typename VertexType, typename EdgeType>
VertexType* Graph<VertexType, EdgeType>::getVertexData(int vid) {
  return &vertexData[vid];
}

template<typename VertexType, typename EdgeType>
EdgeType* Graph<VertexType, EdgeType>::getEdgeData(int eid) {
  return &edgeData[eid];
}

template<typename VertexType, typename EdgeType>
int Graph<VertexType, EdgeType>::getVertexColor(int vid) {
  return vertexColors[vid];
}

template<typename VertexType, typename EdgeType>
int Graph<VertexType, EdgeType>::getOutDegree(int vid) {
  return outDegree[vid];
}

template<typename VertexType, typename EdgeType>
int Graph<VertexType, EdgeType>::getInDegree(int vid) {
  return inDegree[vid];
}

template<typename VertexType, typename EdgeType>
int Graph<VertexType, EdgeType>::num_vertices() {
  return vertexCount;
}

template<typename VertexType, typename EdgeType>
int Graph<VertexType, EdgeType>::num_edges() {
  return nextEdgeId;
}

template<typename VertexType, typename EdgeType>
void Graph<VertexType, EdgeType>::validate_coloring() {
  for (int v = 0; v < vertexCount; v++) {
    struct edge_info* inEdges = getInEdges(v);
    struct edge_info* outEdges = getOutEdges(v);

    for (int i = 0; i < inDegree[v]; i++) {
      if (getVertexColor(inEdges[i].out_vertex) == getVertexColor(v)) {
        printf("Invalid coloring \n");
        return;
      }
    }
    for (int i = 0; i < outDegree[v]; i++) {
      if (getVertexColor(outEdges[i].in_vertex) == getVertexColor(v)) {
        printf("Invalid coloring \n");
        return;
      }
    }
  }
  printf("Valid coloring \n");
}

// returns true if this vertex is in the first root set.
template<typename VertexType, typename EdgeType>
void Graph<VertexType, EdgeType>::partition(int v, int* order,
    int* partitionIndexIn, int* partitionIndexOut) {
  struct edge_info* inEdges = getInEdges(v);
  struct edge_info* outEdges = getOutEdges(v);

  int j = 0;
  for (int i = 0; i < inDegree[v]; i++) {
    if (order[inEdges[i].out_vertex] < order[v]) {
      struct edge_info tmp = inEdges[j];
      inEdges[j] = inEdges[i];
      inEdges[i] = tmp;
      j++;
    }
  }
  partitionIndexIn[v] = j;

  j = 0;
  for (int i = 0; i < outDegree[v]; i++) {
    if (order[outEdges[i].in_vertex] < order[v]) {
      struct edge_info tmp = outEdges[j];
      outEdges[j] = outEdges[i];;
      outEdges[i] = tmp;
      j++;
    }
  }
  partitionIndexOut[v] = j;
}


template<typename VertexType, typename EdgeType>
void Graph<VertexType, EdgeType>::colorVertex(int v) {
  std::set<int> neighbor_colors;
  struct edge_info * inEdges = getInEdges(v);
  struct edge_info * outEdges = getOutEdges(v);
  for (int i = 0; i < inDegree[v]; i++) {
    int u = inEdges[i].out_vertex;
    neighbor_colors.insert(getVertexColor(u));
  }
  for (int i = 0; i < outDegree[v]; i++) {
    int u = outEdges[i].in_vertex;
    neighbor_colors.insert(getVertexColor(u));
  }
  int color = 0;
  while (neighbor_colors.find(color) != neighbor_colors.end()) {
    color++;
  }
  vertexColors[v] = color;
}

template<typename VertexType, typename EdgeType>
void Graph<VertexType, EdgeType>::asyncColor(int v, int* order, int* counters,
    cilk::holder< std::set<int> >* neighbor_set_holder) {
  std::set<int>& neighbor_colors = (*neighbor_set_holder)();
  neighbor_colors.clear();

  struct edge_info * inEdges = getInEdges(v);
  struct edge_info * outEdges = getOutEdges(v);

  int* neighborColors = reinterpret_cast<int*>(
      calloc(sizeof(int), inDegree[v] + outDegree[v]));

  cilk_for (int i = 0; i < inDegree[v] + outDegree[v]; i++) {
    int u;
    if (i < inDegree[v]) {
      u = inEdges[i].out_vertex;
    } else {
      u = outEdges[i - inDegree[v]].in_vertex;
    }
    int color = vertexColors[u];
    if (order[u] < order[v] && color < inDegree[v] + outDegree[v]) {
      neighborColors[color] = 1;
    }
  }

  int color = -1;
  for (int i = 0; i < inDegree[v] + outDegree[v]; i++) {
    if (neighborColors[i] == 0) {
      color = i;
      break;
    }
  }
  if (color == -1) {
    color = inDegree[v] + outDegree[v];
  }
  vertexColors[v] = color;

  // decrement all bigger neighbors
  cilk_for (int i = 0; i < inDegree[v] + outDegree[v]; i++) {
    int u;
    if (i < inDegree[v]) {
      u = inEdges[i].out_vertex;
    } else {
      u = outEdges[i - inDegree[v]].in_vertex;
    }
    if (order[u] > order[v]) {
      __sync_fetch_and_sub(&counters[u], 1);
      if (__sync_bool_compare_and_swap(&counters[u], 0, -1)) {
        asyncColor(u, order, counters, neighbor_set_holder);
      }
    }
  }
}

template<typename VertexType, typename EdgeType>
int Graph<VertexType, EdgeType>::compute_coloring() {
  return compute_coloring_atomiccounter();
}

// Uses the prefix based method to compute a valid coloring.
template<typename VertexType, typename EdgeType>
int Graph<VertexType, EdgeType>::compute_coloring_prefix() {
  // perform a parallel coloring of the graph.
  std::vector<std::pair<int, int> > permutation(vertexCount);

  cilk_for (int v = 0; v < vertexCount; ++v) {
    permutation[v] = std::make_pair(-(inDegree[v] + outDegree[v]), v);
  }
  std::sort(permutation.begin(), permutation.end());

  int* order = reinterpret_cast<int*>(malloc(sizeof(int) * vertexCount));
  vertexColors = reinterpret_cast<int*>(malloc(sizeof(int) * vertexCount));

  cilk_for (int i = 0; i < vertexCount; i++) {
    order[permutation[i].second] = i;
    vertexColors[i] = -1;
  }

  int max_success = 0;
  int prefix_length = 256;
  int num_iterations = 0;
  cilk::reducer_max<int> max_color(-1);
  while (true) {
    bool done = true;
    int stop;
    if (max_success + prefix_length > vertexCount) {
      stop = vertexCount;
    } else {
      stop = max_success + prefix_length;
    }

    cilk::reducer_min<int> min(max_success + prefix_length);
    cilk::holder< std::set<int> > neighbor_colors_holder;
    cilk_for (int v = max_success; v < stop; v++) {
      int vid = permutation[v].second;
      if (vertexColors[vid] == -1) {
        min.calc_min(v);
        bool skip = false;
        std::set<int>& neighbor_colors = neighbor_colors_holder();
        neighbor_colors.clear();
        struct edge_info* in_edges = getInEdges(vid);
        struct edge_info* out_edges = getOutEdges(vid);

        for (int i = 0; i < inDegree[vid]; i++) {
          if (order[in_edges[i].out_vertex] < v &&
              vertexColors[in_edges[i].out_vertex] == -1) {
            skip = true;
            break;
          }
          neighbor_colors.insert(vertexColors[in_edges[i].out_vertex]);
        }

        if (!skip) {
          for (int i = 0; i < outDegree[vid]; i++) {
            if (order[out_edges[i].in_vertex] < v &&
                vertexColors[out_edges[i].in_vertex] == -1) {
              skip = true;
              break;
            }
            neighbor_colors.insert(vertexColors[out_edges[i].in_vertex]);
          }
        }

        if (skip) {
          done = false;
        } else {
          int chosen_color = -1;
          for (int j = 0; j < vertexCount; j++) {
            if (neighbor_colors.find(j) == neighbor_colors.end()) {
              chosen_color = j;
              break;
            }
          }
          max_color.calc_max(chosen_color);
          vertexColors[vid] = chosen_color;
        }
      }
    }
    max_success = min.get_value();
    if (done && max_success > vertexCount) {
      break;
    }
    num_iterations++;
  }

  int the_max_color = max_color.get_value();
  return the_max_color+1;
}

// Compute coloring using the asynchronous rootset method /w atomic counters.
template<typename VertexType, typename EdgeType>
int Graph<VertexType, EdgeType>::compute_coloring_atomiccounter() {
  // Step 1: Give the vertices an order.
  int* counters = reinterpret_cast<int*>(malloc(sizeof(int)*vertexCount));
  std::vector<std::pair<int, int> > permutation(vertexCount);
  cilk_for (int v = 0; v < vertexCount; v++) {
    permutation[v] = std::make_pair(-(inDegree[v] + outDegree[v]), v);
  }
  std::sort(permutation.begin(), permutation.end());

  int* order = reinterpret_cast<int*>(malloc(sizeof(int) * vertexCount));
  vertexColors = reinterpret_cast<int*>(malloc(sizeof(int) * vertexCount));

  cilk_for (int i = 0; i < vertexCount; i++) {
    order[permutation[i].second] = i;
    vertexColors[i] = -1;
    counters[i] = 0;
  }

  // Step 2: Initialize the counters, given the order.
  cilk_for (int i = 0; i < vertexCount; i++) {
    struct edge_info* in_edges = getInEdges(i);
    struct edge_info* out_edges = getOutEdges(i);
    int count = 0;
    for (int j = 0; j < inDegree[i]; j++) {
      int neighbor = in_edges[j].out_vertex;
      if (order[neighbor] < order[i]) {
        count++;
      }
    }
    for (int j = 0; j < outDegree[i]; j++) {
      int neighbor = out_edges[j].in_vertex;
      if (order[neighbor] < order[i]) {
        count++;
      }
    }
    counters[i] = count;
  }
  cilk::holder< std::set<int> > neighbor_set_holder;
  // Step 3: Start the computation
  cilk_for (int v = 0; v < vertexCount; v++) {
    if (counters[v] == 0) {
      asyncColor(v, order, counters, &neighbor_set_holder);
    }
  }
  // compute the maximum color.
  int maxColor = -1;
  for (int i = 0; i < vertexCount; i++) {
    maxColor = vertexColors[i] > maxColor ? vertexColors[i] : maxColor;
  }
  return maxColor + 1;
}

