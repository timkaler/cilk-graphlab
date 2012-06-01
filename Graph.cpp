#include "Graph.h"
template<typename VertexType, typename EdgeType>
Graph<VertexType, EdgeType>::Graph< VertexType,  EdgeType>() {
  // do nothing.
}

template<typename VertexType, typename EdgeType>
void Graph< VertexType,  EdgeType>::addEdge(int vid1, int vid2, EdgeType edgeInfo){
  // increment vertex degrees.
  temp_outDegree[vid1]++;
  temp_inDegree[vid2]++;

  // store the edge data.
  temp_edgeData[nextEdgeId] = edgeInfo;

  // create in and out edges.
  struct edge_info out_edge;
  struct edge_info in_edge;
  out_edge.neighbor_id = vid2;
  out_edge.edge_id = nextEdgeId;
  in_edge.neighbor_id = vid1;
  in_edge.edge_id = nextEdgeId;
 
  // store edges in temporary adjacency lists. 
  temp_out_edges[vid1].push_back(out_edge);
  temp_in_edges[vid2].push_back(in_edge); 

  nextEdgeId++; 
}

template<typename VertexType, typename EdgeType>
void Graph<VertexType, EdgeType>::addVertex(int vid, VertexType vdata){
  assert(vertexCount > vid);
  //vertexCount++;
  temp_vertexData[vid] = vdata;
}

template<typename VertexType, typename EdgeType>
void Graph<VertexType, EdgeType>::resize(int size){
  vertexCount = size;
}

template<typename VertexType, typename EdgeType>
void Graph<VertexType, EdgeType>::finalize(){
  out_edges = (struct edge_info*) malloc(sizeof(struct edge_info) * nextEdgeId);
  in_edges = (struct edge_info*) malloc(sizeof(struct edge_info) * nextEdgeId);

  out_edge_index = (int*) malloc(sizeof(int) * nextEdgeId);
  in_edge_index = (int*) malloc(sizeof(int) * nextEdgeId);

  vertexData = (VertexType*) malloc(sizeof(VertexType) * vertexCount); 
  edgeData = (EdgeType*) malloc(sizeof(EdgeType) * nextEdgeId);
  outDegree = (int*) malloc(sizeof(int) * vertexCount);
  inDegree = (int*) malloc(sizeof(int) * vertexCount);

  // read in edge data.
  for (int i = 0; i < nextEdgeId; i++) {
    edgeData[i] = temp_edgeData[i];
  }

  int outEdgePosition = 0;
  int inEdgePosition = 0;
  for (int i = 0; i < vertexCount; i++) {
    // record the start index of this vertex into the two arrays.
    out_edge_index[i] = outEdgePosition;
    in_edge_index[i] = inEdgePosition;
   
    vertexData[i] = temp_vertexData[i];
    outDegree[i] = temp_outDegree[i];
    inDegree[i] = temp_inDegree[i];   
 
    // read temporary out edges out into finalized array.
    for (int j = 0; j < outDegree[i]; j++) {
      out_edges[outEdgePosition++] = temp_out_edges[i][j];
    }

    // read temporary in edges out into finalized array.
    for (int j = 0; j < inDegree[i]; j++) {
      in_edges[inEdgePosition++] = temp_in_edges[i][j];
    }
  }
}

template<typename VertexType, typename EdgeType>
struct edge_info* Graph< VertexType,  EdgeType>::getOutEdges(int vid){
  return &out_edges[out_edge_index[vid]];
}

template<typename VertexType, typename EdgeType>
struct edge_info* Graph< VertexType,  EdgeType>::getInEdges(int vid){
  return &in_edges[in_edge_index[vid]];
}

template<typename VertexType, typename EdgeType>
VertexType* Graph< VertexType,  EdgeType>::getVertexData(int vid){
  return &vertexData[vid];
}

template<typename VertexType, typename EdgeType>
EdgeType* Graph< VertexType,  EdgeType>::getEdgeData(int eid){
  return &edgeData[eid];
}

template<typename VertexType, typename EdgeType>
int Graph< VertexType,  EdgeType>::getVertexColor(int vid){
  return vertexColors[vid];
}

template<typename VertexType, typename EdgeType>
int Graph< VertexType,  EdgeType>::getOutDegree(int vid){
  return outDegree[vid];
}

template<typename VertexType, typename EdgeType>
int Graph< VertexType,  EdgeType>::getInDegree(int vid){
  return inDegree[vid];
}

template<typename VertexType, typename EdgeType>
int Graph< VertexType,  EdgeType>::num_vertices(){
  return vertexCount;
}

template<typename VertexType, typename EdgeType>
int Graph< VertexType,  EdgeType>::num_edges(){
  return nextEdgeId;
}

template<typename VertexType, typename EdgeType>
int Graph< VertexType,  EdgeType>::compute_coloring(){
  // perform a parallel coloring of the graph.
  std::vector<std::pair<int, int> > permutation(vertexCount);
  
  cilk_for(int v = 0; v < vertexCount; ++v){ 
    permutation[v] = std::make_pair(-(inDegree[v] + outDegree[v]), v);
  }
  std::sort(permutation.begin(), permutation.end());
   
  int* order = (int*) malloc(sizeof(int) * vertexCount);
  vertexColors = (int*) malloc(sizeof(int) * vertexCount);

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
    if (max_success + prefix_length > vertexCount){
      stop = vertexCount;
    } else {
      stop = max_success + prefix_length;
    }

    cilk::reducer_min<int> min(max_success + prefix_length);
    cilk::holder< std::set<int> > neighbor_colors_holder;
    cilk_for(int v = max_success; v < stop; v++) {
      int vid = permutation[v].second;
      if (vertexColors[vid] == -1) {
        min.calc_min(v);
        bool skip = false;
        std::set<int> neighbor_colors = neighbor_colors_holder();

        struct edge_info* in_edges = getInEdges(vid);
        struct edge_info* out_edges = getOutEdges(vid);

        for (int i = 0; i < inDegree[vid]; i++) {
          if (order[in_edges[i].neighbor_id] < v && vertexColors[in_edges[i].neighbor_id] == -1) {
            skip = true;
            break;
          }
          neighbor_colors.insert(vertexColors[in_edges[i].neighbor_id]); 
        }
        
        if (!skip) {
          for (int i = 0; i < outDegree[vid]; i++) {
            if (order[out_edges[i].neighbor_id] < v && vertexColors[out_edges[i].neighbor_id] == -1) {
              skip = true;
              break;
            }
            neighbor_colors.insert(vertexColors[out_edges[i].neighbor_id]); 
          }
        }

        if (skip) {
          done = false;
        } else {
          int chosen_color = -1;
          for (int j = 0; j < vertexCount; j++ ) { 
            if (neighbor_colors.find(j) == neighbor_colors.end()){
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

  printf("the max color %d, num_iterations: %d \n", the_max_color, num_iterations);
  return the_max_color+1;
}


