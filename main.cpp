#include <assert.h>
#include <sys/time.h>
#include <sys/resource.h>

#include "bag.cpp"
#include "Graph.cpp"
#include "engine.cpp"
void pagerank_update(int vid,
                     void* scheduler_void);

#include "scheduler.cpp"

// Parameters for the pagerank demo application
double termination_bound = 1e-10;
double random_reset_prob = 0.15;   // PageRank random reset probability

// Simple [0,1] RNG
double tfkRand(double fMin, double fMax)
{
    double f = (double)rand() / RAND_MAX;
    return fMin + f * (fMax - fMin);
}

// Simple timer for benchmarks
double tfk_get_time()
{
    struct timeval t;
    struct timezone tzp;
    gettimeofday(&t, &tzp);
    return t.tv_sec + t.tv_usec*1e-6;
}

// Programmer defined vertex and edge structures.
struct vdata{
  double value;
  double self_weight; // GraphLab does not support edges from vertex to itself, so
  // we save weight of vertex's self-edge in the vertex data
  vdata(double value = 1) : value(value), self_weight(0) { }
};
struct edata {
  double weight;
  double old_source_value;
  edata(double weight = 1) :
    weight(weight), old_source_value(0) { } 
};

// The scheduler and graph objects.
Scheduler* scheduler;
Graph<vdata, edata>* graph;

/**
 * Load a graph file specified in the format:
 *
 *   source_id <tab> target_id <tab> weight
 *   source_id <tab> target_id <tab> weight
 *   source_id <tab> target_id <tab> weight
 *               ....
 *
 * The file should not contain repeated edges.
 */
bool load_graph_from_file(const std::string& filename) {
  std::ifstream fin(filename.c_str());
  if(!fin.good()) return false;
  // Loop through file reading each line
  while(fin.good() && !fin.eof()) {
    size_t source = 0;
    size_t target = 0;
    float weight = -1;
    fin >> source;
    if(!fin.good()) break;
    //  fin.ignore(1); // skip comma
    fin >> target;
    assert(fin.good());
    //  fin.ignore(1); // skip comma
    fin >> weight;
    assert(fin.good());
    // Ensure that the number of vertices is correct
    if(source >= graph->num_vertices() ||
       target >= graph->num_vertices())
      graph->resize(std::max(source, target) + 1);
    if(source != target) {
      // Add the edge
      edata e(weight);
      graph->addEdge(source, target, e);
    } else {
      // add the self edge by updating the vertex weight
      graph->getVertexData(source)->self_weight = weight;
    }       
  }

  std::cout 
    << "Finished loading graph with: " << std::endl
    << "\t Vertices: " << graph->num_vertices() << std::endl
    << "\t Edges: " << graph->num_edges() << std::endl;
  graph->finalize();

//  std::cout << "Normalizing out edge weights." << std::endl;
  // This could be done in graphlab but the focus of this app is
  // demonstrating pagerank
  for(int vid = 0; 
      vid < graph->num_vertices(); ++vid) {
    vdata* vd = graph->getVertexData(vid);
    // Initialze with self out edge weight
    double sum = vd->self_weight;
    struct edge_info* out_edges = graph->getOutEdges(vid);

    // Sum up weight on out edges
    for(size_t i = 0; i < graph->getOutDegree(vid); ++i) {
      sum += graph->getEdgeData(out_edges[i].edge_id)->weight;
    }
    if (sum == 0) {
        vd->self_weight = 1.0;
        sum = 1.0; //Dangling page
    }

    assert(sum > 0);
    // divide everything by sum
    vd->self_weight /= sum;
    for(size_t i = 0; i < graph->getOutDegree(vid); ++i) {
      graph->getEdgeData(out_edges[i].edge_id)->weight /= sum;
    } 
  }
  //std::cout << "Finished normalizing edes." << std::endl;

  /*std::cout 
    << "Finalizing graph." << std::endl
    << "\t This is required for the locking protocol to function correctly"
    << std::endl;

  std::cout << "Finished finalization!" << std::endl;
*/
  return true;
} // end of load graph

/**
 * The Page rank update function
 */
void pagerank_update(int vid,
                     void* scheduler_void) {
  Scheduler* scheduler = (Scheduler*) scheduler_void;
  // Get the data associated with the vertex
  vdata* vd = graph->getVertexData(vid);
  
  // Sum the incoming weights; start by adding the 
  // contribution from a self-link.
  double sum = vd->value*vd->self_weight;

  struct edge_info* in_edges = graph->getInEdges(vid);
  int in_degree = graph->getInDegree(vid);

  for (int i = 0; i < in_degree; i++) {
    // Get the neighbor vertex value.
    vdata* neighbor_vdata = graph->getVertexData(in_edges[i].out_vertex);
    double neighbor_value = neighbor_vdata->value;

    // Get the edge data for the neighbor.
    edata* ed = graph->getEdgeData(in_edges[i].edge_id);
    
    // Compute the contribution for the neighbor, and add it to the sum.
    double contribution = ed->weight * neighbor_value;
    sum += contribution;

    // Remember this value as last read from the neighbor.
    ed->old_source_value = neighbor_value;
  }

  // compute the jumpweight
  sum = random_reset_prob/graph->num_vertices() + 
    (1-random_reset_prob)*sum;
  vd->value = sum;
  

  struct edge_info* out_edges = graph->getOutEdges(vid);
  int out_degree = graph->getOutDegree(vid);

  for (int i = 0; i < out_degree; i++) {
    edata* ed = graph->getEdgeData(out_edges[i].edge_id);
    
    // Compute edge-specific residual by comparing the new value of this
    // vertex to the previous value seen by the neighbor vertex.
    double residual = 
        ed->weight * std::fabs(ed->old_source_value - vd->value);

    // If the neighbor changed sufficiently add to scheduler.
    if (residual > termination_bound) {
      scheduler->add_task(out_edges[i].in_vertex, &pagerank_update);
    }
  }
} // end of pagerank update function

int main(int argc, char **argv)
{
  graph = new Graph<vdata, edata>();

  double load_start = tfk_get_time();
  load_graph_from_file(std::string(argv[1]));
  double load_end = tfk_get_time();


  double sum = 0;
  srand(1);
  for (int i = 0; i < graph->num_vertices(); i++) {
    graph->getVertexData(i)->value = 1 + tfkRand(0, 1);
    sum += graph->getVertexData(i)->value;
  }

  for (int i = 0; i < graph->num_vertices(); i++) {
    graph->getVertexData(i)->value = graph->getVertexData(i)->value / sum; 
  } 

  double color_start = tfk_get_time();
  int colorCount = graph->compute_coloring();
  double color_end = tfk_get_time();

  scheduler = new Scheduler(graph->vertexColors, colorCount, graph->num_vertices());
  for (int i = 0; i < graph->num_vertices(); i++){ 
    scheduler->add_task(i, &pagerank_update);
  }
  
  engine<vdata, edata>* e = new engine<vdata, edata>(graph, scheduler);
   

  double start = tfk_get_time();
  e->run();
  double end = tfk_get_time();
  printf("\n Graph Coloring: %d colors used \n", colorCount);

  printf("\n*** Benchmark Results ***\n");
  printf("Time to load graph %g \n", load_end - load_start);
  printf("Time spent coloring %f \n", (color_end-color_start));
  printf("Time spent iterating %f \n", (end-start));
  printf("Total runtime (including loading graph) %f \n", load_end + color_end + end - load_start - color_start - start);
  printf("Total runtime (not including loading graph) %f \n", color_end + end - color_start - start);
  
  printf("\n*** First 5 pagerank values ***\n");
  for (int i = 0; i < 5; i++) {
    printf("vertex %d value is %g \n", i, graph->getVertexData(i)->value);
  }
  return 0;
}

