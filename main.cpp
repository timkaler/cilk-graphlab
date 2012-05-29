#include <assert.h>
#include "bag.cpp"
#include "Graph.cpp"
void pagerank_update(int vid,
                     void* scheduler);

#include "scheduler.cpp"


double termination_bound = 1e-10;
double random_reset_prob = 0.15;   // PageRank random reset probability

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
  std::cout << "Normalizing out edge weights." << std::endl;
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
  std::cout << "Finished normalizing edes." << std::endl;

  std::cout 
    << "Finalizing graph." << std::endl
    << "\t This is required for the locking protocol to function correctly"
    << std::endl;

  std::cout << "Finished finalization!" << std::endl;
  return true;
} // end of load graph

  void process_update_task(Scheduler::update_task task) {
    // code to process the update task here.
    task.update_fun(task.vid, scheduler); 
  } 

  void process_update_tasks(const Scheduler::update_task* tasks, int taskCount) {
    for (int i = 0; i < taskCount; i++) {
      process_update_task(tasks[i]);
    }
  }

  void parallel_process_pennant(Pennant<Scheduler::update_task>* p, int fillSize) {
    if (p->getLeft() != NULL) {
      cilk_spawn parallel_process_pennant(p->getLeft(), BLK_SIZE);
    }
    if (p->getRight() != NULL){
      cilk_spawn parallel_process_pennant(p->getRight(), BLK_SIZE);
    }
    process_update_tasks(p->getElements(), fillSize);
    
    cilk_sync;
    if (p->getLeft() != NULL) {
      delete p->getLeft();
    }

    if (p->getRight() != NULL) {
      delete p->getRight();
    }
  }
  void parallel_process(Bag<Scheduler::update_task>* bag) {
    Pennant<Scheduler::update_task>* p = NULL;
    if (bag->getFill() > 0) {
      bag->split(&p);
      cilk_spawn parallel_process(bag);
      parallel_process_pennant(p, BLK_SIZE);
    } else {
      process_update_tasks(bag->getFilling(), bag->getFillingSize());
    }
    cilk_sync;
    if (p != NULL) {
      delete p;
    }
  }

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
    vdata* neighbor_vdata = graph->getVertexData(in_edges[i].neighbor_id);
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
      //printf("scheduling task from wtihin pagerank function\n");
      scheduler->add_task(out_edges[i].neighbor_id, &pagerank_update);
    }
  }
} // end of pagerank update function


/*  void update_function (int vid, Scheduler* scheduler) {
    scheduler->add_task(vid, &update_function);
  }
*/

int main(int argc, char **argv)
{
  graph = new Graph<vdata, edata>();
  load_graph_from_file(std::string(argv[1]));
  double sum = 0;
  for (int i = 0; i < graph->num_vertices(); i++) {
    graph->getVertexData(i)->value = 1.5;
    sum += 1.5;
  }
  for (int i = 0; i < graph->num_vertices(); i++) {
    graph->getVertexData(i)->value = graph->getVertexData(i)->value / sum; 
  } 

  int colorCount = graph->compute_coloring(); 
  scheduler = new Scheduler(graph->vertexColors, colorCount, graph->num_vertices());
  for (int i = 0; i < graph->num_vertices(); i++){ 
    scheduler->add_task(i, &pagerank_update);
  }
  int iterationCount = 0;
  Bag<Scheduler::update_task>* b = scheduler->get_task_bag();
  while (b->numElements() > 0 /*&& iterationCount < 40*/) {
    iterationCount++;
    parallel_process(b); 
    b = scheduler->get_task_bag(); 
  }
  for (int i = 0; i < 5; i++) {
    printf("vertex %d value is %f \n", i, graph->getVertexData(i)->value);
  }
  printf("total iterations %d \n", iterationCount);
  return 0;
}
