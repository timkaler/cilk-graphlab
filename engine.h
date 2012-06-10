#include <iostream>
#include <fstream>
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
#include "scheduler.h"
#ifndef ENGINE_H
#define ENGINE_H



template<typename VertexType, typename EdgeType>
class engine {
  private:
    Graph<VertexType, EdgeType>* graph;
    Scheduler* scheduler;
  public:
    engine(Graph<VertexType, EdgeType>* graph, Scheduler* scheduler);
    void run();
    void process_update_task(Scheduler::update_task task);
    void process_update_tasks(const Scheduler::update_task* tasks, int taskCount);
    void parallel_process_pennant(Pennant<Scheduler::update_task>* p, int fillSize);
    void parallel_process(Bag<Scheduler::update_task>* bag); 
};

#endif
