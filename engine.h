// Copyright (c) 2013, Tim Kaler - MIT License

#include <cilk/cilk.h>
#include <cilk/reducer_list.h>
#include <cilk/reducer_min.h>
#include <cilk/reducer_max.h>
#include <cilk/holder.h>
#include <stdio.h>
#include <cstdlib>
#include <vector>
#include <list>
#include <map>
#include <string>
#include <set>
#include <cmath>
#include <algorithm>
#include "./scheduler.h"

#ifndef ENGINE_H_
#define ENGINE_H_

template<typename VertexType, typename EdgeType>
class engine {
 private:
  Graph<VertexType, EdgeType>* graph;
  Scheduler* scheduler;
 public:
  engine(Graph<VertexType, EdgeType>* graph, Scheduler* scheduler);
  void run();
  void process_update_task(Scheduler::update_task task);
  void parallel_process(
      std::vector<std::vector<Scheduler::update_task>*> subbags);
};

#endif  // ENGINE_H_
