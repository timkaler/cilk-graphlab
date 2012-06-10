#include "engine.h"

template<typename VertexType, typename EdgeType>
engine<VertexType, EdgeType>::engine< VertexType,  EdgeType>(
    Graph<VertexType, EdgeType>* graph, Scheduler* scheduler) {
  this->graph = graph;
  this->scheduler = scheduler;
}

template<typename VertexType, typename EdgeType>
void engine<VertexType, EdgeType>::run(){
  int iterationCount = 0;
  Bag<Scheduler::update_task>* b = scheduler->get_task_bag();
  while (b->numElements() > 0 /*&& iterationCount < 40*/) {
    iterationCount++;
    parallel_process(b); 
    b = scheduler->get_task_bag(); 
  }
}

template<typename VertexType, typename EdgeType>
  void engine<VertexType, EdgeType>::process_update_task(Scheduler::update_task task) {
    // code to process the update task here.
    task.update_fun(task.vid, scheduler); 
  } 

template<typename VertexType, typename EdgeType>
  void engine<VertexType, EdgeType>::process_update_tasks(const Scheduler::update_task* tasks, int taskCount) {
    // prefetch data
/*
    for (int i = 0; i < taskCount; i++) {
      graph->prefetch_vertex(tasks[i].vid);
    }
*/
    for (int i = 0; i < taskCount; i++) {
      process_update_task(tasks[i]);
    }
  }

template<typename VertexType, typename EdgeType>
  void engine<VertexType, EdgeType>::parallel_process_pennant(Pennant<Scheduler::update_task>* p, int fillSize) {
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

template<typename VertexType, typename EdgeType>
  void engine<VertexType, EdgeType>::parallel_process(Bag<Scheduler::update_task>* bag) {
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

