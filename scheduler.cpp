#include "scheduler.h"
#ifndef SCHEDULER_CPP
#define SCHEDULER_CPP
Scheduler::Scheduler(int* vertexColors, int colorCount, int vertexCount, int num_functions) {
  // create bags for each color. 
  currentBags = (Bag_reducer<update_task>**) malloc (sizeof(Bag_reducer<update_task>*) * colorCount * num_functions);
  nextBags = (Bag_reducer<update_task>**) malloc (sizeof(Bag_reducer<update_task>*) * colorCount * num_functions);
  currentColor = colorCount;
  this->vertexColors = vertexColors;
  this->colorCount = colorCount;
  this->numVertices = vertexCount;
  this->num_functions = num_functions;
  for (int i = 0; i < colorCount * num_functions; i++) {
    currentBags[i] = new Bag_reducer<update_task>();
    nextBags[i] = new Bag_reducer<update_task>();
  }
  vertex_task_added = (int*) calloc(numVertices * num_functions, sizeof(int));
}

void Scheduler::add_task(int vid, void (*the_update_function) (int, void*),
    int num_function) {
  if (vertex_task_added[vid + numVertices * num_function] == 0 &&
      __sync_bool_compare_and_swap(&vertex_task_added[vid + numVertices * num_function], 0, 1)) {
    update_task t;
    t.vid = vid;
    t.update_fun = the_update_function;
    assert(vertexColors[vid] < colorCount && vertexColors[vid]>=0);
    nextBags[vertexColors[vid] + colorCount * num_function]->insert(t);
  }
}

void Scheduler::collect_tasks() {
// do nothing
}
Bag<Scheduler::update_task>* Scheduler::get_task_bag(bool* terminated) {
  if (currentColor >= colorCount * num_functions) {
    *terminated = true;
    // this->collect_tasks();
    bool empty = true;
    delete vertex_task_added;
    vertex_task_added = (int*) calloc(numVertices * num_functions, sizeof(int));
    // delete all the current bags
    for (int i = 0; i < colorCount * num_functions; i++) {
      delete currentBags[i];
      currentBags[i] = nextBags[i];
      nextBags[i] = new Bag_reducer<update_task>();
      if (currentBags[i]->numElements() > 0){
        empty = false;
      }
    }
    if (empty) {
      return &currentBags[0]->get_reference();
    }
    currentColor = 0;
  }
  
  while (currentBags[currentColor]->numElements() == 0) {
    currentColor++; 
    if (currentColor >= colorCount) {
      return get_task_bag(terminated);
    }
  }
  currentColor++;
  return &(currentBags[currentColor-1]->get_reference());   
}
#endif
