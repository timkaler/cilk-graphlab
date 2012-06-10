#include "scheduler.h"
#ifndef SCHEDULER_CPP
#define SCHEDULER_CPP
Scheduler::Scheduler(int* vertexColors, int colorCount, int vertexCount) {
  // create bags for each color. 
  currentBags = (Bag_reducer<update_task>**) malloc (sizeof(Bag_reducer<update_task>*) * colorCount);
  nextBags = (Bag_reducer<update_task>**) malloc (sizeof(Bag_reducer<update_task>*) * colorCount);
  currentColor = colorCount;
  this->vertexColors = vertexColors;
  this->colorCount = colorCount;
  this->numVertices = vertexCount;
  for (int i = 0; i < colorCount; i++) {
    currentBags[i] = new Bag_reducer<update_task>();
    nextBags[i] = new Bag_reducer<update_task>();
  }
  vertex_task_added = (int*) calloc(numVertices, sizeof(int));
}

void Scheduler::add_task(int vid, void (*the_update_function) (int, void*)) {
  if (vertex_task_added[vid] == 0 &&
      __sync_bool_compare_and_swap(&vertex_task_added[vid], 0, 1)) {
    update_task t;
    t.vid = vid;
    t.update_fun = the_update_function;
    assert(vertexColors[vid] < colorCount && vertexColors[vid]>=0);
    nextBags[vertexColors[vid]]->insert(t);
  }
}

void Scheduler::collect_tasks() {
// do nothing
}
Bag<Scheduler::update_task>* Scheduler::get_task_bag() {
  if (currentColor >= colorCount) {
    // this->collect_tasks();
    bool empty = true;
    delete vertex_task_added;
    vertex_task_added = (int*) calloc(numVertices, sizeof(int));
    // delete all the current bags
    for (int i = 0; i < colorCount; i++) {
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
      return get_task_bag();
    }
  }
  currentColor++;
  return &(currentBags[currentColor-1]->get_reference());   
}
#endif
