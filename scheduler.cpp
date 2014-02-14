// Copyright (c) 2013, Tim Kaler - MIT License

#include <vector>
#include "./scheduler.h"
#include "./multibag.h"

#ifndef SCHEDULER_CPP
#define SCHEDULER_CPP

Scheduler::Scheduler(int* vertexColors, int colorCount, int vertexCount) {
  this->vertexColors = vertexColors;
  this->colorCount = colorCount;
  this->numVertices = vertexCount;
  this->currentColor = colorCount;

  // Create the multibag
  Q = new Multibag<update_task>(colorCount);

  vertex_task_added = reinterpret_cast<int*>(calloc(numVertices, sizeof(int)));
}

void Scheduler::add_task(int vid, void (*the_update_function)(int, void*)) {
  if (vertex_task_added[vid] == 0 &&
      __sync_bool_compare_and_swap(&vertex_task_added[vid], 0, 1)) {
    update_task t;
    t.vid = vid;
    t.update_fun = the_update_function;
    assert(vertexColors[vid] < colorCount && vertexColors[vid] >= 0);
    Q->insert(t, vertexColors[vid]);
  }
}

std::vector<std::vector<Scheduler::update_task>*> Scheduler::get_task_bag() {
  if (currentColor >= colorCount) {
    Q->collect();
    bool empty = true;
    free(vertex_task_added);
    vertex_task_added =
        reinterpret_cast<int*>(calloc(numVertices, sizeof(int)));
    // delete all the current bags
    for (int i = 0; i < colorCount; i++) {
      if (Q->get_vector_list(i).size() > 0) {
        empty = false;
      }
    }
    if (empty) {
      return Q->get_vector_list(0);
    }
    currentColor = 0;
  }

  while (Q->get_vector_list(currentColor).size() == 0) {
    currentColor++;
    if (currentColor >= colorCount) {
      return get_task_bag();
    }
  }
  currentColor++;
  return Q->get_vector_list(currentColor-1);
}
#endif
