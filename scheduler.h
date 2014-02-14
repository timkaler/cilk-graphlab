// Copyright Tim Kaler 2013

#include <stdio.h>
#include <cstdlib>
#include <vector>
#include <list>
#include <map>
#include "./multibag.h"
#ifndef SCHEDULER_H_
#define SCHEDULER_H_

class Scheduler {
 public:
    struct update_task{
      int vid;
      void (*update_fun)(int, void*);
    };
 private:
    Multibag<update_task>* Q;
    int currentColor;
    int* vertexColors;
    int colorCount;
    int* vertex_task_added;
    int numVertices;
 public:
    Scheduler(int* vertexColors, int colorCount, int vertexCount);
    void add_task(int vid, void (*update_function)(int, void*));
    std::vector<std::vector<update_task>*> get_task_bag();
    void collect_tasks();
};
#endif  // SCHEDULER_H_
