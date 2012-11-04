#include <iostream>
#include <stdio.h>
#include <cstdlib> 
#include <vector> 
#include <list>
#include <map>
#include "bag.h"
#include "bag.cpp"
#ifndef SCHEDULER_H
#define SCHEDULER_H
//using namespace std; 

class Scheduler {
  public:
    struct update_task{
      int vid;
      void (*update_fun) (int, void*);
    };
  private:
    Bag_reducer<update_task>** currentBags;
    Bag_reducer<update_task>** nextBags;
    int currentColor;
    int* vertexColors;
    int colorCount;
    int* vertex_task_added;
    int numVertices;
    int num_functions;
  public:
    Scheduler(int* vertexColors, int colorCount, int vertexCount, int num_functions);
    void add_task(int vid, void (*update_function) (int, void*), int num_function);
    Bag<update_task>* get_task_bag(bool* terminated);
    void collect_tasks();
};
#endif
