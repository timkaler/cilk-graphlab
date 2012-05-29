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
  public:
    Scheduler(int* vertexColors, int colorCount, int vertexCount);
    void add_task(int vid, void (*update_function) (int, void*));
    Bag<update_task>* get_task_bag();
};
#endif
