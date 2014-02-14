// Copyright (c) 2013, Tim Kaler - MIT License

#include <vector>

#ifndef SPA_TABLE_H_
#define SPA_TABLE_H_

template <typename T>
struct SPA_table {
  std::vector<std::vector<T>*> view_array;
  std::vector<std::vector<T>*> log_array;
  std::vector<int> log_colors_array;
  uint64_t buffer[64];  // make sure this goes off a cache line.
};

struct SPA_reference {
  int log_index;
  int worker_id;
  int size;
  int num_elements;
  uint64_t buffer[64];  // make sure this goes off a cache line.
};

#endif  // SPA_TABLE_H_
