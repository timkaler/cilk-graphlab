// Copyright (c) 2013, Tim Kaler - MIT License

#ifndef MULTIBAG_H_
#define MULTIBAG_H_

#include <cilk/cilk.h>
#include <cilk/reducer.h>
#include <cilk/reducer_list.h>
#include <cilk/cilk_api.h>
#include <stdlib.h>
#include <assert.h>
#include <vector>
#include <list>
#include "./spa_table.h"

template <typename T>
class Multibag {
 public:
  explicit Multibag(int _ncolors);
  void insert(T t, int c);
  std::vector<std::vector<T>*> get_vector_list(int color);
  inline uint64_t numElements(int color);
  void collect();
  inline bool isEmpty(int color);
  inline void clear();
  std::vector<std::vector<std::vector<T>*> > collected_vectors;
  std::vector<SPA_table<T> > wl_vectors;
  int total_size;
  int ncolors;
};

template <typename T>
Multibag<T>::Multibag(int _ncolors) {
  ncolors = _ncolors;
  total_size = 0;
  wl_vectors.resize(__cilkrts_get_nworkers());
  for (int i = 0; i < __cilkrts_get_nworkers(); i++) {
    wl_vectors[i] = SPA_table<T>();
    wl_vectors[i].view_array.resize(ncolors);
    for (int j = 0; j < ncolors; j++) {
      wl_vectors[i].view_array[j] = NULL;
    }
    wl_vectors[i].log_array.clear();
  }
  collected_vectors.resize(ncolors);
}

template <typename T>
void Multibag<T>::insert(T el, int color) {
  int wid = __cilkrts_get_worker_number();
  if (wl_vectors[wid].view_array[color] == NULL) {
    std::vector<T>* local_v = new std::vector<T>();
    wl_vectors[wid].log_array.push_back(local_v);
    wl_vectors[wid].log_colors_array.push_back(color);
    wl_vectors[wid].view_array[color] = local_v;
  }
  wl_vectors[wid].view_array[color]->push_back(el);
}

template <typename T>
void Multibag<T>::collect() {
  collected_vectors.clear();
  collected_vectors.resize(ncolors);

  int* count_array = reinterpret_cast<int*>(calloc(ncolors, sizeof(int)));

  // Count the vectors for each color.
  cilk_for (int i = 0; i < __cilkrts_get_nworkers(); i++) {
    cilk_for (int j = 0; j < wl_vectors[i].log_array.size(); j++) {
      int c = wl_vectors[i].log_colors_array[j];
      __sync_fetch_and_add(&(count_array[c]), 1);
    }
  }

  // Count array now tells us how many of each color there is.
  cilk_for (int i = 0; i < ncolors; i++) {
    if (count_array[i] > 0) {
      collected_vectors[i].resize(count_array[i]);
    }
  }

  cilk_for (int i = 0; i < __cilkrts_get_nworkers(); i++) {
    cilk_for (int j = 0; j < wl_vectors[i].log_array.size(); j++) {
      int c = wl_vectors[i].log_colors_array[j];
      int index = __sync_sub_and_fetch(&(count_array[c]), 1);
      std::vector<T>* local_v = wl_vectors[i].log_array[j];
      collected_vectors[c][index]=(local_v);

      // reset this worker's view entry
      wl_vectors[i].view_array[c] = NULL;
    }
    // reset this worker's log.
    wl_vectors[i].log_array.clear();
    wl_vectors[i].log_colors_array.clear();
  }

  free(count_array);
}

template <typename T>
std::vector<std::vector<T>*> Multibag<T>::get_vector_list(int color) {
  return collected_vectors[color];
}

template <typename T>
inline uint64_t
Multibag<T>::numElements(int color) {
  int total_size = 0;
  for (int i = 0; i < collected_vectors[color].size(); i++) {
    total_size += collected_vectors[color][i]->size();
  }
  return total_size;
}

template <typename T>
inline bool Multibag<T>::isEmpty(int color) {
  return numElements(color) == 0;
}

template <typename T>
inline void Multibag<T>::clear() {
  cilk_for (int c = 0; c < ncolors; c++) {
    std::vector<std::vector<T>*> list =  get_vector_list(c);
    cilk_for (int i = 0; i < list.size(); i++) {
      delete list[i];
    }
    collected_vectors[c].clear();
  }
}

#endif  // MULTIBAG_H_
