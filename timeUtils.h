// Copyright (c) 2013, Tim Kaler - MIT License

#ifndef TIMEUTILS_H_
#define TIMEUTILS_H_
// Simple timer for benchmarks
double tfk_get_time() {
    struct timeval t;
    struct timezone tzp;
    gettimeofday(&t, &tzp);
    return t.tv_sec + t.tv_usec*1e-6;
}
#endif  // TIMEUTILS_H_
