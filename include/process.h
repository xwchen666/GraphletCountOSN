#ifndef PROCESS_H_
#define PROCESS_H_

#include <stdlib.h>
#include <sys/time.h>
#include <vector>

// return current time in seconds
inline double GetCurrentTimeSec() {
  struct timeval tv;
  gettimeofday(&tv, NULL);
  return tv.tv_sec + tv.tv_usec * 1e-6;
}

// read graph from file
void readGraph(const char *filename, std::vector<std::vector<int>> &edges,
               int format = 2);

// preprocess the graph, get the largest connected component and remove other
// smaller component.
void preprocess(std::vector<std::vector<int>> &edges, int &nprime, int &mprime);

std::vector<unsigned long long>
count3NodeCISNaive(const std::vector<std::vector<int>> &edges, int vecsize = 2);

// count the exact number of triangles in the graph
// Paper: Finding, Counting and Listing all Triangles in Large Graphs, An
// Experimental Study
// time complexity: O(E^3/2)
std::vector<unsigned long long>
count3NodeCISForward(const std::vector<std::vector<int>> &edges,
                     int vecsize = 2);

double localClusteringCoeff(const std::vector<std::vector<int>> &edges);

#endif
