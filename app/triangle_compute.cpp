#include "process.h"
#include <stdio.h>
#include <stdlib.h>

#include <algorithm>

using namespace std;

int main(int argc, char *argv[]) {
  if (argc < 2) {
    printf("Usage: %s filename\n", argv[0]);
    exit(-1);
  }
  const char *filename = argv[1];
  int n, m;
  vector<vector<int>> edges;
  readGraph(filename, edges); // read graph from file
  preprocess(edges, n, m);    // get the largest connected component
  printf("compute number of triangles\n");
  double time_count = 0.0;
  time_count -= GetCurrentTimeSec();
  vector<unsigned long long> motifTypeCounter;
  if (n < 3e7) {
    motifTypeCounter = count3NodeCISForward(edges, 2);
  } else {
    printf("use the naive algorithm \n");
    motifTypeCounter = count3NodeCISNaive(edges, 2);
  }
  time_count += GetCurrentTimeSec();
  unsigned long long triangle = motifTypeCounter[0];
  unsigned long long line = motifTypeCounter[1];
  printf("EdgeIterator: triangle = %llu, line = %llu, time = %.6lf\n", triangle,
         line, time_count);

  FILE *f_ground_truth = fopen("./data/motif3_ground_truth.txt", "a+");
  fprintf(f_ground_truth, "%s\t", filename);
  // print to file
  fprintf(f_ground_truth, " %llu\t%llu\t\n", line, triangle);
  fclose(f_ground_truth);
  return 0;
}
