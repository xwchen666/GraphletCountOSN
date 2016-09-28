#include "graphlet.h"
#include "process.h"
#include <stdio.h>
#include <stdlib.h>

#include <algorithm>
#include <vector>

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
  printf("compute graphlet statistics\n");
  double time_count = 0.0;
  time_count -= GetCurrentTimeSec();
  graphlet_core G_core(edges);
  G_core.graphlet_decomposition(8);
  time_count += GetCurrentTimeSec();
  G_core.print_connected_GFD();
  vector<unsigned long long> motifTypeCounter =
    G_core.get_graphlet_size4_stats();
  vector<string> g4_names = G_core.get_graphlet_size4_names();
  for (size_t i = 0; i < motifTypeCounter.size(); ++i) {
    printf("%s\t%lu\n", g4_names[i].c_str(), motifTypeCounter[i]);
  }
  printf("computation time %lfs\n", time_count);
  std::reverse(motifTypeCounter.begin(), motifTypeCounter.end());
  printf("true concentration \n");
  FILE *f_ground_truth = fopen("./data/motif4_ground_truth.txt", "a+");
  unsigned long long sum = 0;
  for (auto v : motifTypeCounter)
    sum += v;
  vector<double> groundtruth(6, 0);
  for (size_t idx = 0; idx < motifTypeCounter.size(); ++idx)
    groundtruth[idx] = 1.0 * motifTypeCounter[idx] / sum;
  fprintf(f_ground_truth, "%s\t", filename);
  for (int id = 0; id < 6; ++id) {
    // for debug, in reverse order
    printf("%d %.8lf\n", id, groundtruth[id]);
    // print to file
    fprintf(f_ground_truth, " %lu\t", motifTypeCounter[id]);
  }
  fprintf(f_ground_truth, "\n");
  fclose(f_ground_truth);

  return 0;
}
