#include "orca.h"
#include "graphlet.h"
#include "process.h"
#include <stdio.h>
#include <stdlib.h>

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
  graphlet_core G_core(edges);
  G_core.graphlet_decomposition(8);
  vector<unsigned long long> g4 = G_core.get_graphlet_size4_stats();
  vector<string> g4_names = G_core.get_graphlet_size4_names();
  for (size_t i = 0; i < g4.size(); ++i) {
    printf("%s\t%lu\n", g4_names[i].c_str(), g4[i]);
  }
  orca G_orca(edges);
  vector<long long> orca_g4 = G_orca.countGraphlet(4);
  for (size_t i = 0; i < orca_g4.size(); ++i) {
    printf("%s\t%lu\n", g4_names[i].c_str(), orca_g4[5 - i]);
  }
  FILE *flog = fopen("data/orca_result.txt", "a+");
  FILE *f_time = fopen("node5_time.txt", "a+");
  double time = 0.0;
  time -= GetCurrentTimeSec();
  vector<long long> orca_g5 = G_orca.countGraphlet(5);
  time += GetCurrentTimeSec();
  printf("%s\t time = %.6lf\n", filename, time);
  fprintf(f_time, "%s\t time = %.6lf\n", filename, time);
  // fprintf(flog, "%s\t", filename);
  for (size_t i = 0; i < orca_g5.size(); ++i) {
    printf("%d\t%lu\n", i, orca_g5[i]);
    // fprintf(flog, "%lu\t", orca_g5[i]);
  }
  // fprintf(flog, "\n");
  fclose(f_time);

  return 0;
}
