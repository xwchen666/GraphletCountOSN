#include "IMPR.h"
#include "PSRW.h"
#include "errorMetric.h"
#include "graphlet.h"
#include "process.h"
#include "utility.h"

#include <iomanip>
#include <iostream>
#include <iterator>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <time.h>
#include <vector>

using std::vector;
using std::string;

void print(vector<double> &nrmes, std::ostringstream &buffer) {
  buffer.str("");
  for (auto error : nrmes)
    buffer << std::setw(10) << error << "\t";
}

int main(int argc, char **argv) {
  if (argc < 3) {
    printf("Usage: %s graphfile size-k [testsize]\n", argv[0]);
    exit(-1);
  }
  const char *graphfile = argv[1];
  int size_k = argc > 2 ? atoi(argv[2]) : 4; // default graphlet size: 4
  int testsize =
    argc > 3 ? atoi(argv[3]) : 1000; // default simulation time: 1000

  vector<string> randomwalkType = {"SRW", "SRWNOE"};
  if (size_k != 3)
    randomwalkType = {"SRW", "SRWIMPR", "SRWNOE", "SRWIMPRNOE"};

  vector<string> errorType = {"RV", "RE", "LB", "UB", "NRMSE"};
  vector<FILE *> flogs(errorType.size(), NULL);

  int n, m; // |V|, |E|
  vector<vector<int>> edges;
  readGraph(graphfile, edges, 3); // read graph from file
  preprocess(edges, n, m);        // get the largest connected component

  // result file
  string truthfilename =
    "data/motif" + std::to_string(size_k) + "_ground_truth.txt";
  // read ground truth from file
  vector<long long> truecount =
    getMotifCounter(truthfilename.c_str(), graphfile);
  assert(!truecount.empty());
  printf("true counts\n");
  for (size_t i = 0; i < truecount.size(); ++i) {
    printf("type i = %-2d, count: %-12llu\n", i, truecount[i]);
  }
  vector<double> truecount_double(truecount.size(), 0);
  for (int i = 0; i < truecount.size(); ++i)
    truecount_double[i] = truecount[i];

  // log files
  for (size_t idx = 0; idx < errorType.size(); ++idx) {
    // relative error
    string resultfilename = "data/counts/node" + std::to_string(size_k) +
                            "RWBCount" + errorType[idx] + ".txt";
    flogs[idx] = fopen(resultfilename.c_str(), "a+");
    if (flogs[idx] == NULL) {
      printf("Cannot open the log RE file\n");
      exit(-1);
    }
  }

  srand(time(NULL));
  int randomwalk_count = randomwalkType.size();
  vector<std::ostringstream> buffers(randomwalk_count);
  for (int count = 2; count <= 20; count += 1) {
    int steps = count * 1000;
    vector<vector<vector<double>>> simulations(
      randomwalk_count, vector<vector<double>>(testsize));
    for (int i = 0; i != testsize; ++i) {
      if (size_k == 3) {
        vector<vector<double>> estimates = node3Count_B(edges, steps, n, m);
        for (int idxz = 0; idxz < randomwalk_count; ++idxz)
          simulations[idxz][i] = estimates[idxz];
      } else if (size_k == 4) {
        vector<vector<double>> estimates = node4Count_B(edges, steps, n, m);
        for (int idxz = 0; idxz < randomwalk_count; ++idxz)
          simulations[idxz][i] = estimates[idxz];
      } else if (size_k == 5) {
        vector<vector<double>> estimates = node5Count_B(edges, steps, n, m);
        for (int idxz = 0; idxz < randomwalk_count; ++idxz)
          simulations[idxz][i] = estimates[idxz];
      }
    }
    vector<vector<double>> nrmses(randomwalk_count,
                                  vector<double>(truecount.size(), 0));
    vector<vector<double>> relative_error(randomwalk_count,
                                          vector<double>(truecount.size(), 0));
    vector<vector<double>> relative_value(randomwalk_count,
                                          vector<double>(truecount.size(), 0));
    vector<vector<double>> lower_bound(randomwalk_count,
                                       vector<double>(truecount.size(), 0));
    vector<vector<double>> upper_bound(randomwalk_count,
                                       vector<double>(truecount.size(), 0));
    for (int idxz = 0; idxz < randomwalk_count; ++idxz) {
      computeError(simulations[idxz], truecount_double, relative_value[idxz],
                   relative_error[idxz], lower_bound[idxz], upper_bound[idxz],
                   nrmses[idxz]);
    }
    vector<vector<vector<double>>> errors = {relative_value, relative_error,
                                             lower_bound, upper_bound, nrmses};
    printf("---------------------------------------------\n");
    printf("Compare random walks, query step = %d\n", steps);
    for (size_t idx = 0; idx < errorType.size(); ++idx) {
      for (int idxz = 0; idxz < randomwalk_count; ++idxz) {
        print(errors[idx][idxz], buffers[idxz]);
        printf("%-10s %s = %s\n", errorType[idx].c_str(),
               randomwalkType[idxz].c_str(), buffers[idxz].str().c_str());
        fprintf(flogs[idx], "%s\t%s\t%d\t%s\n", argv[1],
                randomwalkType[idxz].c_str(), steps,
                buffers[idxz].str().c_str());
        fflush(flogs[idx]);
      }
    }
  }
  for (size_t idx = 0; idx < errorType.size(); ++idx) {
    fclose(flogs[idx]);
  }
  return 0;
}
