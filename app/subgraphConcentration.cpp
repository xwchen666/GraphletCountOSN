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

  vector<string> randomwalkType = {"Basic", "PSRW"};

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
  unsigned long long sum = 0;
  for (auto value : truecount)
    sum += value;
  vector<double> trueconcentraion(truecount.size(), 0);
  for (int i = 0; i < truecount.size(); ++i)
    trueconcentraion[i] = 1.0 * truecount[i] / sum;

  for (size_t i = 0; i < truecount.size(); ++i) {
    printf("type i = %-2d, concentration: %-12f\n", i, trueconcentraion[i]);
  }

  // log files
  for (size_t idx = 0; idx < errorType.size(); ++idx) {
    // relative error
    string resultfilename = "data/concentration/node" + std::to_string(size_k) +
                            "Concentration" + errorType[idx] + ".txt";
    flogs[idx] = fopen(resultfilename.c_str(), "a+");
    if (flogs[idx] == NULL) {
      printf("Cannot open the %s file\n", resultfilename.c_str());
      exit(-1);
    }
  }

  double time_count_IMPR = 0.0;
  double time_count_PSRW = 0.0;

  srand(time(NULL));
  int randomwalk_count = randomwalkType.size();
  vector<std::ostringstream> buffers(randomwalk_count);
  for (int count = 20; count <= 20; count += 1) {
    int steps = count * 1000;
    vector<vector<vector<double>>> simulations(
      randomwalk_count, vector<vector<double>>(testsize));
    for (int i = 0; i != testsize; ++i) {
      vector<vector<double>> estimates(randomwalk_count);
      if (size_k == 3) {
        time_count_IMPR -= GetCurrentTimeSec();
        estimates[0] = node3Concentration_B(edges, steps); // proposed
        time_count_IMPR += GetCurrentTimeSec();
        time_count_PSRW -= GetCurrentTimeSec();
        estimates[1] = node3SRW2(edges, steps); // PSRW
        time_count_PSRW += GetCurrentTimeSec();
      } else if (size_k == 4) {
        time_count_IMPR -= GetCurrentTimeSec();
        estimates[0] = node4Concentration_B(edges, steps); // proposed
        time_count_IMPR += GetCurrentTimeSec();
        time_count_PSRW -= GetCurrentTimeSec();
        estimates[1] = node4SRW3(edges, steps); // PSRW
        time_count_PSRW += GetCurrentTimeSec();
      } else if (size_k == 5) {
        time_count_IMPR -= GetCurrentTimeSec();
        estimates[0] = node5Concentration_B(edges, steps); // proposed
        time_count_IMPR += GetCurrentTimeSec();
        time_count_PSRW -= GetCurrentTimeSec();
        estimates[1] = node5SRW4(edges, steps); // PSRW
        time_count_PSRW += GetCurrentTimeSec();
      }
      for (int idxz = 0; idxz < randomwalk_count; ++idxz)
        simulations[idxz][i] = estimates[idxz];
    }
    // NRMSE
    vector<vector<double>> nrmses(randomwalk_count,
                                  vector<double>(truecount.size(), 0));
    // MRE
    vector<vector<double>> relative_error(randomwalk_count,
                                          vector<double>(truecount.size(), 0));
    // Mean of multiple runs
    vector<vector<double>> relative_value(randomwalk_count,
                                          vector<double>(truecount.size(), 0));
    // lower bound of 1000 runs
    vector<vector<double>> lower_bound(randomwalk_count,
                                       vector<double>(truecount.size(), 0));
    // upper bound of 1000 runs
    vector<vector<double>> upper_bound(randomwalk_count,
                                       vector<double>(truecount.size(), 0));
    // compute all of these error metric
    for (int idxz = 0; idxz < randomwalk_count; ++idxz) {
      computeError(simulations[idxz], trueconcentraion, relative_value[idxz],
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
  printf("%s k = %d, IMPR = %lf\n, PSRW = %lf\n", argv[1], size_k,
         time_count_IMPR / (1.0 * testsize),
         time_count_PSRW / (1.0 * testsize));
  return 0;
}
