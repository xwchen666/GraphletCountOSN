#ifndef ERROR_METRIC_H_
#define ERROR_METRIC_H_

#include <algorithm>
#include <assert.h>
#include <math.h>
#include <vector>

using std::vector;

void computeBound(std::vector<double> &simulation, double &percentile05,
                  double &percentile95) {
  int simulationSize = simulation.size();
  int idx05 = int(0.05 * simulationSize);
  int idx95 = int(0.95 * simulationSize);
  std::sort(simulation.begin(), simulation.end());
  percentile05 = simulation[idx05];
  percentile95 = simulation[idx95];
}

double NRMSE(std::vector<double> &simulation) {
  double merror = 0.0;
  for (double s : simulation) {
    merror += pow((s - 1), 2);
  }
  return sqrt(merror / simulation.size());
}

// square error for each graphlet type
std::vector<double> SE(const std::vector<double> &simulation,
                       std::vector<double> &groundtruth) {
  std::vector<double> serror(groundtruth.size(), 0.0);
  for (size_t idx = 0; idx < groundtruth.size(); ++idx) {
    serror[idx] = pow(simulation[idx] / groundtruth[idx] - 1, 2);
  }
  return serror;
}

std::vector<double> NRMSE(std::vector<std::vector<double>> &simulation,
                          std::vector<double> &groundtruth) {
  std::vector<double> merror(simulation[0].size());
  for (size_t id = 0; id < simulation.size(); ++id) {
    std::vector<double> squareerror = SE(simulation[id], groundtruth);
    for (size_t idx = 0; idx < simulation[id].size(); ++idx) {
      merror[idx] += squareerror[idx];
    }
  }
  // normalize and square root
  for (size_t idx = 0; idx < simulation[0].size(); ++idx) {
    merror[idx] /= simulation.size();
    merror[idx] = sqrt(merror[idx]);
  }
  return merror;
}

// write relative error, lower bound, upper bound
void computeError(const vector<vector<double>> &simulation,
                  vector<double> &true_count, vector<double> &relative_value,
                  vector<double> &relative_error, vector<double> &percentile05,
                  vector<double> &percentile95, vector<double> &nrmse) {
  // transpose
  vector<vector<double>> transpose(true_count.size(),
                                   vector<double>(simulation.size(), 0));
  vector<double> average_estimate(true_count.size(), 0);
  vector<double> average_error(true_count.size(), 0);

  for (size_t idx = 0; idx < simulation.size(); ++idx) {
    vector<double> squareerror = SE(simulation[idx], true_count);
    for (size_t idy = 0; idy < true_count.size(); ++idy) {
      transpose[idy][idx] = simulation[idx][idy];
      average_estimate[idy] += simulation[idx][idy];
      average_error[idy] += fabs(simulation[idx][idy] / true_count[idy] - 1);
      nrmse[idy] += squareerror[idy];
    }
  }
  for (size_t idy = 0; idy < true_count.size(); ++idy) {
    average_estimate[idy] /= simulation.size();
    average_error[idy] /= simulation.size();
    relative_value[idy] = average_estimate[idy];
    relative_error[idy] = average_error[idy];
    computeBound(transpose[idy], percentile05[idy], percentile95[idy]);
    nrmse[idy] /= simulation.size();
    nrmse[idy] = sqrt(nrmse[idy]);
  }
}

#endif
