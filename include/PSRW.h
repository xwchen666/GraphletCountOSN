#ifndef PSRW_H_
#define PSRW_H_

#include "utility.h"
#include <vector>

using std::vector;
// PSRW
// Paper: "Efficiently Estimating Motif Statistics of Large Networks", Wang
// PingHui et al. TKDD'14

/*
 * 3 nodes
 */

vector<double> node3SRW2(const vector<vector<int>> &edges, int steps);

/*
 * 4 nodes
 */

vector<double> node4SRW3(const vector<vector<int>> &edges, int steps);

/*
 * 5 nodes
 */
vector<double> node5SRW4(const vector<vector<int>> &edges, int steps);

#endif
