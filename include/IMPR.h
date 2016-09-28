#ifndef IMPR_H_
#define IMPR_H_

#include "process.h"
#include "utility.h"
#include <vector>

using std::vector;
/*
 * 3 nodes
 */

vector<vector<double>> node3Count_B(const vector<vector<int>> &edges, int steps,
                                    long long V, long long E);

vector<double> node3Concentration_B(vector<vector<int>> &edges, int steps);

/*
 * 4 nodes
 */

vector<vector<double>> node4Count_B(const vector<vector<int>> &edges, int steps,
                                    long long V, long long E);

vector<double> node4Concentration_B(const vector<vector<int>> &edges,
                                    int steps);
/*
 * 5 nodes
 */

vector<vector<double>> node5Count_B(const vector<vector<int>> &edges, int steps,
                                    long long V, long long E);

vector<double> node5Concentration_B(const vector<vector<int>> &edges,
                                    int steps);

// helper functions
// we put these functions here just for test
// 4 node graphlets
void update_bit_count_4_motif(const vector<vector<int>> &edges,
                              vector<int> &ind,     // visited flag
                              vector<int> &counter, // bitCounter
                              int oldu, int w, int upos);

vector<int> bitCounter_motifCounter4(const vector<int> &ind,
                                     const vector<int> &counter,
                                     const vector<int> &nodes);

void count_4_motif(const vector<vector<int>> &edges, int u, int v, int w,
                   vector<double> &concentration);

// 5 node graphlets
void write_bitCounter(const vector<vector<int>> &edges, vector<int> &ind,
                      vector<int> &nodes);

void update_bit_count_5_motif(const vector<vector<int>> &edges,
                              vector<int> &ind,     // visited flag
                              vector<int> &counter, // bitCounter
                              int oldu, int z, int upos);

vector<int> bitCounter_motifCounter5(const vector<int> &ind,
                                     const vector<int> &counter,
                                     const vector<int> &nodes);

void count_5_motif(const vector<vector<int>> &edges, const vector<int> &nodes,
                   vector<int> &tmp_count);

#endif
