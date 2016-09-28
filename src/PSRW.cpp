// pairwise subgraph random walk (PSRW, Wang et al.)
#include "PSRW.h"

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <time.h>

#include <algorithm>
#include <iostream>
#include <string>
#include <unordered_set>
#include <vector>

using std::unordered_set;
using std::vector;

// PSRW
// Paper: "Efficiently Estimating Motif Statistics of Large Networks", Wang
// PingHui et al. TKDD'14

/*
 * Three Node: Triangle and Lines
 */

/******************** SRW(2) ********************************/
vector<double> node3SRW2(const vector<vector<int>> &edges, int steps) {
  const int n = edges.size();
  vector<double> concentration(2, 0);
  // randomly choose a start point
  int u = rand() % n;
  int v = edges[u][rand() % edges[u].size()];
  Node2CIS currentNode(u, v);
  double triangle = 0, line = 0;
  for (int s = 0; s < steps; ++s) {
    Node2CIS nextNode;
    int another = currentNode.getRandomNeighbor(edges, nextNode);
    int w = nextNode.v;
    // w(Xs) = 1
    if (std::binary_search(edges[w].begin(), edges[w].end(), another)) {
      ++triangle;
    } else
      ++line;
    currentNode = nextNode;
  }
  // alpha31 = 6, alpha32 = 2
  concentration[0] = 3 * line / (triangle + 3 * line);
  concentration[1] = triangle / (triangle + 3 * line);
  return concentration;
}

/*
 * 4 nodes
 */

/******************** SRW(3) ********************************/

// PSRW
// regular random walk on CIS relationship graph
// Paper: "Efficiently Estimating Motif Statistics of Large Networks", Wang
// PingHui et al. TKDD'14
/* -------------------------------------------*/
/**
 * @brief subgraph random walk on G(3)
 *
 * @param edges edgelist of original graph
 * @param steps steps of random walk
 *
 * @return  estimated concentration vector
 */
/* -------------------------------------------*/
vector<double> node4SRW3(const vector<vector<int>> &edges, int steps) {
  const int n = edges.size();
  // get the first three nodes u-v-w, v is the centered nodes
  int u = 0, v = 0, w = 0;
  while (u == v || u == w || v == w) {
    u = rand() % n;
    while (edges[u].size() < 3)
      u = rand() % n;
    v = edges[u][rand() % edges[u].size()];
    w = edges[u][rand() % edges[u].size()];
  }
  bool isTriangle = std::binary_search(edges[v].begin(), edges[v].end(), w);
  Node3CIS current_cis(v, u, w, isTriangle);

  // for statistics
  vector<double> reweight_vec = {1, 3, 6, 3, 6, 6};
  vector<double> concentration(6, 0);

  for (int i = 0; i < steps; ++i) {
    current_cis.populate_nbs(edges);
    // walk to next triangle
    current_cis = current_cis.getRandomNeighbor(edges, concentration);
  }
  // Horvitz-Thompson estimator
  double sum = 0.0;
  for (int i = 0; i < 6; ++i) {
    concentration[i] /= reweight_vec[i];
    sum += concentration[i];
  }
  for (int i = 0; i < 6; ++i)
    concentration[i] /= sum;
  return concentration;
}

/*
 * 5 nodes
 */

/******************** SRW(4) ********************************/

// regular random walk on CIS relationship graph
// Paper: "Efficiently Estimating Motif Statistics of Large Networks", Wang
// PingHui et al. TKDD'14
/* -------------------------------------------*/
/**
 * @brief random walk on 4 node CIS to estimate motif
 *
 * @param edges edges of the original graph
 * @param steps number of steps for the random walk
 *
 * @return concentration vector
 */
/* -------------------------------------------*/
vector<double> node5SRW4(const vector<vector<int>> &edges, int steps) {
  const int n = edges.size();
  // get the first three nodes u-v-w, v is the centered nodes
  int u = rand() % n;
  while (edges[u].size() < 2)
    u = rand() % n;
  int v = edges[u][rand() % edges[u].size()];
  int w = edges[v][rand() % edges[v].size()];
  int z = edges[w][rand() % edges[w].size()];
  std::unordered_set<int> node_set({u, v, w, z});
  int reals = node_set.size();
  while (reals < 4) {
    int newnode = edges[z][rand() % edges[z].size()];
    u = v;
    v = w;
    w = z;
    z = newnode;
    std::unordered_set<int> node_set1({u, v, w, z});
    reals = node_set1.size();
  }
  // for statistics
  vector<double> includeNum = {2, 3, 4, 3, 3, 4, 5, 4, 4, 4, 4,
                               5, 5, 5, 4, 5, 5, 5, 5, 5, 5};
  vector<double> reweight_vec = includeNum;
  for (size_t i = 0; i < reweight_vec.size(); ++i) {
    reweight_vec[i] = (reweight_vec[i]) * (reweight_vec[i] - 1) / 2;
  }
  vector<int> ind(n, 0);
  vector<double> concentration(21, 0);

  Node4CIS current_cis(u, v, w, z);
  for (int i = 0; i < steps; ++i) {
    current_cis.populate_nbs(edges, ind);
    // walk to next triangle
    current_cis = current_cis.getRandomNeighbor(concentration);
  }
  // Horvitz-Thompson estimator
  double sum = 0.0;
  for (int i = 0; i < 21; ++i) {
    concentration[i] /= reweight_vec[i];
    sum += concentration[i];
  }
  for (int i = 0; i < 21; ++i)
    concentration[i] /= sum;
  return concentration;
}
