// Our proposed algorithm framework
#include "IMPR.h"

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

using std::pair;
using std::unordered_set;
using std::vector;

/*
 * 3 node: Triangle and Line
 */

vector<vector<double>> node3Count_B(const vector<vector<int>> &edges, int steps,
                                    long long V, long long E) {
  const int n = edges.size();
  vector<vector<double>> counts(2, vector<double>(2, 0));
  int xk0 = rand() % n; // randomly choose a start point
  int xk = edges[xk0][rand() % edges[xk0].size()];
  double tri = 0.0, li = 0.0;
  double reciprocal = 0.0;
  reciprocal += 1.0 / edges[xk0].size();
  for (int s = 0; s < steps; ++s) {
    // now we are at edge xk0->xk
    // we count the number of triangles
    reciprocal += 1.0 / edges[xk].size();
    int triangle = 0, line = 0;
    {
      auto first1 = edges[xk0].begin(), last1 = edges[xk0].end();
      auto first2 = edges[xk].begin(), last2 = edges[xk].end();
      while (first1 != last1 && first2 != last2) {
        if (*first1 < *first2) { // 1 < 2
          ++first1;
        } else {
          if (!(*first2 < *first1)) { // 1 == 2
            ++triangle;
            ++first1;
          }
          ++first2;
        }
      }
    }
    line = edges[xk0].size() + edges[xk].size() - 2 - 2 * triangle;
    // beta31 = 3, beta32 = 2
    tri += triangle / 3.0;
    li += line / 2.0;
    xk0 = xk;
    xk = edges[xk][rand() % edges[xk].size()];
  }
  counts[0][0] = li * E / steps, counts[0][1] = tri * E / steps;
  double estimated_e = 1.0 * V * (steps + 1) / reciprocal / 2;
  counts[1][0] = li * estimated_e / steps,
  counts[1][1] = tri * estimated_e / steps;
  return counts;
}

vector<double> node3Concentration_B(vector<vector<int>> &edges, int steps) {
  vector<vector<double>> tmp = node3Count_B(edges, steps, edges.size(), 1);
  double sum = 0.0;
  for (auto value : tmp[0])
    sum += value;
  vector<double> concentration(tmp[0].size(), 0);
  for (size_t idx = 0; idx < tmp[0].size(); ++idx) {
    concentration[idx] = tmp[0][idx] / sum;
  }
  return concentration;
}

/*
 * 4 nodes
 */

// helper function area
/* -------------------------------------------*/
/**
 * @brief random walk on edgegraph to estimate four-node motif concentration
 *
 * @param edges original graph edges
 * @param u u-v-w are nodes
 * @param v
 * @param w
 * @param concentration u-v-w and their neighbors's contribution to the
 * concentration
 */
/* -------------------------------------------*/
void count_4_motif(const vector<vector<int>> &edges, int u, int v, int w,
                   vector<double> &concentration) {
  static vector<int> ind(edges.size(), 0);
  int C_side = 0, C_center = 0;
  int C_sidecenter = 0, C_sideside = 0, C_all = 0;
  bool uw_connected = false;
  auto first1 = edges[u].begin(), last1 = edges[u].end();
  auto first2 = edges[w].begin(), last2 = edges[w].end();
  int uppper_bound = std::max(edges[u].back(), edges[w].back());
  uppper_bound = std::min(uppper_bound, edges[v].back());
  for (auto itr = edges[v].begin(); itr != edges[v].end(); ++itr) { // set to 1
    if (*itr > uppper_bound)
      break;
    if (*itr == w || *itr == u)
      continue;
    ind[*itr] = 1;
  } // set value
  while (true) {
    if (first1 == last1) {
      while (first2 != last2) {
        if (*first2 > uppper_bound)
          break;
        if (ind[*first2] == 1)
          ++C_sidecenter;
        ++first2;
      }
      break; // jump out of the forloop
    }
    if (first2 == last2) {
      while (first1 != last1) {
        if (*first1 > uppper_bound)
          break;
        if (ind[*first1] == 1)
          ++C_sidecenter;
        ++first1;
      }
      break; // jump out of the forloop
    }
    // mergesort-like scan of the two vectors
    if (*first1 < *first2) { // 1 < 2
      if (ind[*first1] == 1)
        ++C_sidecenter;
      if (*first1 == w)
        uw_connected = true;
      ++first1;
    } else {
      if (!(*first2 < *first1)) { // 1 == 2
        if (*first1 != v) {
          if (ind[*first1] == 1)
            ++C_all;
          else
            ++C_sideside;
        }
        ++first1;
      } else { // 1 > 2
        if (ind[*first2] == 1)
          ++C_sidecenter;
        if (*first2 == u)
          uw_connected = true;
      }
      ++first2;
    }
  } // end of while
  C_side = edges[u].size() + edges[w].size() - C_sidecenter - 2 * C_all -
           2 * C_sideside - 2 - 2 * uw_connected;
  C_center = edges[v].size() - C_sidecenter - C_all - 2;

  for (auto itr = edges[v].begin(); itr != edges[v].end(); ++itr) {
    if (*itr > uppper_bound)
      break;
    ind[*itr] = 0;
  } // clear
  if (uw_connected) {
    concentration[M4_TAILEDTRIANGLE] += C_side + C_center; // tailedtriangle
    concentration[M4_CHORDALCYCLE] += C_sidecenter + C_sideside; // chordalcycle
    concentration[M4_CLIQUE] += C_all;                           // clique
  } else {
    concentration[M4_PATH] += C_side;                 // path
    concentration[M3_STAR] += C_center;               // star
    concentration[M4_CYCLE] += C_sideside;            // cycle
    concentration[M4_TAILEDTRIANGLE] += C_sidecenter; // tailedtriangle
    concentration[M4_CHORDALCYCLE] += C_all;          //  chordalcycle
  }
}

/* -------------------------------------------*/
/**
 * @brief update the ind and counter
 *
 * @param edges original edges
 * @param ind visited flag
 * @param counter counter for different 4-bit integer
 * @param oldu old u that should be replaced
 * @param w new node w
 * @param upos old position of old u, this position should be
 * overwritten by node w
 */
/* -------------------------------------------*/
void update_bit_count_4_motif(const vector<vector<int>> &edges,
                              vector<int> &ind,     // visited flag
                              vector<int> &counter, // bitCounter
                              int oldu, int w, int upos) {
  // update the bitCounter
  // clear bit of u's neighbor node
  for (auto e : edges[oldu]) {
    assert(ind[e] != 0);
    --counter[ind[e]];
    ind[e] &= (~(1 << upos)) & (0x7); // reset bit for node u's neighbors
    if (ind[e] != 0)
      ++counter[ind[e]];
  }
  // set bit of node w's neighbor node
  for (auto e : edges[w]) {
    if (ind[e] != 0)
      --counter[ind[e]];
    ind[e] |= (1 << upos); // set bit for node z's neighbors
    ++counter[ind[e]];
  }
}

/* -------------------------------------------*/
/**
 * @brief convert bitCounter to motifCounter
 *
 * @param ind visited bits flag
 * @param tmp_count bitCounter
 * @param u node u
 * @param v node v
 * @param w node w
 * @param rfs right shift steps
 *
 * @return motif counter
 */
/* -------------------------------------------*/
vector<int> bitCounter_motifCounter4(const vector<int> &ind,
                                     const vector<int> &counter,
                                     const vector<int> &nodes) {
  // once we have tmp counter,
  // we can convert it to the motif counter
  vector<int> tmp_count = counter;

  vector<int> subgraph_3_node(3);
  for (int i = 0; i < 3; ++i) {
    // right shift to get the neighborhood nodes
    subgraph_3_node[i] = ind[nodes[i]];
    // remove redundant count! (remove count of u, v, w)
    --tmp_count[subgraph_3_node[i]];
  }
  vector<int> motif_4_counter(6, 0);
  for (int i = 1; i < 8; ++i) {
    vector<int> degree_signature(4, 0);
    for (int j = 0; j < 3; ++j) { // partial degree of node u, v, w
      degree_signature[j] = __builtin_popcount(subgraph_3_node[j]);
      // builtin function of gcc (# of one bits)
    }
    for (int j = 0; j < 3; ++j) { // update degree info
      bool bit = (i >> j) & 1;
      if (bit) {
        ++degree_signature[3];
        ++degree_signature[j];
      }
    }
    // check the motif type
    int type = determine_motif_type(hash_signature(degree_signature), 4);
    assert(tmp_count[i] >= 0);
    motif_4_counter[type] += tmp_count[i];
  }
  return motif_4_counter;
}

vector<vector<double>> node4Count_B(const vector<vector<int>> &edges, int steps,
                                    long long V, long long E) {
  const int n = edges.size();
  int xk0 = rand() % n; // randomly choose a start point
  int xk = edges[xk0][rand() % edges[xk0].size()];
  double reciprocal = 0;
  reciprocal += 1.0 / edges[xk0].size() + 1.0 / edges[xk].size();
  // beta / 2
  vector<double> reweight_vec = {2, 3, 4, 5, 8, 12};
  vector<vector<double>> counts(4, vector<double>(6, 0));
  for (int i = 0; i < steps; ++i) {
    // now we are at edge xk0->xk
    vector<double> tmp_count(6, 0);
    // propose a new node
    int xk1 = edges[xk][rand() % edges[xk].size()];
    reciprocal += 1.0 / edges[xk1].size();

    if (xk1 != xk0) {
      count_4_motif(edges, xk0, xk, xk1, tmp_count);
      double weight = edges[xk].size();
      if (std::binary_search(edges[xk0].begin(), edges[xk0].end(), int(xk1))) {
        weight = 3.0 / (1.0 / edges[xk0].size() + 1.0 / edges[xk].size() +
                        1.0 / edges[xk1].size());
      }
      // reweight
      for (int idx = 0; idx < 6; ++idx) {
        counts[0][idx] += tmp_count[idx] * edges[xk].size();
        counts[1][idx] += tmp_count[idx] * weight;
      }
    }
    xk0 = xk;
    xk = xk1;
  }
  counts[2] = counts[0];
  counts[3] = counts[1];
  double estimated_e = 1.0 * V * (steps + 2) / reciprocal / 2;
  for (int r = 0; r < int(counts.size()); ++r) {
    for (int i = 0; i < 6; ++i) {
      if (r < 2) {
        counts[r][i] = counts[r][i] * E / steps / reweight_vec[i];
      } else {
        counts[r][i] = counts[r][i] * estimated_e / steps / reweight_vec[i];
      }
    }
  }
  return counts;
}

vector<double> node4Concentration_B(const vector<vector<int>> &edges,
                                    int steps) {
  vector<vector<double>> tmp = node4Count_B(edges, steps, edges.size(), 1);
  double sum = 0;
  for (auto value : tmp[1])
    sum += value;
  vector<double> concentration(tmp[1].size(), 0);
  for (size_t idx = 0; idx < tmp[1].size(); ++idx) {
    concentration[idx] = tmp[1][idx] / sum;
  }
  return concentration;
}

// helper function area
/* -------------------------------------------*/
/**
 * @brief reweight for the four-node motif if we run simple random walk
 *
 * @param edges edgelist of original graph
 * @param v1 node 1
 * @param v2 node 2
 * @param v3 node 3
 * @param v4 node 4
 * @param type motif type
 *
 * @return  weight coefficient
 */
/* -------------------------------------------*/
double reweight(const vector<vector<int>> &edges, int v1, int v2, int v3,
                int v4, Motif_4_type type) {
  double weight = 0.0;
  double d1 = edges[v1].size(), d2 = edges[v2].size(), d3 = edges[v3].size(),
         d4 = edges[v4].size();
  switch (type) {
  case M4_PATH:
    weight = 1.0 / (1.0 / d2 + 1.0 / d3);
    break;
  case M3_STAR:
    weight = d2 / 3.0;
    break;
  case M4_CYCLE:
    weight = 1.0 / (1.0 / d1 + 1.0 / d2 + 1.0 / d3 + 1.0 / d4);
    break;
  case M4_TAILEDTRIANGLE:
    weight = 1.0 / (1.0 / d1 + 1.0 / d4 + 3.0 / d2);
    break;
  case M4_CHORDALCYCLE:
    weight = 1.0 / (3.0 / d1 + 3.0 / d3 + 1.0 / d2 + 1.0 / d4);
    break;
  case M4_CLIQUE:
    weight = 1.0 / (3.0 / d1 + 3.0 / d2 + 3.0 / d3 + 3.0 / d4);
    break;
  default:
    break;
  }
  return weight;
}

double reweight(const vector<vector<int>> &edges, vector<int> &tmpNodes,
                Motif_4_type type) {
  int v1 = tmpNodes[0], v2 = tmpNodes[1], v3 = tmpNodes[2], v4 = tmpNodes[3];
  return reweight(edges, v1, v2, v3, v4, type);
}

/*
 * 5 nodes
 */

// helper function area
void write_bitCounter(const vector<vector<int>> &edges, vector<int> &ind,
                      vector<int> &nodes) {
  for (int i = 0; i < nodes.size(); ++i) {
    for (auto e : edges[nodes[i]]) {
      ind[e] |= (1 << i);
    }
  }
}

/* -------------------------------------------*/
/**
 * @brief  update bitCounter for 5 node motifs
 * update ind and counter
 *
 * @param edges original edges of graphs
 * @param ind visited flag
 * @param counter counter for all 5-bit integers
 * @param oldu old node u
 * @param z latest node z
 * @param upos position of old node u
 */
/* -------------------------------------------*/
void update_bit_count_5_motif(const vector<vector<int>> &edges,
                              vector<int> &ind,     // visited flag
                              vector<int> &counter, // bitCounter
                              int oldu, int z, int upos) {
  // update the bitCounter
  // clear bit of u's neighbor node
  for (auto e : edges[oldu]) {
    assert(ind[e] != 0);
    --counter[ind[e]];
    ind[e] &= (~(1 << upos)) & (0xF); // reset bit for node u's neighbors
    if (ind[e] != 0)
      ++counter[ind[e]];
  }
  // set bit of node w's neighbor node
  for (auto e : edges[z]) {
    if (ind[e] != 0)
      --counter[ind[e]];
    assert(ind[e] == 0 || counter[ind[e]] >= 0);
    ind[e] |= (1 << upos); // set bit for node z's neighbors
    ++counter[ind[e]];
  }
}

/* -------------------------------------------*/
/**
 * @brief convert bitCounter to motifCounter
 *
 * @param ind visited flag
 * @param counter counter for all 5-bit integers
 *
 * @return  count of all 5-node motifs
 */
/* -------------------------------------------*/
vector<int> bitCounter_motifCounter5(const vector<int> &ind,
                                     const vector<int> &counter,
                                     const vector<int> &nodes) {
  vector<int> tmp_count = counter;

  vector<int> subgraph_4_node(4);
  for (int i = 0; i < 4; ++i) {
    subgraph_4_node[i] = ind[nodes[i]];
    // remove redundant count! (remove count of u, v, w, z)
    assert(tmp_count[subgraph_4_node[i]] > 0);
    --tmp_count[subgraph_4_node[i]];
  }
  vector<int> motif_5_counter(21, 0);
  for (int i = 1; i < 16; ++i) {
    vector<int> degree_signature(5, 0);
    for (int j = 0; j < 4; ++j) { // partial degree of node u, v, w, z
      degree_signature[j] = __builtin_popcount(subgraph_4_node[j]);
      // builtin function of gcc (# of one bits)
    }
    vector<int> subgraph_5_node(5);
    for (int j = 0; j < 4; ++j) { // update degree info
      subgraph_5_node[j] = subgraph_4_node[j];
      bool bit = (i >> j) & 1;
      if (bit) {
        ++degree_signature[4];
        ++degree_signature[j];
        subgraph_5_node[j] |= 1 << 4;
      }
    }
    subgraph_5_node[4] = i;
    // check the motif type
    int type = determine_motif_type(hash_signature(degree_signature), 5);
    assert(type != -1);
    if (type == 4 || type == 7) {
      int k = -1; // node with degree 1
      // find node with degree 1
      for (int j = 0; j < 5; ++j) {
        if (degree_signature[j] == 1) {
          k = j;
          break;
        }
      }
      assert(k != -1);
      // find neighbor node of degree 1 node
      unsigned int x = subgraph_5_node[k];
      assert(x == 1 || x == 2 || x == 4 || x == 8 || x == 16);
      // the least significant bit that is set
      int nb =
        MultiplyDeBruijnBitPosition[((uint32_t)((x & -x) * 0x077CB531U)) >> 27];
      if (degree_signature[nb] == 2)
        type = 4;
      else {
        assert(degree_signature[nb] == 3);
        type = 7;
      }
    }
    if (type == 11 || type == 12) {
      int k = -1; // node with degree 3
      // find node with degree 3
      for (int j = 0; j < 5; ++j) {
        if (degree_signature[j] == 3) {
          k = j;
          break;
        }
      }
      assert(k != -1);
      type = 11;
      unsigned int x = subgraph_5_node[k];
      for (int j = 0; j < 5; ++j) {
        bool bit = (x >> j) & 1;
        // if the jth bit is 1 and degree of node j is 3
        // then this must be type 12
        if (bit && degree_signature[j] == 3) {
          type = 12;
          break;
        }
      }
    }
    assert(tmp_count[i] >= 0);
    motif_5_counter[type] += tmp_count[i];
  }
  return motif_5_counter;
}

/* -------------------------------------------*/
/**
 * @brief count number of different motifs containing node (u, v, w, z)
 *
 * @param edges edges of original graph
 * @param u
 * @param v
 * @param w
 * @param z current nodes (u, v, w, z)
 * @param tmp_count count of different motifs
 */
/* -------------------------------------------*/
void count_5_motif(const vector<vector<int>> &edges, const vector<int> &nodes,
                   vector<int> &tmp_count) {
  static vector<int> ind; // avoid reallocating
  vector<int> subgraph_4_node(4, 0);
  if (ind.size() == 0) {
    ind.resize(edges.size());
    std::fill(ind.begin(), ind.end(), 0);
  }
  // update ind
  for (int i = 0; i < 4; ++i) {
    for (auto e : edges[nodes[i]]) {
      ind[e] |= (1 << i);
    }
  }
  // write subgraph_4_node
  for (int i = 0; i < 4; ++i) {
    subgraph_4_node[i] = ind[nodes[i]];
    ind[nodes[i]] = 0;
  }
  vector<int> counter(16, 0);
  // update counter
  for (int i = 0; i < 4; ++i) {
    for (auto e : edges[nodes[i]]) {
      ++counter[ind[e]];
      ind[e] = 0;
    }
  }
  counter[0] = 0;
  // convert it to real motif counter
  for (int i = 1; i < 16; ++i) {
    // create a new adjacent matrix
    vector<int> subgraph_5_node(5, 0);
    std::copy(subgraph_4_node.begin(), subgraph_4_node.end(),
              subgraph_5_node.begin());
    subgraph_5_node[4] = i;
    // update 5 node subgraph
    for (int j = 0; j < 4; ++j) {
      if (((i >> j) & 1)) {
        subgraph_5_node[j] |= (1 << 4);
      }
    }
    int type = decide_motif_type<5>(subgraph_5_node);
    tmp_count[type] += counter[i];
  }
}

bool valid_permutations_SRW1_B(const vector<int> &adjacentBitSet,
                               const vector<int> &node_permutation) {
  for (size_t idx = 1; idx < node_permutation.size(); ++idx) {
    int u = node_permutation[idx - 1]; // last node
    int v = node_permutation[idx];     // current node
    // if u and v are not neighbors, we cannot walk from u to v
    if (not((adjacentBitSet[v] >> u) & 1))
      return false; // invalid sequence
  }
  return true; // valid sequence
}

template <int num>
double compute_alpha_SRW1_B(const vector<vector<int>> &edges,
                            const vector<int> &adjacentBitSet,
                            const vector<int> &nodeSet) {
  double alpha = 0;
  vector<int> nodeid;
  for (int i = 0; i < num; ++i)
    nodeid.push_back(i);
  static vector<vector<int>> permutations;
  // generate all the permutations using next_permutations in <algorithm>
  // we set the type of permutations5 as static to avoid multiple generation
  // of all 5! permutations
  if (permutations.empty()) {
    do {
      if (nodeid[0] < nodeid.back())
        permutations.push_back(nodeid);
    } while (std::next_permutation(nodeid.begin(), nodeid.end()));
  }
  // get real degrees
  vector<int> real_d;
  for (auto u : nodeSet) {
    real_d.push_back(edges[u].size());
  }
  int count = 0;
  for (auto p : permutations) {
    if (valid_permutations_SRW1_B(adjacentBitSet, p)) {
      double local_sum = 1.0;
      for (int idx = 1; idx < num - 1; ++idx)
        local_sum /= real_d[p[idx]];
      alpha += local_sum;
      ++count;
    }
  }
  return 1.0 * count / alpha;
}

vector<vector<double>> node5Count_B(const vector<vector<int>> &edges, int steps,
                                    long long V, long long E) {
  const int n = edges.size();

  // reweight coefficient
  vector<double> reweight_vec = {2,  2,  0,  5,  4,  4,  5,  6,  10, 8, 10,
                                 12, 10, 18, 18, 17, 18, 28, 28, 42, 60};

  // equation 1
  vector<int> phi1(21, 0);
  phi1[2] = 1, phi1[5] = 1, phi1[8] = 1, phi1[9] = 1, phi1[13] = 2,
  phi1[14] = 1;
  phi1[15] = 1, phi1[17] = 2, phi1[18] = 1, phi1[19] = 3, phi1[20] = 5;
  double S1 = 0;
  double reciprocal = 0.0;

  // walk 4 steps to get the initial 4 nodes
  int u = rand() % n;
  while (edges[u].size() < 2)
    u = rand() % n;
  int v = edges[u][rand() % edges[u].size()];
  int w = edges[v][rand() % edges[v].size()];
  int z = edges[w][rand() % edges[w].size()];

  // initialization: write ind and update counter
  vector<int> ind(n, 0);
  vector<int> counter(16, 0);
  vector<int> np = {u, v, w, z};
  write_bitCounter(edges, ind, np);
  for (int idx = 0; idx < n; ++idx)
    if (ind[idx] != 0)
      ++counter[ind[idx]];

  for (auto node : np) {
    int du = edges[node].size();
    S1 += du > 3 ? 1.0 * (du - 1) * (du - 2) * (du - 3) / 12.0 : 0;
    reciprocal += 1.0 / du;
  }

  int upos = 0;
  vector<vector<double>> counts(4, vector<double>(21, 0));
  vector<int> adjacentBitSet4(4, 0);
  // start random walk
  for (int s = 0; s < steps; ++s) {
    std::unordered_set<int> node_set(np.begin(), np.end());
    if (node_set.size() == 4) {
      // counter number of different motifs that contains nodeset np
      vector<int> tmp_count = bitCounter_motifCounter5(ind, counter, np);
      for (int i = 0; i < 4; ++i) {
        adjacentBitSet4[i] = ind[np[i]];
      }
      // w(X)
      int dv = edges[v].size() * edges[w].size();
      double weight = compute_alpha_SRW1_B<4>(edges, adjacentBitSet4, np);
      for (int i = 0; i < 21; ++i) {
        counts[0][i] += 1.0 * tmp_count[i] * dv;
        counts[1][i] += 1.0 * tmp_count[i] * weight;
      }
    }
    int oldu = u;
    int newnode = edges[z][rand() % edges[z].size()];
    assert(np[upos] == u);
    update_bit_count_5_motif(edges, ind, counter, oldu, newnode, upos);
    np[upos] = newnode;
    u = v, v = w, w = z, z = newnode;
    upos = (upos + 1) % 4;
    int dnewnode = edges[z].size();
    S1 += dnewnode > 3
            ? 1.0 * (dnewnode - 1) * (dnewnode - 2) * (dnewnode - 3) / 12.0
            : 0;
    reciprocal += 1.0 / dnewnode;
  }
  assert(counts[0][2] == 0);
  S1 = S1 / (steps + 4);
  for (int r = 0; r < 2; ++r) {
    for (int i = 0; i < 21; ++i) {
      if (i == 2)
        continue;
      counts[r][i] = counts[r][i] / steps / reweight_vec[i];
    }
    counts[r][2] = S1;
    for (size_t idx = 0; idx < phi1.size(); ++idx) {
      if (idx == 2)
        continue;
      counts[r][2] -= phi1[idx] * counts[r][idx];
    }
  }
  counts[2] = counts[0];
  counts[3] = counts[1];
  double estimated_e = 1.0 * V * (steps + 4) / reciprocal / 2;
  for (int r = 0; r < int(counts.size()); ++r) {
    for (int i = 0; i < 21; ++i) {
      if (r < 2) {
        counts[r][i] = counts[r][i] * E;
      } else {
        counts[r][i] = counts[r][i] * estimated_e;
      }
    }
  }
  return counts;
}

vector<double> node5Concentration_B(const vector<vector<int>> &edges,
                                    int steps) {
  vector<vector<double>> tmp = node5Count_B(edges, steps, edges.size(), 1);
  double sum = 0;
  for (auto value : tmp[1])
    sum += value;
  vector<double> concentration(tmp[1].size(), 0);
  for (size_t idx = 0; idx < tmp[1].size(); ++idx) {
    concentration[idx] = tmp[1][idx] / sum;
  }
  return concentration;
}
