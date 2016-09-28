#ifndef UTILITY_H_
#define UTILITY_H_

#include "process.h"
#include <algorithm>
#include <assert.h>
#include <fstream>
#include <sstream>
#include <stdlib.h>
#include <tuple>
#include <vector>

using std::vector;
using std::pair;
using std::tuple;

enum ConType { CONCENTRATION, CLUSTERINGCOEFF, LOCALCLUSTERING };

/***
// 1. M4_PATH o---o---o---o
// 2. M3_STAR
// o o
// |/
// o--o
// 3. M4_CYCLE
// o-----o
// |     |
// |     |
// o-----o
// 4. M4_TAILEDTRIANGLE
// o   o
// |\ /
// | /
// |/ \
// o---o
// 5. M4_CHORDALCYCLE
// o-----o
// |\    |
// | \   |
// |  \  |
// |   \ |
// |    \|
// o-----o
// 6. M4_CLIQUE
// o-----o
// |\   /|
// | \ / |
// |  \  |
// | / \ |
// |/   \|
// o-----o
**/
// four nodes motif types
enum Motif_4_type {
  M4_PATH,
  M3_STAR,
  M4_CYCLE,
  M4_TAILEDTRIANGLE,
  M4_CHORDALCYCLE,
  M4_CLIQUE
};
// four nodes motif name
static const char *Motif_4_name[] = {"4-PATH",         "3-STAR",
                                     "4-CYCLE",        "4-TAILEDTRIANGLE",
                                     "4-CHORDALCYCLE", "4-CLIQUE"};
// four nodes motif shape
static const char *Motif_4_shape[] = {
  "o---o---o---o\n",
  "o o\n|/\no--o\n",
  "o-----o\n|     |\n|     |\no-----o\n",
  "o   o\n|\\ /\n| /\n|/ \\\no---o\n",
  "o-----o\n|\\    |\n| \\   |\n|  \\  |\n|   \\ |\n|    \\|\no-----o\n",
  "o-----o\n|\\   /|\n| \\ / |\n|  \\  |\n| / \\ |\n|/   \\|\no-----o\n"};
// degree signature of 3 node motif
static const vector<vector<int>> Motif_3_Signature = {{1, 1, 2}, {2, 2, 2}};
// degree signature of 4 node motif
static const vector<vector<int>> Motif_4_Signature = {
  {1, 1, 2, 2}, {1, 1, 1, 3}, {2, 2, 2, 2},
  {1, 2, 2, 3}, {2, 2, 3, 3}, {3, 3, 3, 3}};
// degree signature of 5 node motif
static const vector<vector<int>> Motif_5_Signature = {
  {1, 1, 2, 2, 2}, // 0
  {1, 1, 1, 2, 3}, // 1
  {1, 1, 1, 1, 4}, // 2
  {1, 1, 2, 3, 3}, // 3
  {1, 2, 2, 2, 3}, // 4*
  {1, 1, 2, 2, 4}, // 5
  {2, 2, 2, 2, 2}, // 6
  {1, 2, 2, 2, 3}, // 7*
  {1, 2, 2, 3, 4}, // 8
  {2, 2, 2, 2, 4}, // 9
  {1, 2, 3, 3, 3}, // 10
  {2, 2, 2, 3, 3}, // 11^
  {2, 2, 2, 3, 3}, // 12^
  {2, 2, 2, 4, 4}, // 13
  {1, 3, 3, 3, 4}, // 14
  {2, 2, 3, 3, 4}, // 15
  {2, 3, 3, 3, 3}, // 16
  {2, 3, 3, 4, 4}, // 17
  {3, 3, 3, 3, 4}, // 18
  {3, 3, 4, 4, 4}, // 19
  {4, 4, 4, 4, 4}  // 20
};

// number of 4-node set the 4-node motif contains
static const vector<int> Motif4_Number_EdgeTriplet = {1, 3, 4, 5, 12, 24};

static const vector<tuple<int, int, int, int, int, int>> Motif5_Vector = {
  std::make_tuple(2, 0, 0, 0, 0, 0), // 0
  std::make_tuple(2, 1, 0, 0, 0, 0), // 1
  std::make_tuple(0, 4, 0, 0, 0, 0), // 2
  std::make_tuple(1, 0, 0, 2, 0, 0), // 3
  std::make_tuple(2, 0, 0, 1, 0, 0), // 4
  std::make_tuple(0, 2, 0, 2, 0, 0), // 5
  std::make_tuple(5, 0, 0, 0, 0, 0), // 6
  std::make_tuple(2, 1, 1, 0, 0, 0), // 7
  std::make_tuple(0, 1, 0, 2, 1, 0), // 8
  std::make_tuple(0, 0, 0, 4, 0, 0), // 9
  std::make_tuple(2, 0, 0, 1, 1, 0), // 10
  std::make_tuple(0, 2, 3, 0, 0, 0), // 11
  std::make_tuple(2, 0, 1, 2, 0, 0), // 12
  std::make_tuple(0, 2, 0, 0, 3, 0), // 13
  std::make_tuple(0, 0, 0, 3, 0, 1), // 14
  std::make_tuple(1, 0, 0, 2, 2, 0), // 15
  std::make_tuple(0, 0, 2, 2, 1, 0), // 16
  std::make_tuple(0, 0, 0, 2, 2, 1), // 17
  std::make_tuple(0, 0, 1, 0, 4, 0), // 18
  std::make_tuple(0, 0, 0, 0, 3, 2), // 19
  std::make_tuple(0, 0, 0, 0, 0, 5)  // 20
};
// 5 node motif reweight vector if we run random walk on G(2)
static const vector<int> Motif5_Reweight_Vec = {2,  5,  12, 11, 7,  16, 5,
                                                9,  25, 20, 19, 18, 16, 42,
                                                39, 35, 30, 58, 52, 84, 120};

// bit operation, position of least significant 1 bit
static const int MultiplyDeBruijnBitPosition[32] = {
  0,  1,  28, 2,  29, 14, 24, 3, 30, 22, 20, 15, 25, 17, 4,  8,
  31, 27, 13, 23, 21, 19, 16, 7, 26, 12, 18, 6,  11, 5,  10, 9};

// position of least signicant 1 bit
static const int lookuptable[32] = {-1, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1,
                                    0,  2, 0, 1, 0, 4, 0, 1, 0, 2, 0,
                                    1,  0, 3, 0, 1, 0, 2, 0, 1, 0};

// permutations of all 4 numbers
static const vector<vector<int>> permutations4 = {
  {0, 1, 2, 3}, {0, 1, 3, 2}, {0, 2, 1, 3}, {0, 2, 3, 1},
  {0, 3, 1, 2}, {0, 3, 2, 1}, {1, 0, 2, 3}, {1, 0, 3, 2},
  {1, 3, 0, 2}, {1, 2, 0, 3}, {2, 0, 1, 3}, {2, 1, 0, 3}};

inline int hash_signature(const vector<int> &signature) {
  int hash = 0;
  vector<int> local_signature(signature);
  std::sort(local_signature.begin(), local_signature.end());
  int base = 1;
  for (size_t i = 0; i < local_signature.size(); ++i) {
    hash += local_signature[i] * base;
    base *= 5;
  }
  return hash;
}

inline int determine_motif_type(int hash_value, int type = 4) {
  if (type == 3) { // three node
    for (int i = 0; i < 2; ++i) {
      if (hash_value == hash_signature(Motif_3_Signature[i]))
        return i;
    }
  }
  if (type == 4) { // four node
    for (int i = 0; i < 6; ++i) {
      if (hash_value == hash_signature(Motif_4_Signature[i]))
        return i;
    }
  }
  if (type == 5) { // five node
    for (int i = 0; i < 21; ++i) {
      if (hash_value == hash_signature(Motif_5_Signature[i]))
        return i;
    }
  }
  return -1; // invalid motif type
}

/* -------------------------------------------*/
/**
 * @brief Given bitSet of all 4 or 5 nodes, decide the motif type
 *
 * @param adjacentBitSet: bitSet of 4 or 5 nodes, each  bitSet has 4 or 5 bits
 *
 * @return type of current motif type
 */
/* -------------------------------------------*/
template <int num> int decide_motif_type(const vector<int> &adjacentBitSet) {
  vector<int> degree_signature(num, 0);
  // compute degree
  for (int j = 0; j < num; ++j) {
    degree_signature[j] = __builtin_popcount(adjacentBitSet[j]);
    assert(degree_signature[j] < num);
  }
  // determine motif type
  int type = determine_motif_type(hash_signature(degree_signature), num);
  assert(type >= 0);
  // extra criteria
  if (num == 5) {
    if (type == 4 || type == 7) {
      // find degree 1 node
      int k = -1;
      for (int j = 0; j < 5; ++j) {
        if (degree_signature[j] == 1) {
          k = j;
          break;
        }
      }
      int nb = lookuptable[adjacentBitSet[k]];
      assert(k != -1);
      if (degree_signature[nb] == 2)
        type = 4;
      else {
        assert(degree_signature[nb] == 3);
        type = 7;
      }
    }
    if (type == 11 || type == 12) {
      // find degree 3 node
      type = 11;
      int k = -1;
      for (int j = 0; j < 5; ++j) {
        if (degree_signature[j] == 3) {
          k = j;
          break;
        }
      }
      for (int j = 0; j < 5; ++j) {
        if ((adjacentBitSet[k] >> j) & 1) {
          if (degree_signature[j] == 3) {
            type = 12;
            break;
          }
        }
      }
    }
  }
  return type;
}

std::string getGraphName(const std::string &fullpath);

vector<long long> getMotifCounter(const char *truthfile, const char *graphfile);

// 2 node CIS
class Node2CIS {
public:
  int u, v;
  // constructor
  Node2CIS(int u = -1, int v = -1) : u(u), v(v) {}
  // get a random neighobor of current node
  int getRandomNeighbor(const std::vector<std::vector<int>> &edges,
                        Node2CIS &nextNode) {
    // ensure that w is adjacent to one
    int one, another;
    int e; // sampled edge
    assert(edges[u].size() + edges[v].size() > 2);
    if (edges[u].size() + edges[v].size() - 2 == 1) {
      one = edges[u].size() == 2 ? u : v;
      another = u + v - one;
      e = edges[one][0] == another ? edges[one][1] : edges[one][0];
    } else {
      // exclude u and v
      do {
        int idxw = rand() % (edges[u].size() + edges[v].size());
        one = idxw < (int)edges[u].size() ? u : v;
        idxw = one == u ? idxw : idxw - edges[u].size();
        e = edges[one][idxw];
        another = u + v - one;
      } while (e == another);
    }
    nextNode.u = one, nextNode.v = e;
    return another;
  }

  bool operator==(const Node2CIS &rhs) const {
    return (this->u == rhs.u && this->v == rhs.v) ||
           (this->v == rhs.u && this->u == rhs.v);
  }
};

class SimpleNode3CIS {
public:
  int u, v, w;
  int degree;
  vector<int> current_tri;
  SimpleNode3CIS(int u, int v, int w, int degree)
    : u(u), v(v), w(w), degree(degree) {
    current_tri = {u, v, w};
  }
  SimpleNode3CIS(vector<int> &nodes, int degree) {
    u = nodes[0], v = nodes[1], w = nodes[2];
    this->degree = degree;
    current_tri = {u, v, w};
  }
};

// 3 node CIS
class Node3CIS {
private:
  vector<vector<pair<int, int>>> nbs; // neighborhood (k-1)-node CIS
  int degree = -1;
  int isTriangle;
  bool alreadyVisited = false;

public:
  vector<int> current_tri; // u, v, w

  // constructor and deconstructor
  Node3CIS() {}
  Node3CIS(vector<int> triplet, bool isTriangle)
    : current_tri(triplet), isTriangle(isTriangle) {}
  // u-v-w, v is the centered node
  Node3CIS(int u, int v, int w, bool isTriangle) {
    current_tri.resize(3);
    current_tri[0] = u;
    current_tri[1] = v;
    current_tri[2] = w;
    this->isTriangle = isTriangle;
  }
  ~Node3CIS() {}

  /* -------------------------------------------*/
  /**
   * @brief populate the neighbors of current CIS
   *
   * @param edges edgelist of original graph
   */
  /* -------------------------------------------*/
  void populate_nbs(const vector<vector<int>> &edges);

  /* -------------------------------------------*/
  /**
   * @brief generate a random neighbor of current CIS
   *
   * @param edges edgelist of original graph
   * @param concentration concentration vector
   *
   * @return nextCIS (make sure that u,v,w and isTriangle are correct)
   */
  /* -------------------------------------------*/
  Node3CIS getRandomNeighbor(const vector<vector<int>> &edges,
                             vector<double> &concentration);

  int getRandomNeighbor(const vector<vector<int>> &edges, Node3CIS &nextCIS);

  void printAllNeighbors(int k = 3);

  int printNBS(int v);

  int getDegree() { return this->degree; }

  bool istriangle() { return this->isTriangle; }
};

// 3 node CIS using bit operation, more effcient than the 3CIS
class Node3CISBIT {
private:
  vector<vector<int>> nbs; // neighborhood (k-1)-node CIS
  vector<vector<int>> nbs_bitSet;
  vector<int> original_3_subgraph;
  vector<int> original_3_degree;
  int isTriangle;
  bool alreadyVisited = false;

public:
  vector<int> current_tri; // u, v, w

  // constructor and deconstructor
  Node3CISBIT() {}
  Node3CISBIT(vector<int> triplet) : current_tri(triplet) {}
  // u-v-w, v is the centered node
  Node3CISBIT(int u, int v, int w) {
    current_tri.resize(3);
    current_tri[0] = u;
    current_tri[1] = v;
    current_tri[2] = w;
  }
  ~Node3CISBIT() {}

  void populate_nbs(const vector<vector<int>> &edges, vector<int> &ind);

  Node3CISBIT getRandomNeighbor(vector<double> &concentration);

  void printAllNeighbors(int k = 3);

  int getDegree() {
    int degree = 0;
    for (size_t i = 0; i < current_tri.size(); ++i)
      degree += nbs[i].size();
    return degree;
  }

  bool istriangle() { return this->isTriangle; }
};

// 4 node CIS
class Node4CIS {
private:
  vector<vector<int>> nbs;         // neighborhood k-node CIS
  vector<vector<int>> nbs_bitSet;  // bitSet of these neighbor nodes
  vector<int> original_4_subgraph; // bitSet of current 4 nodes
  vector<int> original_4_degree;   // degree of current 4 nodes
  int degree;                      // degree of current cis
  int type;                        // what is the current motif type
  bool alreadyVisited = false;

public:
  vector<int> current_nodes; // u, v, w, z

  // constructor and deconstructor
  Node4CIS() {}
  Node4CIS(vector<int> nodes) : current_nodes(nodes) {}
  // u-v-w, v is the centered node
  Node4CIS(int u, int v, int w, int z) {
    current_nodes.resize(4);
    current_nodes[0] = u;
    current_nodes[1] = v;
    current_nodes[2] = w;
    current_nodes[3] = z;
  }
  ~Node4CIS() {}

  void populate_nbs(const vector<vector<int>> &edges, vector<int> &ind);

  Node4CIS getRandomNeighbor(vector<double> &concentration);

  int getDegree() {
    int degree = 0;
    for (size_t i = 0; i < nbs.size(); ++i) {
      degree += nbs[i].size();
    }
    return degree;
  }

  bool whichtype() { return this->type; }
};

#endif
