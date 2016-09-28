#ifndef ORCA_H_
#define ORCA_H_
#include "process.h"
#include <algorithm>
#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <fstream>
#include <iostream>
#include <omp.h>
#include <set>
#include <unordered_map>
#include <vector>
using namespace std;

typedef long long int64;
typedef pair<int, int> PII;

struct PAIR {
  int a, b;
  PAIR(int a0, int b0) : a(min(a0, b0)), b(max(a0, b0)) {}
  bool operator<(const PAIR &rhs) const {
    return ((this->a < rhs.a) || (this->a == rhs.a && this->b < rhs.b));
  }
  bool operator==(const PAIR &rhs) const {
    return this->a == rhs.a && this->b == rhs.b;
  }
};

struct hash_PAIR {
  size_t operator()(const PAIR &x) const { return (x.a << 8) ^ (x.b << 0); }
};

struct TRIPLE {
  int a, b, c;
  TRIPLE(int a0, int b0, int c0) {
    a = a0;
    b = b0;
    c = c0;
    if (a > b)
      swap(a, b);
    if (b > c)
      swap(b, c);
    if (a > b)
      swap(a, b);
  }
  bool operator<(const TRIPLE &rhs) const {
    return std::tie(a, b, c) < std::tie(rhs.a, rhs.b, rhs.c);
  }
  bool operator==(const TRIPLE &rhs) const {
    return a == rhs.a && b == rhs.b && c == rhs.c;
  }
};

struct hash_TRIPLE {
  size_t operator()(const TRIPLE &x) const {
    return (x.a << 16) ^ (x.b << 8) ^ (x.c << 0);
  }
};

class orca {
private:
  int n, m;           // n = number of nodes, m = number of egdes
  vector<int> deg;    // degrees of individual nodes
  vector<PAIR> edges; // list of edges

  vector<vector<int>> adj; // adj[x] - adjacency list of node x
  vector<vector<PII>> inc; // inc[x] - incidence list of node x : (y, edge id)

  int64 *
    *orbit; // orbit[x][o] - how many times does node x participate in orbit o

  int *adj_matrix; // compressed adjacency matrix (each entry uses 1 bit)
  const int adj_chunk = 8 * sizeof(int);

  unordered_map<PAIR, int, hash_PAIR> common2;
  unordered_map<TRIPLE, int, hash_TRIPLE> common3;
  unordered_map<PAIR, int, hash_PAIR>::iterator common2_it;
  unordered_map<TRIPLE, int, hash_TRIPLE>::iterator common3_it;

  bool adjacent_list(int x, int y) {
    return binary_search(adj[x].begin(), adj[x].end(), y);
  }

  bool adjacent_matrix(int x, int y) {
    return adj_matrix[(x * n + y) / adj_chunk] &
           (1 << ((x * n + y) % adj_chunk));
  }
  bool (orca::*adjacent)(int, int);

  inline int common3_get(TRIPLE x) {
    return ((common3_it = common3.find(x)) != common3.end())
             ? (common3_it->second)
             : 0;
  }
  inline int common2_get(PAIR x) {
    return ((common2_it = common2.find(x)) != common2.end())
             ? (common2_it->second)
             : 0;
  }

public:
  orca(vector<vector<int>> &origin_edges) { init(origin_edges); }

  void count4();

  void count5();

  int init(vector<vector<int>> &origin_edges);

  vector<long long> countGraphlet(int g = 5);
};

#endif
