#include "process.h"

#include <algorithm>
#include <assert.h>
#include <numeric>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <vector>

using std::vector;

// read graph from file
void readGraph(const char *filename, vector<vector<int>> &edges, int format) {
  FILE *f = fopen(filename, "r");
  if (!f) {
    printf("Cannot open file %s\n", filename);
    exit(-1);
  }
  int n, m;
  if (fscanf(f, "%d%d", &n, &m) < 2) {
    printf("Read number of nodes and edges error\n");
    exit(-1);
  }
  edges.resize(n);
  char buffer[128];
  // skip new line
  if (!fgets(buffer, sizeof(buffer), f)) {
    printf("Read error\n");
    exit(-1);
  }
  int realm = 0;
  int self_loop = 0;
  while (fgets(buffer, sizeof(buffer), f)) {
    int u, v, w = 0;
    if (format == 2) {
      sscanf(buffer, "%d%d%*lf", &u, &v);
    } else {
      sscanf(buffer, "%d%d%d", &u, &v, &w);
    }
    if (u == v) {
      ++self_loop;
      continue; // no self loop
    }
    edges[u].push_back(v);
    edges[v].push_back(u);
    ++realm;
  }
  printf("Read file finished, n = %d, m = %d\nreal m = %d, self loop = %d\n", n,
         m, realm, self_loop);
  fclose(f);
}

// nprime: size of largest connected component
// mprime: number of edges in the largest connected component
// preprocess the graph, get the largest connected component and reid the nodes
void preprocess(vector<vector<int>> &edges, int &nprime, int &mprime) {
  // find largest connected component first
  vector<int> idOfNode(edges.size(), -1);
  vector<int> comp_size;
  vector<int> que(edges.size());
  int current_id = 0;
  // connected component
  for (size_t u = 0; u < edges.size(); ++u) {
    // bfs
    if (idOfNode[u] == -1) {
      int que_h = 0, que_t = 0;
      que[que_t++] = u;
      idOfNode[u] = current_id;
      while (que_h < que_t) {
        int s = que[que_h++];
        for (auto w : edges[s]) {
          // if not visited yet
          if (idOfNode[w] == -1) {
            que[que_t++] = w;
            idOfNode[w] = current_id;
          }
        }
      }
      comp_size.push_back(que_t);
      ++current_id;
    }
  }
  assert(current_id == (int)comp_size.size());
  auto it = std::max_element(comp_size.begin(), comp_size.end());
  int largest_id = it - comp_size.begin();
  printf("largest component size = %d, largest_id = %d\n",
         comp_size[largest_id], largest_id);
  // remove other smaller components
  vector<int> newid(edges.size());
  size_t ptr = 0;
  for (size_t u = 0; u < edges.size(); ++u) {
    if (idOfNode[u] == largest_id) {
      newid[u] = ptr;
      if (ptr != u)
        edges[ptr] = edges[u];
      ++ptr;
    }
  }
  assert(ptr == (size_t)comp_size[largest_id]);
  edges.resize(ptr);
  // reid the graph and sort
  nprime = 0, mprime = 0;
  for (size_t u = 0; u < edges.size(); ++u) {
    for (size_t idxw = 0; idxw < edges[u].size(); ++idxw) {
      int oldw = edges[u][idxw];
      edges[u][idxw] = newid[oldw];
      assert(edges[u][idxw] != (int)u);
    }
    std::sort(edges[u].begin(), edges[u].end());
    // remove duplicate nodes
    auto itr = std::unique(edges[u].begin(), edges[u].end());
    edges[u].resize(itr - edges[u].begin());
    mprime += edges[u].size();
  }
  assert(mprime % 2 == 0);
  nprime = ptr;
  mprime /= 2;
}

// assume the edge list is sorted
// count triangles and lines O(n^3)
// is long long enough for triangle number?
vector<unsigned long long> count3NodeCISNaive(const vector<vector<int>> &edges,
                                              int vecsize) {
  vector<unsigned long long> motifCounter(vecsize, 0);
  const int n = edges.size();
  // similar to NodeIterator
  // make sure node u is the lowest rank node
  for (int u = 0; u < n; ++u) {
    for (size_t idxv = 0; idxv < edges[u].size(); ++idxv) {
      int v = edges[u][idxv];
      // count triangles
      if (v < u)
        continue;
      for (size_t idxw = idxv + 1; idxw < edges[u].size(); ++idxw) {
        // is {u, v, w} a triangle?
        int w = edges[u][idxw];
        assert(u != w && u != v);
        auto itr =
          std::lower_bound(edges[v].begin(), edges[v].end(), edges[u][idxw]);
        if (itr != edges[v].end() && *itr == w) {
          motifCounter[0] += 1;
        }
      }
    }
    motifCounter[1] += (edges[u].size() - 1) * (edges[u].size()) / 2;
  }
  motifCounter[1] -= 3 * motifCounter[0];
  return motifCounter;
}

// time complexity O(E^3/2)
vector<unsigned long long>
count3NodeCISForward(const vector<vector<int>> &edges, int vecsize) {
  vector<unsigned long long> motifCounter(2, 0);
  // similar to edge iterator
  // start from high degree node
  vector<int> order(edges.size());
  std::iota(order.begin(), order.end(), 0);
  std::sort(order.begin(), order.end(), [&](const int &lhs, const int &rhs) {
    return edges[lhs].size() > edges[rhs].size() ||
           (edges[lhs].size() == edges[rhs].size() && lhs < rhs);
  });
  vector<vector<int>> aux(edges.size()); // need aux memory
  for (size_t idx = 0; idx < order.size(); ++idx) {
    int s = order[idx];
    for (auto t : edges[s]) {
      int pos = 0;
      // find lower rank node (lower degree or higher id)
      if (edges[s].size() < edges[t].size() ||
          (edges[s].size() == edges[t].size() && s > t))
        continue;
      assert(s != t);
      // find common nodes edge list, mergesort like operation
      size_t idx1 = 0, idx2 = 0;
      size_t end1 = aux[s].size(), end2 = aux[t].size();
      while (idx1 != end1 && idx2 != end2) {
        if /**/ (aux[s][idx1] < aux[t][idx2])
          ++idx1;
        else if (aux[t][idx2] < aux[s][idx1])
          ++idx2;
        else {
          ++motifCounter[0];
          ++idx1;
          ++idx2;
        }
      }
      aux[t].push_back(idx);
    }
    motifCounter[1] += (edges[s].size() - 1) * (edges[s].size()) / 2;
  }
  motifCounter[1] -= 3 * motifCounter[0];
  return motifCounter;
}

// compute the local clustering coefficient
double localClusteringCoeff(const vector<vector<int>> &edges) {
  vector<int> counter(edges.size(), 0);
  long long sum = 0;
  // start from high degree node
  vector<int> order(edges.size());
  std::iota(order.begin(), order.end(), 0);
  std::sort(order.begin(), order.end(), [&](const int &lhs, const int &rhs) {
    return edges[lhs].size() > edges[rhs].size() ||
           (edges[lhs].size() == edges[rhs].size() && lhs < rhs);
  });
  vector<vector<int>> aux(edges.size()); // need aux memory
  for (size_t idx = 0; idx < order.size(); ++idx) {
    int s = order[idx];
    for (auto t : edges[s]) {
      int pos = 0;
      // find lower rank node (lower degree or higher id)
      if (edges[s].size() < edges[t].size() ||
          (edges[s].size() == edges[t].size() && s > t))
        continue;
      assert(s != t);
      // find common nodes edge list, mergesort like operation
      size_t idx1 = 0, idx2 = 0;
      size_t end1 = aux[s].size(), end2 = aux[t].size();
      while (idx1 != end1 && idx2 != end2) {
        if /**/ (aux[s][idx1] < aux[t][idx2])
          ++idx1;
        else if (aux[t][idx2] < aux[s][idx1])
          ++idx2;
        else {
          ++counter[s], ++counter[t], ++counter[order[aux[s][idx1]]];
          ++idx1;
          ++idx2;
        }
      }
      aux[t].push_back(idx);
    }
  }
  double localClusteringSum = 0.0;
  for (int i = 0; i < edges.size(); ++i) {
    sum += counter[i];
    if (edges[i].size() == 1)
      continue;
    localClusteringSum +=
      2.0 * counter[i] / edges[i].size() / (edges[i].size() - 1);
  }
  return localClusteringSum / edges.size();
}
