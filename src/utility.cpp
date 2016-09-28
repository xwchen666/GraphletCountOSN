#include "utility.h"
#include "process.h"

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <time.h>

#include <iostream>
#include <string>
#include <unordered_set>
#include <utility>
#include <vector>

using std::vector;
using std::pair;
using std::unordered_set;

std::string getGraphName(const std::string &fullpath) {
  size_t found = fullpath.find_last_of("/\\");
  return fullpath.substr(found + 1);
}

// truthfile: ground truth filename
// graphfile: input graph file
vector<long long> getMotifCounter(const char *truthfile,
                                  const char *graphfile) {
  std::string graphname = getGraphName(std::string(graphfile));
  std::ifstream infile(truthfile);
  std::string line;
  std::vector<long long> dataarea;

  while (std::getline(infile, line)) {
    std::istringstream iss(line);
    std::string fullpath;
    std::getline(iss, fullpath, '\t');
    std::string localgraphname = getGraphName(fullpath);
    if (localgraphname.compare(graphname) != 0)
      continue;
    long long data;
    while (iss >> data) {
      dataarea.push_back(data);
    }
    break;
  }
  infile.close();
  return dataarea;
}

//
// three nodes connected induced subgraph
// u: left node, v: right node
// result: first is the node id, second is the location indicator
// -1: left, 0: common, 1: right
int merge(const vector<int> &edgesOfU, const vector<int> &edgesOfV, const int u,
          const int v, const int w, vector<pair<int, int>> &result) {
  result.clear();
  auto first1 = edgesOfU.begin(), last1 = edgesOfU.end();
  auto first2 = edgesOfV.begin(), last2 = edgesOfV.end();
  while (true) {
    std::pair<int, int> ele = {0, -1};
    if (first1 == last1 && first2 == last2)
      break;
    if /**/ (first1 == last1)
      ele = {*first2++, 1};
    else if (first2 == last2)
      ele = {*first1++, -1};
    else {
      if (*first1 < *first2)
        ele = {*first1++, -1}; // 1 < 2
      else {
        if (*first1 == *first2)
          ele = {*first1++, 0}; // 1 == 2
        else
          ele = {*first2, 1}; // 1 > 2
        ++first2;
      }
    }
    if (result.empty() || ele.first != result.back().first) {
      if (ele.first != u && ele.first != v && ele.first != w)
        result.push_back(ele);
    }
  }
  return result.size();
}

// edgesOfU and edgesOfV are sorted
// ensure u is the left node and v is the right node
int intersection(const vector<int> &edgesOfU, const vector<int> &edgesOfV,
                 const int u, const int v, const int w,
                 vector<pair<int, int>> &result) {
  result.clear();
  auto first1 = edgesOfU.begin(), last1 = edgesOfU.end();
  auto first2 = edgesOfV.begin(), last2 = edgesOfV.end();
  while (true) {
    if /**/ (first1 == last1 || first2 == last2)
      break;
    if /**/ (*first1 < *first2)
      ++first1;
    else if (*first2 < *first1)
      ++first2;
    else {
      if (*first1 != u && *first1 != v && *first1 != w)
        result.push_back({*first1, 0});
      ++first1;
      ++first2;
    }
  }
  return result.size();
}

/* -------------------------------------------*/
/**
 * @brief populate the neighbors of current CIS
 *
 * @param edges edgelist of original graph
 */
/* -------------------------------------------*/
void Node3CIS::populate_nbs(const vector<vector<int>> &edges) {
  if (alreadyVisited)
    return;
  nbs.resize(3); // for node u, v, w (the neighborhood if we remove u, v, w
                 // correspondingly)
  int totalSize = 0;
  int u = current_tri[0], v = current_tri[1], w = current_tri[2];
  totalSize += merge(edges[v], edges[w], u, v, w, nbs[0]); // u
  totalSize += merge(edges[u], edges[v], u, v, w, nbs[2]); // w
  // is {u, v, w} a triangle
  if (isTriangle) {
    // {u, v, w} is a triangle
    totalSize += merge(edges[u], edges[w], u, v, w, nbs[1]); // v
  } else {
    // {u, v, w} is a line
    // ensure that we will get a connected CIS
    totalSize += intersection(edges[u], edges[w], u, v, w, nbs[1]); // v
  }
  this->degree = totalSize;
  alreadyVisited = true;
}

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
Node3CIS Node3CIS::getRandomNeighbor(const vector<vector<int>> &edges,
                                     vector<double> &concentration) {

  vector<int> next_tri = current_tri; // nodes for next CIS
  bool next_is_triangle = false;      // is next CIS a triangle?
  int index = rand() % this->degree;
  int j = 0;
  // locate
  for (j = 0; j < 3; ++j) {
    if ((index - (int)nbs[j].size()) < 0)
      break;
    else
      index -= nbs[j].size();
  }
  Motif_4_type motif_type;
  bool is_connected = false;
  // discussion @@ dazzle me @@
  // if we remove node v and u-v-w is not a triangle
  int oldnode, newnode = nbs[j][index].first;
  if (j == 1 && !this->isTriangle) {
    next_tri[j] = nbs[j][index].first;
    // is newnode connected to v
    oldnode = current_tri[1];
    auto itr =
      std::lower_bound(edges[oldnode].begin(), edges[oldnode].end(), newnode);
    if (itr != edges[oldnode].end() && *itr == newnode) {
      is_connected = true;
      motif_type = M4_CHORDALCYCLE;
    } else
      motif_type = M4_CYCLE;
  } else {
    // replace oldnode with newnode
    oldnode = current_tri[j];
    // is oldnode and newnode are connected
    auto itr =
      std::lower_bound(edges[oldnode].begin(), edges[oldnode].end(), newnode);
    if (itr != edges[oldnode].end() && *itr == newnode)
      is_connected = true;
    // remove the oldnode
    next_tri.erase(next_tri.begin() + j);
    if (nbs[j][index].second == -1) // adj to left node
      std::swap(next_tri[0], next_tri[1]);
    next_tri.push_back(nbs[j][index].first);
    next_is_triangle = nbs[j][index].second == 0;
    // which type motif ?
    if (this->isTriangle) { // u-v-w is a triangle
      if (next_is_triangle && is_connected)
        motif_type = M4_CLIQUE;
      if (next_is_triangle && !is_connected)
        motif_type = M4_CHORDALCYCLE;
      if (!next_is_triangle && is_connected)
        motif_type = M4_CHORDALCYCLE;
      if (!next_is_triangle && !is_connected)
        motif_type = M4_TAILEDTRIANGLE;
    } else { // u-v-w is line
      if (next_is_triangle && is_connected)
        motif_type = M4_CHORDALCYCLE;
      if (next_is_triangle && !is_connected)
        motif_type = M4_TAILEDTRIANGLE;
      if (!next_is_triangle) {
        bool leanatcenter =
          (oldnode == current_tri[0] && nbs[j][index].second == -1) ||
          (oldnode == current_tri[2] && nbs[j][index].second == 1);
        if (leanatcenter && !is_connected)
          motif_type = M3_STAR;
        if (leanatcenter && is_connected)
          motif_type = M4_TAILEDTRIANGLE;
        if (!leanatcenter && is_connected)
          motif_type = M4_CYCLE;
        if (!leanatcenter && !is_connected)
          motif_type = M4_PATH;
      }
    }
  }
  concentration[motif_type] += 1;
  return Node3CIS(next_tri, next_is_triangle);
}

int Node3CIS::getRandomNeighbor(const vector<vector<int>> &edges,
                                Node3CIS &nextCIS) {
  if (!alreadyVisited)
    populate_nbs(edges);

  vector<int> next_tri = current_tri; // nodes for next CIS
  bool next_is_triangle = false;      // is next CIS a triangle?
  int index = rand() % this->degree;
  int j = 0;
  // locate
  for (j = 0; j < 3; ++j) {
    if ((index - (int)nbs[j].size()) < 0)
      break;
    else
      index -= nbs[j].size();
  }
  // if we remove node v and u-v-w is not a triangle
  int oldnode, newnode = nbs[j][index].first;
  if (j == 1 && !this->isTriangle) {
    next_tri[j] = nbs[j][index].first;
  } else {
    // replace oldnode with newnode
    oldnode = current_tri[j];
    // remove the oldnode
    next_tri.erase(next_tri.begin() + j);
    if (nbs[j][index].second == -1) // adj to left node
      std::swap(next_tri[0], next_tri[1]);
    next_tri.push_back(nbs[j][index].first);
    next_is_triangle = nbs[j][index].second == 0;
  }
  nextCIS = Node3CIS(next_tri, next_is_triangle);
  return newnode;
}

/* -------------------------------------------*/
/**
 * @brief print all the neigbor CIS of current CIS (for debug)
 */
/* -------------------------------------------*/
void Node3CIS::printAllNeighbors(int k) {
  printf("<<< u = %d, v = %d, w = %d >>>\n", current_tri[0], current_tri[1],
         current_tri[2]);
  for (int j = 0; j < k; ++j) {
    vector<int> next_tri(current_tri);
    if (j != 1 || isTriangle) {
      next_tri.erase(next_tri.begin() + j);
      for (auto nb : nbs[j]) {
        int u = next_tri[0], v = next_tri[1];
        assert(u != v);
        if (nb.second == -1)
          std::swap(u, v);
        int w = nb.first;
        printf("(%d %d %d [%d])\n", u, v, w, nb.second == 0);
      }
    } else {
      for (auto nb : nbs[j]) {
        int u = next_tri[0], v = next_tri[2];
        printf("(%d %d %d [0])\n", u, nb.first, v);
      }
    }
  }
}

int Node3CIS::printNBS(int v) {
  int k = 0;
  for (; k < 3; ++k) {
    if (current_tri[k] == v)
      break;
  }
  for (auto t : nbs[k])
    printf("%d  ", t);
  printf("\n");
  return nbs[k].size();
}

// implementation of class member function
void Node3CISBIT::populate_nbs(const vector<vector<int>> &edges,
                               vector<int> &ind) {
  if (alreadyVisited)
    return;
  vector<int> node_set;
  nbs.resize(current_tri.size());
  nbs_bitSet.resize(current_tri.size());
  // marked the ind vector
  for (size_t i = 0; i < current_tri.size(); ++i) {
    for (auto e : edges[current_tri[i]]) {
      if (ind[e] == 0 && (std::find(current_tri.begin(), current_tri.end(),
                                    e) == current_tri.end())) {
        node_set.push_back(e);
      }
      ind[e] |= (1 << i);
    }
  }
  std::sort(node_set.begin(), node_set.end());
  // determine the motif type of current 3-nodes
  original_3_subgraph.resize(3);
  original_3_degree.resize(3);
  for (size_t i = 0; i < current_tri.size(); ++i) {
    original_3_subgraph[i] = ind[current_tri[i]];
    original_3_degree[i] = __builtin_popcount(ind[current_tri[i]]);
    ind[current_tri[i]] = 0;
  }
  this->isTriangle = determine_motif_type(hash_signature(original_3_degree), 3);
  assert(this->isTriangle != -1);
  // populate nbs
  for (size_t i = 0; i < current_tri.size(); ++i) {
    // replace the ith node
    for (auto t : node_set) {
      vector<int> local_subgraph_3_node(3, 0);
      vector<int> local_subgraph_3_degree(3, 0);
      // write bitSet of current nodes
      for (size_t j = 0, idx = 0; j < current_tri.size(); ++j) {
        if (j == i)
          continue;
        local_subgraph_3_node[idx] =
          original_3_subgraph[j] & ((~(1 << i)) & (0x7));
        local_subgraph_3_degree[idx] =
          __builtin_popcount(local_subgraph_3_node[idx]);
        ++idx;
      }
      // update bitSet
      for (size_t j = 0, idx = 0; j < current_tri.size(); ++j) {
        if (j == i)
          continue;
        if (((ind[t] >> j) & 1)) {
          ++local_subgraph_3_degree[2];
          ++local_subgraph_3_degree[idx];
        }
        ++idx;
      }
      local_subgraph_3_node[2] = ind[t];
      assert(local_subgraph_3_degree[2] ==
             __builtin_popcount((ind[t]) & ((~(1 << i)) & (0x7))));
      // check whether the new motif is connected
      if (determine_motif_type(hash_signature(local_subgraph_3_degree), 3) ==
          -1)
        continue;
      assert(ind[t] > 0 && ind[t] < 8);
      nbs[i].push_back(t);
      nbs_bitSet[i].push_back(ind[t]);
    }
  }
  this->alreadyVisited = true;
  // clear
  for (auto t : node_set)
    ind[t] = 0;
}

Node3CISBIT Node3CISBIT::getRandomNeighbor(vector<double> &concentration) {
  // get a random neighbor
  // determine the motif type plus this new node
  // return a 3 node cis
  int degree = getDegree();
  int ridx = rand() % degree;
  int nb = -1;
  size_t replace_node = 0;
  assert(nbs.size() == 3);
  for (; replace_node < nbs.size(); ++replace_node) {
    if (nbs[replace_node].size() == 0)
      continue;
    if (ridx >= nbs[replace_node].size()) {
      ridx -= nbs[replace_node].size();
    } else {
      break;
    }
  }
  nb = nbs[replace_node][ridx]; // get random neighbor
  vector<int> subgraph_4_node(4, 0);
  for (size_t i = 0; i < current_tri.size(); ++i) {
    subgraph_4_node[i] = original_3_subgraph[i];
  }
  subgraph_4_node[3] = nbs_bitSet[replace_node][ridx];
  for (int j = 0; j < 3; ++j) {
    if (((subgraph_4_node[3] >> j) & 1)) {
      subgraph_4_node[j] |= 1 << 3;
    }
  }
  vector<int> subgraph_4_degree(4, 0);
  for (size_t i = 0; i < subgraph_4_degree.size(); ++i) {
    subgraph_4_degree[i] = __builtin_popcount(subgraph_4_node[i]);
  }
  int type = determine_motif_type(hash_signature(subgraph_4_degree), 4);
  concentration[type] += 1;
  vector<int> next_nodes(3, 0);
  for (size_t i = 0, idx = 0; i < current_tri.size(); ++i) {
    if (i == replace_node)
      continue;
    next_nodes[idx++] = current_tri[i];
  }
  next_nodes[2] = nb;
  return Node3CISBIT(next_nodes);
}

void Node3CISBIT::printAllNeighbors(int k) {
  printf("<<< u = %d, v = %d, w = %d >>>\n", current_tri[0], current_tri[1],
         current_tri[2]);
  for (int j = 0; j < k; ++j) {
    vector<int> nodes(2, 0);
    for (int i = 0, idx = 0; i < 3; ++i) {
      if (i == j)
        continue;
      nodes[idx++] = current_tri[i];
    }
    for (auto nb : nbs[j]) {
      printf("(%d %d %d)\n", nodes[0], nodes[1], nb);
    }
  }
}

// implementation of class member function
void Node4CIS::populate_nbs(const vector<vector<int>> &edges,
                            vector<int> &ind) {
  if (alreadyVisited)
    return;
  vector<int> node_set;
  nbs.resize(current_nodes.size());
  nbs_bitSet.resize(current_nodes.size());
  // marked the ind vector
  for (size_t i = 0; i < current_nodes.size(); ++i) {
    for (auto e : edges[current_nodes[i]]) {
      if (ind[e] == 0 && (std::find(current_nodes.begin(), current_nodes.end(),
                                    e) == current_nodes.end())) {
        node_set.push_back(e);
      }
      ind[e] |= (1 << i);
    }
  }
  std::sort(node_set.begin(), node_set.end());
  // determine the motif type of current 4-nodes
  original_4_subgraph.resize(4);
  original_4_degree.resize(4);
  for (size_t i = 0; i < current_nodes.size(); ++i) {
    original_4_subgraph[i] = ind[current_nodes[i]];
    original_4_degree[i] = __builtin_popcount(ind[current_nodes[i]]);
    ind[current_nodes[i]] = 0;
  }
  this->type = determine_motif_type(hash_signature(original_4_degree), 4);
  assert(this->type != -1);
  // populate nbs
  for (size_t i = 0; i < current_nodes.size(); ++i) {
    // replace the ith node
    for (auto t : node_set) {
      vector<int> local_subgraph_4_node(4, 0);
      vector<int> local_subgraph_4_degree(4, 0);
      // write bitSet of current nodes
      for (size_t j = 0, idx = 0; j < current_nodes.size(); ++j) {
        if (j == i)
          continue;
        local_subgraph_4_node[idx] =
          original_4_subgraph[j] & ((~(1 << i)) & (0xF));
        local_subgraph_4_degree[idx] =
          __builtin_popcount(local_subgraph_4_node[idx]);
        ++idx;
      }
      // update bitSet
      for (size_t j = 0, idx = 0; j < current_nodes.size(); ++j) {
        if (j == i)
          continue;
        if (((ind[t] >> j) & 1)) {
          ++local_subgraph_4_degree[3];
          ++local_subgraph_4_degree[idx];
        }
        ++idx;
      }
      // check whether the new motif is connected
      if (determine_motif_type(hash_signature(local_subgraph_4_degree), 4) ==
          -1)
        continue;
      assert(ind[t] > 0 && ind[t] < 16);
      nbs[i].push_back(t);
      nbs_bitSet[i].push_back(ind[t]);
    }
  }
  this->alreadyVisited = true;
  // clear
  for (auto t : node_set)
    ind[t] = 0;
}

/* -------------------------------------------*/
/**
 * @brief get random neighbor of current CIS
 *
 * @param concentration
 *
 * @return  a new neighborhood 4-node cis
 */
/* -------------------------------------------*/
Node4CIS Node4CIS::getRandomNeighbor(vector<double> &concentration) {
  // get a random neighbor
  // determine the motif type plus this new node
  // return a 4 node cis
  int degree = getDegree();
  int ridx = rand() % degree;
  int nb = -1;
  size_t replace_node = 0;
  for (; replace_node < nbs.size(); ++replace_node) {
    if (nbs[replace_node].size() == 0)
      continue;
    if (ridx >= nbs[replace_node].size()) {
      ridx -= nbs[replace_node].size();
    } else {
      break;
    }
  }
  nb = nbs[replace_node][ridx]; // get random neighbor
  vector<int> subgraph_5_node(5, 0);
  for (size_t i = 0; i < current_nodes.size(); ++i) {
    subgraph_5_node[i] = original_4_subgraph[i];
  }
  subgraph_5_node[4] = nbs_bitSet[replace_node][ridx];
  for (int j = 0; j < 4; ++j) {
    if (((subgraph_5_node[4] >> j) & 1)) {
      subgraph_5_node[j] |= 1 << 4;
    }
  }
  int type = decide_motif_type<5>(subgraph_5_node);
  assert(type != -1);
  concentration[type] += 1;
  vector<int> next_nodes(4, 0);
  for (size_t i = 0, idx = 0; i < current_nodes.size(); ++i) {
    if (i == replace_node)
      continue;
    next_nodes[idx++] = current_nodes[i];
  }
  next_nodes[3] = nb;
  return Node4CIS(next_nodes);
}
