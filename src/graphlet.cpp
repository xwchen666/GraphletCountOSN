#include "graphlet.h"
#include <algorithm>
#include <iostream>
#include <sstream>
using std::cout;
using std::endl;

// use the same format and variable name as PGD to avoid confusion
void graphlet_core::reformat(const vector<vector<int>> &origin_edges) {
  edges.clear();
  vertices.clear();
  vertices.push_back(edges.size());
  for (size_t i = 0; i < origin_edges.size(); ++i) {
    for (auto e : origin_edges[i]) {
      edges.push_back(e);
    }
    vertices.push_back(edges.size());
  }
}

// compute the degree for each vertices
// return max_degree
void graphlet_core::vertex_degrees() {
  long long n = vertices.size() - 1;
  degree.resize(n);
  printf("n = %lu\n", n);
  int max_degree_tmp = vertices[1] - vertices[0];
#pragma omp parallel for schedule(schedule_type,                               \
                                  block_size) reduction(max : max_degree_tmp)
  for (long long v = 0; v < n; ++v) {
    degree[v] = vertices[v + 1] - vertices[v];
    if (max_degree_tmp < degree[v])
      max_degree_tmp = degree[v];
  }
  max_degree = max_degree_tmp;
  avg_degree = 1.0 * vertices[n - 1] / n;
}

// initialize
void graphlet_core::initialize() {
  max_degree = 0;
  avg_degree = 0;
  block_size = 64;
}

graphlet_core::~graphlet_core() {}

// creare edge list arrays
void graphlet_core::create_edge_list_arrays() {
  long long n = vertices.size() - 1;
  long long m = edges.size();
  long long unique_edges = m / 2;
  e_v.reserve(unique_edges + 1);
  e_u.reserve(unique_edges + 1);
  printf("n = %lu, m = %lu\n", n, m);
  for (long long v = 0; v < n; v++) {
    for (long long j = vertices[v]; j < (long long)vertices[v + 1]; j++) {
      long long u = edges[j];
      if (v < u) {
        e_v.push_back(v);
        e_u.push_back(u);
      }
    }
  }
}

// sort edges according to degree
void graphlet_core::sort_edges() {
  int m = e_v.size();
  E_ordered.reserve(m);
  for (int e = 0; e < m; ++e) {
    long long v = e_v[e], u = e_u[e];
    long long val = degree[u] + degree[v];
    E_ordered.push_back(Vertex(e, val));
  }
  // large_to_small
  std::sort(E_ordered.begin(), E_ordered.end(), decr_bound);
}

inline void graphlet_core::mark_neighbors(long long &v, long long &u,
                                          vector<long long> &ind) {
  for (long long i = vertices[v]; i < (long long)vertices[v + 1]; ++i) {
    if (edges[i] == u)
      continue;
    ind[edges[i]] = 1;
  }
}

inline void graphlet_core::reset_perfect_hash(long long &v,
                                              vector<long long> &ind) {
  for (long long i = vertices[v]; i < (long long)vertices[v + 1]; ++i) {
    ind[edges[i]] = 0;
  }
}

// compute number of triangles and 2-stars for vertex u
inline void graphlet_core::triangles_and_wedges(
  long long &v, long long &u, vector<long long> &T_vu,
  unsigned long long &tri_count, vector<long long> &W_u,
  unsigned long long &w_local_count, vector<long long> &ind) {
  for (long long j = vertices[u]; j < (long long)vertices[u + 1]; ++j) {
    long long w = edges[j];
    if (w == v) {
      continue;
    }
    if (ind[w] == 1) {
      ind[w] = 3; // triangle
      T_vu[tri_count++] = w;
    } else {
      W_u[w_local_count++] = w;
      ind[w] = 2; // star for vertex u
    }
  }
}

// compute the number of cycles from star vertices of node u
inline void graphlet_core::cycle(unsigned long long &w_local_count,
                                 vector<long long> &W_u,
                                 unsigned long long &cycle4_count, long long &v,
                                 vector<long long> &ind) {
  for (unsigned long long j = 0; j < w_local_count; ++j) {
    long long w = W_u[j];
    for (long long i = vertices[w]; i < (long long)vertices[w + 1]; ++i) {
      // ind[edges[i]] == 1, star for vertex v
      if (ind[edges[i]] == 1)
        ++cycle4_count;
    }
    W_u[j] = 0;
  }
}

// Compute the number of cliques centered at edge u-v
inline void graphlet_core::clique(unsigned long long &tri_count,
                                  vector<long long> &T_vu,
                                  unsigned long long &clique4_count,
                                  long long &v, vector<long long> &ind) {
  for (unsigned long long tr_i = 0; tr_i < tri_count; ++tr_i) {
    long long w = T_vu[tr_i];
    for (long long i = vertices[w]; i < (long long)vertices[w + 1]; ++i) {
      if (ind[edges[i]] == 3)
        ++clique4_count;
    }
    ind[w] = 0;
    T_vu[tr_i] = 0;
  }
}

//  Reset the counts of all graphlet motifs to zero
inline void graphlet_core::reset_graphlet_counts() {
  // connected k=3 motifs
  total_3_tris = 0.0;
  total_2_star = 0.0;
  // connected k=4 motifs
  total_4_clique = 0.0;
  total_4_chordcycle = 0.0;
  total_4_tailed_tris = 0.0;
  total_4_cycle = 0.0;
  total_3_star = 0.0;
  total_4_path = 0.0;
}

// Check correctness of graphlet counts
// thread_n_count is an array containing the counts for each graphlet motif
inline bool graphlet_core::test_graphlet_counts(
  vector<vector<unsigned long long>> &thread_n_count) {
  unsigned long long ver_n1 = 3 * total_3_star + 3 * total_4_tailed_tris +
                              4 * total_4_cycle + total_4_path +
                              5 * total_4_chordcycle + 6 * total_4_clique;
  unsigned long long ver_n2 =
    3 * total_3_star + total_4_tailed_tris + 4 * total_4_cycle + total_4_path;
  unsigned long long ver_n3 = total_4_chordcycle + 6 * total_4_clique;
  unsigned long long ver_n4 = 2 * total_4_tailed_tris + 4 * total_4_chordcycle;
  unsigned long long ver_n5 = total_4_path + 4 * total_4_cycle;
  unsigned long long ver_n6 = 3 * total_3_star + total_4_tailed_tris;
  int total_verify = 0;
  total_verify = ((ver_n1 - thread_n_count[0][1]) > 0) +
                 ((ver_n2 - thread_n_count[0][2]) > 0) +
                 ((ver_n3 - thread_n_count[0][3]) > 0);
  total_verify += ((ver_n4 - thread_n_count[0][4]) > 0) +
                  ((ver_n5 - thread_n_count[0][5]) > 0) +
                  ((ver_n6 - thread_n_count[0][6]) > 0);
  if (total_verify != 0) {
    cout << "\tver_n1 = " << ver_n1 - thread_n_count[0][1] << std::endl;
    cout << "\tver_n2 = " << ver_n2 - thread_n_count[0][2] << std::endl;
    cout << "\tver_n3 = " << ver_n3 - thread_n_count[0][3] << std::endl;
    cout << "\tver_n4 = " << ver_n4 - thread_n_count[0][4] << std::endl;
    cout << "\tver_n5 = " << ver_n5 - thread_n_count[0][5] << std::endl;
    cout << "\tver_n6 = " << ver_n6 - thread_n_count[0][6] << std::endl;
    return false;
  }
  return true;
}

/**
 * @brief Direct computation of numerous graphlet counts by leveraging
 * relationships
 * and solutions to a number of graphlet equations.
 * Resulting counts are obtained in O(1) time.
 *
 * @param thread_n_count
 * @param n_count
 * @param thread_tmp_triangles
 * @param thread_tmp_3_star
 * @param deg_v
 * @param deg_u
 * @param tri_count
 * @param w_local_count
 * @param m
 * @param n
 */
inline void graphlet_core::solve_graphlet_equations(
  vector<vector<unsigned long long>> &thread_n_count,
  vector<unsigned long long> &n_count,
  vector<unsigned long long> &thread_tmp_triangles,
  vector<unsigned long long> &thread_tmp_3_star, long long &deg_v,
  long long &deg_u, unsigned long long &tri_count,
  unsigned long long &w_local_count, long long &m, long long &n) {

  unsigned long long star3_count = 0;
  thread_tmp_triangles[omp_get_thread_num()] += tri_count;

  star3_count = deg_v - tri_count - 1;
  star3_count = star3_count + deg_u - tri_count - 1;
  thread_tmp_3_star[omp_get_thread_num()] += star3_count;

  n_count[1] = (tri_count + star3_count) * (tri_count + star3_count - 1) / 2.0;
  n_count[2] = (star3_count * (star3_count - 1) / 2.0);
  n_count[3] = (tri_count * (tri_count - 1) / 2.0);
  n_count[4] = (tri_count * star3_count);
  n_count[5] = (deg_v - tri_count - 1) * (deg_u - tri_count - 1);
  n_count[6] = (deg_v - tri_count - 1) * (deg_v - tri_count - 2) / 2;
  n_count[6] =
    n_count[6] + ((deg_u - tri_count - 1) * (deg_u - tri_count - 2) / 2);

  thread_n_count[omp_get_thread_num()][1] += n_count[1];
  thread_n_count[omp_get_thread_num()][2] += n_count[2];
  thread_n_count[omp_get_thread_num()][3] += n_count[3];
  thread_n_count[omp_get_thread_num()][4] += n_count[4];
  thread_n_count[omp_get_thread_num()][5] += n_count[5];
  thread_n_count[omp_get_thread_num()][6] += n_count[6];
}

void graphlet_core::graphlet_decomposition(int max_num_workers) {
  omp_set_num_threads(max_num_workers);
  long long e, v, u, w, r, i, j, tr_i, n = num_vertices(), m = num_edges();
  vector<long long> T_vu(max_degree + 1, 0), W_u(max_degree + 1, 0);
  vector<unsigned long long> thread_tmp_triangles(max_num_workers, 0);
  vector<unsigned long long> thread_tmp_3_star(max_num_workers, 0);
  vector<unsigned long long> thread_tmp_cycles(max_num_workers, 0);
  vector<unsigned long long> thread_tmp_cliques(max_num_workers, 0);
  vector<vector<unsigned long long>> thread_n_count(
    max_num_workers, vector<unsigned long long>(7, 0));
  vector<unsigned long long> n_count(7, 0);

  if (E_ordered.size() == 0) {
    sort_edges();
  }
  double sec = tic();

  vector<long long> ind(n, 0);
  sec = tic();
#pragma omp parallel for schedule(schedule_type, block_size) firstprivate(     \
  ind, T_vu, W_u, n_count) private(v, u, w, r, i, j, tr_i, e)
  for (long long e = 0; e < (long long)E_ordered.size(); ++e) {
    int edge_id = E_ordered[e].get_id();
    long long v = e_v[edge_id], u = e_u[edge_id];
    long long deg_v = vertices[v + 1] - vertices[v];
    long long deg_u = vertices[u + 1] - vertices[u];
    unsigned long long w_local_count = 0, tri_count = 0, clique4_count = 0,
                       cycle4_count = 0;
    mark_neighbors(v, u, ind);
    triangles_and_wedges(v, u, T_vu, tri_count, W_u, w_local_count, ind);
    solve_graphlet_equations(thread_n_count, n_count, thread_tmp_triangles,
                             thread_tmp_3_star, deg_v, deg_u, tri_count,
                             w_local_count, m, n);
    cycle(w_local_count, W_u, cycle4_count, v, ind);
    clique(tri_count, T_vu, clique4_count, v, ind);
    thread_tmp_cliques[omp_get_thread_num()] += clique4_count;
    thread_tmp_cycles[omp_get_thread_num()] += cycle4_count;
    reset_perfect_hash(v, ind);
  }
  toc(sec);
  ind.clear();

  reset_graphlet_counts();

  for (int tid = 1; tid < max_num_workers; ++tid) {
    thread_tmp_triangles[0] += thread_tmp_triangles[tid];
    thread_tmp_3_star[0] += thread_tmp_3_star[tid];
    thread_tmp_cliques[0] += thread_tmp_cliques[tid];
    thread_tmp_cycles[0] += thread_tmp_cycles[tid];

    thread_n_count[0][1] += thread_n_count[tid][1];
    thread_n_count[0][2] += thread_n_count[tid][2];
    thread_n_count[0][3] += thread_n_count[tid][3];
    thread_n_count[0][5] += thread_n_count[tid][5];
    thread_n_count[0][4] += thread_n_count[tid][4];
    thread_n_count[0][6] += thread_n_count[tid][6];
  }
  total_3_tris = thread_tmp_triangles[0] / 3.0;
  total_4_clique = thread_tmp_cliques[0] / 6.0;
  total_4_chordcycle = thread_n_count[0][3] - (6 * total_4_clique);
  total_4_cycle = thread_tmp_cycles[0] / 4.0;
  total_4_path = thread_n_count[0][5] - (4 * total_4_cycle);
  total_4_tailed_tris = (thread_n_count[0][4] - (4 * total_4_chordcycle)) / 2.0;
  total_3_star = (thread_n_count[0][6] - total_4_tailed_tris) / 3.0;
  test_graphlet_counts(thread_n_count);
}

string graphlet_core::compute_connected_GFD() {
  vector<unsigned long long> G4;
  G4.push_back(total_4_clique);
  G4.push_back(total_4_chordcycle);
  G4.push_back(total_4_tailed_tris);
  G4.push_back(total_4_cycle);
  G4.push_back(total_3_star);
  G4.push_back(total_4_path);

  const char *G4_name_tmp[] = {"4-clique     ", "4-chordal-cycle",
                               "4-tailed-tri",  "4-cycle       ",
                               "3-star       ", "4-path       "};
  vector<string> G4_names(G4_name_tmp, G4_name_tmp + G4.size());
  std::ostringstream os;
  unsigned long long sum = 0;
  for (size_t i = 0; i < G4.size(); i++) {
    sum += G4[i];
  }
  vector<double> gfd_tmp(G4.size(), 0);
  for (size_t i = 0; i < G4.size(); i++) {
    if (sum > 0)
      gfd_tmp[i] = ((double)G4[i] / (double)sum);
    else
      gfd_tmp[i] = 0;
    os << G4_names[i] << "\t" << gfd_tmp[i] << "\n";
  }
  os << "\n";
  connected_GFD = gfd_tmp;
  return os.str();
}

void graphlet_core::print_connected_GFD() {
  print_line(80);
  std::cout << "Connected Motif Graphlet Frequency Distribution (GFD)\n";
  print_line(80);
  std::cout << compute_connected_GFD() << std::endl;
}
