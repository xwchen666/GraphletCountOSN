#ifndef GRAPHLET_H_
#define GRAPHLET_H_

// parallel
#define schedule_type dynamic

#ifdef _OPENMP
#include <omp.h>
#else
int omp_get_max_threads() { return 1; }
void omp_set_num_threads(int) {}
int omp_get_thread_num() { return 0; }
#endif

// normal headers
#include "process.h"
#include <cstddef>
#include <iostream>
#include <limits>
#include <stdlib.h>
#include <string>
#include <sys/time.h>
#include <sys/timeb.h>
#include <sys/types.h>
#include <unistd.h>
#include <vector>
using std::vector;
using std::string;

// utility
inline double get_time() {
  struct timeval t;
  gettimeofday(&t, 0);
  return (t.tv_sec * 1.0 + t.tv_usec / 1000000.0);
}

inline void print_line(int n = 80, string sym = "-") {
  for (int i = 0; i < n; ++i)
    std::cout << sym;
  std::cout << std::endl;
}

inline double tic() { return get_time(); }

inline void toc(double &start) { start = get_time() - start; }

// vertex class
class Vertex {
private:
  int id, b;

public:
  Vertex(int vertex_id, int bound) : id(vertex_id), b(bound){};

  void set_id(int vid) { id = vid; }
  int get_id() { return id; }

  void set_bound(int value) { b = value; }
  int get_bound() { return b; }
};

/**
 * @brief Order from largest to smallest value
 * Note that ties are broken by vertex id
 *
 * @param v is a vertex object containing a vertex id and its value
 * @param u is a vertex object containing a vertex id and its value
 * @return Returns true if the value of v is larger than the value of u
 *
 */
static bool decr_bound(Vertex v, Vertex u) {
  return (v.get_bound() > u.get_bound() ||
          (v.get_bound() == u.get_bound() && v.get_id() > u.get_id()));
}
/**
 * @brief Order from smallest to largest value
 * Note that ties are broken by vertex id
 *
 * @param v is a vertex object containing a vertex id and its value
 * @param u is a vertex object containing a vertex id and its value
 * @return Returns true if the value of v is smaller than the value of u
 */
static bool incr_bound(Vertex v, Vertex u) {
  return (v.get_bound() < u.get_bound() ||
          (v.get_bound() == u.get_bound() && v.get_id() > u.get_id()));
};

// graphlet_core class
class graphlet_core {
private:
  void reformat(const vector<vector<int>> &origin_edges);

public:
  // number of jobs each worker is assigned (at a time).
  // Must be a positive integer.
  int block_size;
  vector<int> edges;

  // stores the position of the neighbors in edges
  vector<size_t> vertices;

  // Vertex degree stats
  vector<int> degree;
  int max_degree, avg_degree;

  void initialize();
  ~graphlet_core();
  graphlet_core() {}
  graphlet_core(const vector<vector<int>> &origin_edges) {
    initialize();              // set variables to zero
    reformat(origin_edges);    // reformat
    vertex_degrees();          // vertex degrees
    create_edge_list_arrays(); // create sorted edge list
    sort_edges();
  }

  void create_edge_list_arrays();

  int num_vertices() { return vertices.size() - 1; }

  int num_edges() { return edges.size() / 2; }

  void vertex_degrees();

  int get_max_degree() { return max_degree; }

  double get_avg_degree() { return avg_degree; }

  vector<Vertex> E_ordered;
  void sort_edges();
  vector<long long> e_v, e_u;

  unsigned long long total_3_tris; // triangle
  unsigned long long total_2_star; // line

  unsigned long long total_4_clique;
  unsigned long long total_4_chordcycle;
  unsigned long long total_4_tailed_tris;
  unsigned long long total_4_cycle;
  unsigned long long total_3_star;
  unsigned long long total_4_path;

  // connected (k=3) motifs
  /** @brief triangle counts (per edge). */
  vector<unsigned long long> tri;

  inline unsigned long long get_2_star(long long &e, long long &deg_v,
                                       long long &deg_u) {
    return (deg_v - tri[e]) + (deg_u - tri[e]) - 2;
  }

  vector<string> get_graphlet_size4_names() {
    vector<string> g4_names;
    g4_names.reserve(6);
    g4_names.push_back("4-clique");
    g4_names.push_back("4-chordal-cycle");
    g4_names.push_back("4-tailed-tri");
    g4_names.push_back("4-cycle");
    g4_names.push_back("3-star");
    g4_names.push_back("4-path");
    return g4_names;
  }

  /**
   * @brief This function gets the stats of the k=4 graphlets.
   *
   * @return vector of the k=4 graphlet stats
   */
  vector<unsigned long long> get_graphlet_size4_stats() {
    vector<unsigned long long> g4_stats;
    g4_stats.reserve(6);
    g4_stats.push_back(total_4_clique);
    g4_stats.push_back(total_4_chordcycle);
    g4_stats.push_back(total_4_tailed_tris);
    g4_stats.push_back(total_4_cycle);
    g4_stats.push_back(total_3_star);
    g4_stats.push_back(total_4_path);
    return g4_stats;
  }

  vector<double> connected_GFD;

  string compute_connected_GFD();

  void print_connected_GFD();

  void solve_graphlet_equations(
    vector<vector<unsigned long long>> &thread_n_count,
    vector<unsigned long long> &n_count,
    vector<unsigned long long> &thread_tmp_triangles,
    vector<unsigned long long> &thread_tmp_3_star, long long &deg_v,
    long long &deg_u, unsigned long long &tri_count,
    unsigned long long &w_local_count, long long &m, long long &n);

  /** @brief CSC SPARSE GRAPH REP. ONLY */
  void triangles_and_wedges(long long &v, long long &u, vector<long long> &T_vu,
                            unsigned long long &tri_count,
                            vector<long long> &W_u,
                            unsigned long long &w_local_count,
                            vector<long long> &ind);

  void clique(unsigned long long &tri_count, vector<long long> &T_vu,
              unsigned long long &clique4_count, long long &v,
              vector<long long> &ind);

  void cycle(unsigned long long &w_local_count, vector<long long> &W_u,
             unsigned long long &cycle4_count, long long &v,
             vector<long long> &ind);

  void mark_neighbors(long long &v, long long &u, vector<long long> &ind);

  void reset_perfect_hash(long long &v, vector<long long> &ind);

  void reset_graphlet_counts();

  bool test_graphlet_counts(vector<vector<unsigned long long>> &thread_n_count);

  void graphlet_decomposition(int max_num_workers);
};
#endif
