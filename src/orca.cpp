#include "orca.h"

/** count graphlets on max 4 nodes */
void orca::count4() {
  clock_t startTime, endTime;
  startTime = clock();
  clock_t startTime_all, endTime_all;
  startTime_all = startTime;
  int frac, frac_prev;

  // precompute triangles that span over edges
  // this step is not neccessary
  printf("stage 1 - precomputing common nodes\n");
  int *tri = (int *)calloc(m, sizeof(int));
  frac_prev = -1;
#pragma omp parallel for schedule(guided, 1)
  for (int i = 0; i < m; ++i) {
    int x = edges[i].a, y = edges[i].b;
    for (int xi = 0, yi = 0; xi < deg[x] && yi < deg[y];) {
      // mergesort-like scan to compute number of triangles
      if (adj[x][xi] == adj[y][yi]) {
        ++tri[i], ++xi, ++yi;
      } else if (adj[x][xi] < adj[y][yi]) {
        ++xi;
      } else {
        ++yi;
      }
    }
  }
  endTime = clock();
  printf("%.2f\n", (double)(endTime - startTime) / CLOCKS_PER_SEC);
  startTime = endTime;

  // count full graphlets
  printf("stage 2 - counting full graphlets\n");
  int64 *C4 = (int64 *)calloc(n, sizeof(int64));
  int *neigh = (int *)malloc(n * sizeof(int));
  int nn;
  frac_prev = -1;
  for (int x = 0; x < n; ++x) {
    frac = 100LL * x / n;
    if (frac != frac_prev) {
      printf("%d%%\r", frac);
      fflush(stdout);
      frac_prev = frac;
    }
    for (int y : adj[x]) {
      if (y >= x)
        break;
      nn = 0;
      for (int z : adj[y]) {
        if (z >= y)
          break;
        if ((this->*adjacent)(x, z) == 0)
          continue;
        neigh[nn++] = z;
      }
      for (int i = 0; i < nn; ++i) {
        int z = neigh[i];
        for (int j = i + 1; j < nn; ++j) {
          int zz = neigh[j];
          if ((this->*adjacent)(z, zz)) {
            ++C4[x], ++C4[y], ++C4[z], ++C4[zz];
          }
        }
      }
    }
  }
  endTime = clock();
  printf("%.2f\n", (double)(endTime - startTime) / CLOCKS_PER_SEC);
  startTime = endTime;

  // set up a system of equations relating orbits for every node
  printf("stage 3 - building systems of equations\n");
  int *common = (int *)calloc(n, sizeof(int));
  int *common_list = (int *)malloc(n * sizeof(int));
  int nc = 0;
  frac_prev = -1;
  for (int x = 0; x < n; ++x) {
    frac = 100LL * x / n;
    if (frac != frac_prev) {
      printf("%d%%\r", frac);
      fflush(stdout);
      frac_prev = frac;
    }

    int64 f_12_14 = 0, f_10_13 = 0;
    int64 f_13_14 = 0, f_11_13 = 0;
    int64 f_7_11 = 0, f_5_8 = 0;
    int64 f_6_9 = 0, f_9_12 = 0, f_4_8 = 0, f_8_12 = 0;
    int64 f_14 = C4[x];

    for (int i = 0; i < nc; ++i)
      common[common_list[i]] = 0;
    nc = 0;

    orbit[x][0] = deg[x];
    // x - middle node
    for (int nx1 = 0; nx1 < deg[x]; ++nx1) {
      int y = inc[x][nx1].first, ey = inc[x][nx1].second;
      for (int ny = 0; ny < deg[y]; ++ny) {
        int z = inc[y][ny].first, ez = inc[y][ny].second;
        if ((this->*adjacent)(x, z)) { // triangle
          if (z < y) {
            f_12_14 += tri[ez] - 1;
            f_10_13 += (deg[y] - 1 - tri[ez]) + (deg[z] - 1 - tri[ez]);
          }
        } else {
          if (common[z] == 0)
            common_list[nc++] = z;
          common[z]++;
        }
      }
      for (int nx2 = nx1 + 1; nx2 < deg[x]; ++nx2) {
        int z = inc[x][nx2].first, ez = inc[x][nx2].second;
        if ((this->*adjacent)(y, z)) { // triangle
          ++orbit[x][3];
          f_13_14 += (tri[ey] - 1) + (tri[ez] - 1);
          f_11_13 += (deg[x] - 1 - tri[ey]) + (deg[x] - 1 - tri[ez]);
        } else { // path
          orbit[x][2]++;
          f_7_11 += (deg[x] - 1 - tri[ey] - 1) + (deg[x] - 1 - tri[ez] - 1);
          f_5_8 += (deg[y] - 1 - tri[ey]) + (deg[z] - 1 - tri[ez]);
        }
      }
    }
    // x - side node
    for (int nx1 = 0; nx1 < deg[x]; ++nx1) {
      int y = inc[x][nx1].first, ey = inc[x][nx1].second;
      for (int ny = 0; ny < deg[y]; ++ny) {
        int z = inc[y][ny].first, ez = inc[y][ny].second;
        if (x == z)
          continue;
        if (!(this->*adjacent)(x, z)) { // path
          orbit[x][1]++;
          f_6_9 += (deg[y] - 1 - tri[ey] - 1);
          f_9_12 += tri[ez];
          f_4_8 += (deg[z] - 1 - tri[ez]);
          f_8_12 += (common[z] - 1);
        }
      }
    }

    // solve system of equations
    orbit[x][14] = (f_14);
    orbit[x][13] = (f_13_14 - 6 * f_14) / 2;
    orbit[x][12] = (f_12_14 - 3 * f_14);
    orbit[x][11] = (f_11_13 - f_13_14 + 6 * f_14) / 2;
    orbit[x][10] = (f_10_13 - f_13_14 + 6 * f_14);
    orbit[x][9] = (f_9_12 - 2 * f_12_14 + 6 * f_14) / 2;
    orbit[x][8] = (f_8_12 - 2 * f_12_14 + 6 * f_14) / 2;
    orbit[x][7] = (f_13_14 + f_7_11 - f_11_13 - 6 * f_14) / 6;
    orbit[x][6] = (2 * f_12_14 + f_6_9 - f_9_12 - 6 * f_14) / 2;
    orbit[x][5] = (2 * f_12_14 + f_5_8 - f_8_12 - 6 * f_14);
    orbit[x][4] = (2 * f_12_14 + f_4_8 - f_8_12 - 6 * f_14);
  }

  endTime = clock();
  printf("%.2f\n", (double)(endTime - startTime) / CLOCKS_PER_SEC);

  endTime_all = endTime;
  printf("total: %.2f\n",
         (double)(endTime_all - startTime_all) / CLOCKS_PER_SEC);
}

/** count graphlets on max 5 nodes */
void orca::count5() {
  clock_t startTime, endTime;
  startTime = clock();
  clock_t startTime_all, endTime_all;
  startTime_all = startTime;

  // precompute common nodes
  printf("stage 1 - precomputing common nodes\n");
  for (int x = 0; x < n; ++x) {
    for (int n1 = 0; n1 < deg[x]; ++n1) {
      int a = adj[x][n1];
      for (int n2 = n1 + 1; n2 < deg[x]; ++n2) {
        int b = adj[x][n2];
        PAIR ab = PAIR(a, b);
        ++common2[ab];
        for (int n3 = n2 + 1; n3 < deg[x]; ++n3) {
          int c = adj[x][n3];
          int st = (this->*adjacent)(a, b) + (this->*adjacent)(a, c) +
                   (this->*adjacent)(b, c);
          if (st < 2)
            continue;
          TRIPLE abc = TRIPLE(a, b, c);
          ++common3[abc];
        }
      }
    }
  }
  // precompute triangles that span over edges
  int *tri = (int *)calloc(m, sizeof(int));
  for (int i = 0; i < m; ++i) {
    int x = edges[i].a, y = edges[i].b;
    tri[i] = common2[PAIR(x, y)];
  }
  endTime = clock();
  printf("%.2f sec\n", (double)(endTime - startTime) / CLOCKS_PER_SEC);
  startTime = endTime;

  // count full graphlets
  printf("stage 2 - counting full graphlets\n");
  int64 *C5 = (int64 *)calloc(n, sizeof(int64));
  int *neigh = (int *)malloc(n * sizeof(int)), nn;
  int *neigh2 = (int *)malloc(n * sizeof(int)), nn2;
  for (int x = 0; x < n; ++x) {
    for (int y : adj[x]) {
      if (y >= x)
        break;
      nn = 0;
      for (int z : adj[y]) {
        if (z >= y)
          break;
        if ((this->*adjacent)(x, z)) {
          neigh[nn++] = z;
        }
      }
      // x > y > z
      for (int i = 0; i < nn; i++) {
        int z = neigh[i];
        nn2 = 0;
        for (int j = i + 1; j < nn; j++) {
          int zz = neigh[j];
          if ((this->*adjacent)(z, zz)) {
            neigh2[nn2++] = zz;
          }
        }
        for (int i2 = 0; i2 < nn2; i2++) {
          int zz = neigh2[i2];
          for (int j2 = i2 + 1; j2 < nn2; j2++) {
            int zzz = neigh2[j2];
            if ((this->*adjacent)(zz, zzz)) {
              C5[x]++;
              C5[y]++;
              C5[z]++;
              C5[zz]++;
              C5[zzz]++;
            }
          }
        }
      }
    }
  }
  endTime = clock();
  printf("%.2f sec\n", (double)(endTime - startTime) / CLOCKS_PER_SEC);
  startTime = endTime;

  int *common_x = (int *)calloc(n, sizeof(int));
  int *common_x_list = (int *)malloc(n * sizeof(int)), ncx = 0;
  int *common_a = (int *)calloc(n, sizeof(int));
  int *common_a_list = (int *)malloc(n * sizeof(int)), nca = 0;

  // set up a system of equations relating orbit counts
  printf("stage 3 - building systems of equations\n");
  for (int x = 0; x < n; ++x) {
    for (int i = 0; i < ncx; i++)
      common_x[common_x_list[i]] = 0;
    ncx = 0;

    // smaller graphlets
    orbit[x][0] = deg[x];
    for (int nx1 = 0; nx1 < deg[x]; nx1++) {
      int a = adj[x][nx1];
      for (int nx2 = nx1 + 1; nx2 < deg[x]; nx2++) {
        int b = adj[x][nx2];
        if ((this->*adjacent)(a, b))
          orbit[x][3]++;
        else
          orbit[x][2]++;
      }
      for (int na = 0; na < deg[a]; na++) {
        int b = adj[a][na];
        if (b != x && !(this->*adjacent)(x, b)) {
          orbit[x][1]++;
          if (common_x[b] == 0)
            common_x_list[ncx++] = b;
          common_x[b]++;
        }
      }
    }

    int64 f_71 = 0, f_70 = 0, f_67 = 0, f_66 = 0, f_58 = 0, f_57 = 0; // 14
    int64 f_69 = 0, f_68 = 0, f_64 = 0, f_61 = 0, f_60 = 0, f_55 = 0, f_48 = 0,
          f_42 = 0, f_41 = 0; // 13
    int64 f_65 = 0, f_63 = 0, f_59 = 0, f_54 = 0, f_47 = 0, f_46 = 0,
          f_40 = 0; // 12
    int64 f_62 = 0, f_53 = 0, f_51 = 0, f_50 = 0, f_49 = 0, f_38 = 0, f_37 = 0,
          f_36 = 0;                                                   // 8
    int64 f_44 = 0, f_33 = 0, f_30 = 0, f_26 = 0;                     // 1 1
    int64 f_52 = 0, f_43 = 0, f_32 = 0, f_29 = 0, f_25 = 0;           // 1 0
    int64 f_56 = 0, f_45 = 0, f_39 = 0, f_31 = 0, f_28 = 0, f_24 = 0; // 9
    int64 f_35 = 0, f_34 = 0, f_27 = 0, f_18 = 0, f_16 = 0, f_15 = 0; // 4
    int64 f_17 = 0;                                                   // 5
    int64 f_22 = 0, f_20 = 0, f_19 = 0;                               // 6
    int64 f_23 = 0, f_21 = 0;                                         // 7

    for (int nx1 = 0; nx1 < deg[x]; nx1++) {
      int a = inc[x][nx1].first, xa = inc[x][nx1].second;

      for (int i = 0; i < nca; i++)
        common_a[common_a_list[i]] = 0;
      nca = 0;
      for (int na = 0; na < deg[a]; na++) {
        int b = adj[a][na];
        for (int nb = 0; nb < deg[b]; nb++) {
          int c = adj[b][nb];
          if (c == a || (this->*adjacent)(a, c))
            continue;
          if (common_a[c] == 0)
            common_a_list[nca++] = c;
          common_a[c]++;
        }
      }

      // x = orbit-14 (tetrahedron) (clique, chordal-cycle)
      for (int nx2 = nx1 + 1; nx2 < deg[x]; nx2++) {
        int b = inc[x][nx2].first, xb = inc[x][nx2].second;
        if ((this->*adjacent)(a, b)) {
          // clique
          for (int nx3 = nx2 + 1; nx3 < deg[x]; nx3++) {
            int c = inc[x][nx3].first, xc = inc[x][nx3].second;
            if (!(this->*adjacent)(a, c) || !(this->*adjacent)(b, c))
              continue;
            orbit[x][14]++;
            f_70 += common3_get(TRIPLE(a, b, c)) - 1;
            f_71 += (tri[xa] > 2 && tri[xb] > 2)
                      ? (common3_get(TRIPLE(x, a, b)) - 1)
                      : 0;
            f_71 += (tri[xa] > 2 && tri[xc] > 2)
                      ? (common3_get(TRIPLE(x, a, c)) - 1)
                      : 0;
            f_71 += (tri[xb] > 2 && tri[xc] > 2)
                      ? (common3_get(TRIPLE(x, b, c)) - 1)
                      : 0;
            f_67 += tri[xa] - 2 + tri[xb] - 2 + tri[xc] - 2;
            f_66 += common2_get(PAIR(a, b)) - 2;
            f_66 += common2_get(PAIR(a, c)) - 2;
            f_66 += common2_get(PAIR(b, c)) - 2;
            f_58 += deg[x] - 3;
            f_57 += deg[a] - 3 + deg[b] - 3 + deg[c] - 3;
          }
          // chordal-cycle
          for (int na = 0; na < deg[a]; na++) {
            int c = inc[a][na].first, ac = inc[a][na].second;
            if (c == x || (this->*adjacent)(x, c) || !(this->*adjacent)(b, c))
              continue;
            orbit[x][12]++;
            f_65 += (tri[ac] > 1) ? common3_get(TRIPLE(a, b, c)) : 0;
            f_63 += common_x[c] - 2;
            f_59 += tri[ac] - 1 + common2_get(PAIR(b, c)) - 1;
            f_54 += common2_get(PAIR(a, b)) - 2;
            f_47 += deg[x] - 2;
            f_46 += deg[c] - 2;
            f_40 += deg[a] - 3 + deg[b] - 3;
          }
          // tailed-triangle
          for (int nx3 = 0; nx3 < deg[x]; nx3++) {
            int c = inc[x][nx3].first, xc = inc[x][nx3].second;
            if (c == a || c == b || (this->*adjacent)(a, c) ||
                (this->*adjacent)(b, c))
              continue;
            orbit[x][11]++;
            f_44 += tri[xc];
            f_33 += deg[x] - 3;
            f_30 += deg[c] - 1;
            f_26 += deg[a] - 2 + deg[b] - 2;
          }
        } else {
          // cycle
          for (int na = 0; na < deg[a]; na++) {
            int c = inc[a][na].first, ac = inc[a][na].second;
            if (c == x || (this->*adjacent)(x, c) || !(this->*adjacent)(b, c))
              continue;
            orbit[x][8]++;
            f_62 += (tri[ac] > 0) ? common3_get(TRIPLE(a, b, c)) : 0;
            f_53 += tri[xa] + tri[xb];
            f_51 += tri[ac] + common2_get(PAIR(c, b));
            f_50 += common_x[c] - 2;
            f_49 += common_a[b] - 2;
            f_38 += deg[x] - 2;
            f_37 += deg[a] - 2 + deg[b] - 2;
            f_36 += deg[c] - 2;
          }
        }
      }

      // x = orbit-7 (claw) (3-star)
      for (int nx2 = nx1 + 1; nx2 < deg[x]; nx2++) {
        int b = inc[x][nx2].first, xb = inc[x][nx2].second;
        if ((this->*adjacent)(a, b))
          continue;
        for (int nx3 = nx2 + 1; nx3 < deg[x]; nx3++) {
          int c = inc[x][nx3].first, xc = inc[x][nx3].second;
          if ((this->*adjacent)(a, c) || (this->*adjacent)(b, c))
            continue;
          orbit[x][7]++;
          f_23 += deg[x] - 3;
          f_21 += deg[a] - 1 + deg[b] - 1 + deg[c] - 1;
        }
      }

      // x = orbit-13 (diamond) (chordal-cycle)
      for (int nx2 = 0; nx2 < deg[x]; nx2++) {
        int b = inc[x][nx2].first, xb = inc[x][nx2].second;
        if ((this->*adjacent)(a, b)) {
          // chordal-cycle
          for (int nx3 = nx2 + 1; nx3 < deg[x]; nx3++) {
            int c = inc[x][nx3].first, xc = inc[x][nx3].second;
            if (!(this->*adjacent)(a, c) || (this->*adjacent)(b, c))
              continue;
            orbit[x][13]++;
            f_69 += (tri[xb] > 1 && tri[xc] > 1)
                      ? (common3_get(TRIPLE(x, b, c)) - 1)
                      : 0;
            f_68 += common3_get(TRIPLE(a, b, c)) - 1;
            f_64 += common2_get(PAIR(b, c)) - 2;
            f_61 += tri[xb] - 1 + tri[xc] - 1;
            f_60 += common2_get(PAIR(a, b)) - 1;
            f_60 += common2_get(PAIR(a, c)) - 1;
            f_55 += tri[xa] - 2;
            f_48 += deg[b] - 2 + deg[c] - 2;
            f_42 += deg[x] - 3;
            f_41 += deg[a] - 3;
          }
          // tailed-triangle
          for (int nb = 0; nb < deg[b]; nb++) {
            int c = inc[b][nb].first, bc = inc[b][nb].second;
            if (c == x || c == a || (this->*adjacent)(a, c) ||
                (this->*adjacent)(x, c))
              continue;
            orbit[x][10]++;
            f_52 += common_a[c] - 1;
            f_43 += tri[bc];
            f_32 += deg[b] - 3;
            f_29 += deg[c] - 1;
            f_25 += deg[a] - 2;
          }
        } else if (b != a) {
          // path
          for (int nb = 0; nb < deg[b]; nb++) {
            int c = inc[b][nb].first, bc = inc[b][nb].second;
            if (c == x || (this->*adjacent)(a, c) || (this->*adjacent)(x, c))
              continue;
            orbit[x][5]++;
            f_17 += deg[a] - 1;
          }
        }
      }

      // x = orbit-9 (paw) (tailed-triangle)
      for (int na1 = 0; na1 < deg[a]; na1++) {
        int b = inc[a][na1].first, ab = inc[a][na1].second;
        if (b == x || (this->*adjacent)(x, b))
          continue;
        for (int na2 = na1 + 1; na2 < deg[a]; na2++) {
          int c = inc[a][na2].first, ac = inc[a][na2].second;
          if (c == x || !(this->*adjacent)(b, c) || (this->*adjacent)(x, c))
            continue;
          orbit[x][9]++;
          f_56 +=
            (tri[ab] > 1 && tri[ac] > 1) ? common3_get(TRIPLE(a, b, c)) : 0;
          f_45 += common2_get(PAIR(b, c)) - 1;
          f_39 += tri[ab] - 1 + tri[ac] - 1;
          f_31 += deg[a] - 3;
          f_28 += deg[x] - 1;
          f_24 += deg[b] - 2 + deg[c] - 2;
        }
        // x = orbit-4 (path)
        for (int nb = 0; nb < deg[b]; nb++) {
          int c = inc[b][nb].first, bc = inc[b][nb].second;
          if (c == a || (this->*adjacent)(a, c) || (this->*adjacent)(x, c))
            continue;
          orbit[x][4]++;
          f_35 += common_a[c] - 1;
          f_34 += common_x[c];
          f_27 += tri[bc];
          f_18 += deg[b] - 2;
          f_16 += deg[x] - 1;
          f_15 += deg[c] - 1;
        }
        // x = orbit-6 (claw) (3-star)
        for (int na2 = na1 + 1; na2 < deg[a]; na2++) {
          int c = inc[a][na2].first, ac = inc[a][na2].second;
          if (c == x || (this->*adjacent)(x, c) || (this->*adjacent)(b, c))
            continue;
          orbit[x][6]++;
          f_22 += deg[a] - 3;
          f_20 += deg[x] - 1;
          f_19 += deg[b] - 1 + deg[c] - 1;
        }
      }

    } // end of for-loop for a

    // solve equations
    orbit[x][72] = C5[x];
    orbit[x][71] = (f_71 - 12 * orbit[x][72]) / 2;
    orbit[x][70] = (f_70 - 4 * orbit[x][72]);
    orbit[x][69] = (f_69 - 2 * orbit[x][71]) / 4;
    orbit[x][68] = (f_68 - 2 * orbit[x][71]);
    orbit[x][67] = (f_67 - 12 * orbit[x][72] - 4 * orbit[x][71]);
    orbit[x][66] =
      (f_66 - 12 * orbit[x][72] - 2 * orbit[x][71] - 3 * orbit[x][70]);
    orbit[x][65] = (f_65 - 3 * orbit[x][70]) / 2;
    orbit[x][64] =
      (f_64 - 2 * orbit[x][71] - 4 * orbit[x][69] - 1 * orbit[x][68]);
    orbit[x][63] = (f_63 - 3 * orbit[x][70] - 2 * orbit[x][68]);
    orbit[x][62] = (f_62 - 1 * orbit[x][68]) / 2;
    orbit[x][61] =
      (f_61 - 4 * orbit[x][71] - 8 * orbit[x][69] - 2 * orbit[x][67]) / 2;
    orbit[x][60] =
      (f_60 - 4 * orbit[x][71] - 2 * orbit[x][68] - 2 * orbit[x][67]);
    orbit[x][59] =
      (f_59 - 6 * orbit[x][70] - 2 * orbit[x][68] - 4 * orbit[x][65]);
    orbit[x][58] =
      (f_58 - 4 * orbit[x][72] - 2 * orbit[x][71] - 1 * orbit[x][67]);
    orbit[x][57] = (f_57 - 12 * orbit[x][72] - 4 * orbit[x][71] -
                    3 * orbit[x][70] - 1 * orbit[x][67] - 2 * orbit[x][66]);
    orbit[x][56] = (f_56 - 2 * orbit[x][65]) / 3;
    orbit[x][55] = (f_55 - 2 * orbit[x][71] - 2 * orbit[x][67]) / 3;
    orbit[x][54] =
      (f_54 - 3 * orbit[x][70] - 1 * orbit[x][66] - 2 * orbit[x][65]) / 2;
    orbit[x][53] =
      (f_53 - 2 * orbit[x][68] - 2 * orbit[x][64] - 2 * orbit[x][63]);
    orbit[x][52] =
      (f_52 - 2 * orbit[x][66] - 2 * orbit[x][64] - 1 * orbit[x][59]) / 2;
    orbit[x][51] =
      (f_51 - 2 * orbit[x][68] - 2 * orbit[x][63] - 4 * orbit[x][62]);
    orbit[x][50] = (f_50 - 1 * orbit[x][68] - 2 * orbit[x][63]) / 3;
    orbit[x][49] =
      (f_49 - 1 * orbit[x][68] - 1 * orbit[x][64] - 2 * orbit[x][62]) / 2;
    orbit[x][48] = (f_48 - 4 * orbit[x][71] - 8 * orbit[x][69] -
                    2 * orbit[x][68] - 2 * orbit[x][67] - 2 * orbit[x][64] -
                    2 * orbit[x][61] - 1 * orbit[x][60]);
    orbit[x][47] = (f_47 - 3 * orbit[x][70] - 2 * orbit[x][68] -
                    1 * orbit[x][66] - 1 * orbit[x][63] - 1 * orbit[x][60]);
    orbit[x][46] = (f_46 - 3 * orbit[x][70] - 2 * orbit[x][68] -
                    2 * orbit[x][65] - 1 * orbit[x][63] - 1 * orbit[x][59]);
    orbit[x][45] =
      (f_45 - 2 * orbit[x][65] - 2 * orbit[x][62] - 3 * orbit[x][56]);
    orbit[x][44] = (f_44 - 1 * orbit[x][67] - 2 * orbit[x][61]) / 4;
    orbit[x][43] =
      (f_43 - 2 * orbit[x][66] - 1 * orbit[x][60] - 1 * orbit[x][59]) / 2;
    orbit[x][42] = (f_42 - 2 * orbit[x][71] - 4 * orbit[x][69] -
                    2 * orbit[x][67] - 2 * orbit[x][61] - 3 * orbit[x][55]);
    orbit[x][41] = (f_41 - 2 * orbit[x][71] - 1 * orbit[x][68] -
                    2 * orbit[x][67] - 1 * orbit[x][60] - 3 * orbit[x][55]);
    orbit[x][40] = (f_40 - 6 * orbit[x][70] - 2 * orbit[x][68] -
                    2 * orbit[x][66] - 4 * orbit[x][65] - 1 * orbit[x][60] -
                    1 * orbit[x][59] - 4 * orbit[x][54]);
    orbit[x][39] =
      (f_39 - 4 * orbit[x][65] - 1 * orbit[x][59] - 6 * orbit[x][56]) / 2;
    orbit[x][38] = (f_38 - 1 * orbit[x][68] - 1 * orbit[x][64] -
                    2 * orbit[x][63] - 1 * orbit[x][53] - 3 * orbit[x][50]);
    orbit[x][37] = (f_37 - 2 * orbit[x][68] - 2 * orbit[x][64] -
                    2 * orbit[x][63] - 4 * orbit[x][62] - 1 * orbit[x][53] -
                    1 * orbit[x][51] - 4 * orbit[x][49]);
    orbit[x][36] = (f_36 - 1 * orbit[x][68] - 2 * orbit[x][63] -
                    2 * orbit[x][62] - 1 * orbit[x][51] - 3 * orbit[x][50]);
    orbit[x][35] =
      (f_35 - 1 * orbit[x][59] - 2 * orbit[x][52] - 2 * orbit[x][45]) / 2;
    orbit[x][34] =
      (f_34 - 1 * orbit[x][59] - 2 * orbit[x][52] - 1 * orbit[x][51]) / 2;
    orbit[x][33] = (f_33 - 1 * orbit[x][67] - 2 * orbit[x][61] -
                    3 * orbit[x][58] - 4 * orbit[x][44] - 2 * orbit[x][42]) /
                   2;
    orbit[x][32] = (f_32 - 2 * orbit[x][66] - 1 * orbit[x][60] -
                    1 * orbit[x][59] - 2 * orbit[x][57] - 2 * orbit[x][43] -
                    2 * orbit[x][41] - 1 * orbit[x][40]) /
                   2;
    orbit[x][31] = (f_31 - 2 * orbit[x][65] - 1 * orbit[x][59] -
                    3 * orbit[x][56] - 1 * orbit[x][43] - 2 * orbit[x][39]);
    orbit[x][30] = (f_30 - 1 * orbit[x][67] - 1 * orbit[x][63] -
                    2 * orbit[x][61] - 1 * orbit[x][53] - 4 * orbit[x][44]);
    orbit[x][29] = (f_29 - 2 * orbit[x][66] - 2 * orbit[x][64] -
                    1 * orbit[x][60] - 1 * orbit[x][59] - 1 * orbit[x][53] -
                    2 * orbit[x][52] - 2 * orbit[x][43]);
    orbit[x][28] = (f_28 - 2 * orbit[x][65] - 2 * orbit[x][62] -
                    1 * orbit[x][59] - 1 * orbit[x][51] - 1 * orbit[x][43]);
    orbit[x][27] =
      (f_27 - 1 * orbit[x][59] - 1 * orbit[x][51] - 2 * orbit[x][45]) / 2;
    orbit[x][26] = (f_26 - 2 * orbit[x][67] - 2 * orbit[x][63] -
                    2 * orbit[x][61] - 6 * orbit[x][58] - 1 * orbit[x][53] -
                    2 * orbit[x][47] - 2 * orbit[x][42]);
    orbit[x][25] = (f_25 - 2 * orbit[x][66] - 2 * orbit[x][64] -
                    1 * orbit[x][59] - 2 * orbit[x][57] - 2 * orbit[x][52] -
                    1 * orbit[x][48] - 1 * orbit[x][40]) /
                   2;
    orbit[x][24] = (f_24 - 4 * orbit[x][65] - 4 * orbit[x][62] -
                    1 * orbit[x][59] - 6 * orbit[x][56] - 1 * orbit[x][51] -
                    2 * orbit[x][45] - 2 * orbit[x][39]);
    orbit[x][23] =
      (f_23 - 1 * orbit[x][55] - 1 * orbit[x][42] - 2 * orbit[x][33]) / 4;
    orbit[x][22] = (f_22 - 2 * orbit[x][54] - 1 * orbit[x][40] -
                    1 * orbit[x][39] - 1 * orbit[x][32] - 2 * orbit[x][31]) /
                   3;
    orbit[x][21] = (f_21 - 3 * orbit[x][55] - 3 * orbit[x][50] -
                    2 * orbit[x][42] - 2 * orbit[x][38] - 2 * orbit[x][33]);
    orbit[x][20] = (f_20 - 2 * orbit[x][54] - 2 * orbit[x][49] -
                    1 * orbit[x][40] - 1 * orbit[x][37] - 1 * orbit[x][32]);
    orbit[x][19] = (f_19 - 4 * orbit[x][54] - 4 * orbit[x][49] -
                    1 * orbit[x][40] - 2 * orbit[x][39] - 1 * orbit[x][37] -
                    2 * orbit[x][35] - 2 * orbit[x][31]);
    orbit[x][18] = (f_18 - 1 * orbit[x][59] - 1 * orbit[x][51] -
                    2 * orbit[x][46] - 2 * orbit[x][45] - 2 * orbit[x][36] -
                    2 * orbit[x][27] - 1 * orbit[x][24]) /
                   2;
    orbit[x][17] = (f_17 - 1 * orbit[x][60] - 1 * orbit[x][53] -
                    1 * orbit[x][51] - 1 * orbit[x][48] - 1 * orbit[x][37] -
                    2 * orbit[x][34] - 2 * orbit[x][30]) /
                   2;
    orbit[x][16] = (f_16 - 1 * orbit[x][59] - 2 * orbit[x][52] -
                    1 * orbit[x][51] - 2 * orbit[x][46] - 2 * orbit[x][36] -
                    2 * orbit[x][34] - 1 * orbit[x][29]);
    orbit[x][15] = (f_15 - 1 * orbit[x][59] - 2 * orbit[x][52] -
                    1 * orbit[x][51] - 2 * orbit[x][45] - 2 * orbit[x][35] -
                    2 * orbit[x][34] - 2 * orbit[x][27]);
  } // end of for-loop for x
  endTime = clock();
  printf("%.2f sec\n", (double)(endTime - startTime) / CLOCKS_PER_SEC);

  endTime_all = endTime;
  printf("total: %.2f sec\n",
         (double)(endTime_all - startTime_all) / CLOCKS_PER_SEC);
}

int orca::init(vector<vector<int>> &origin_edges) {
  this->n = origin_edges.size();
  this->m = 0;
  adj.resize(n);
  deg.resize(n);
  for (int v = 0; v < this->n; ++v) {
    deg[v] = origin_edges[v].size();
    // we assume the original edges are sorted.
    for (auto e : origin_edges[v]) {
      if (e >= v)
        break;
      edges.push_back(PAIR(e, v));
      this->m = this->m + 1;
    }
    for (auto e : origin_edges[v])
      adj[v].push_back(e);
  }

  // set up adjacency matrix if it's smaller than 100MB
  if ((int64)n * n < 100LL * 1024 * 1024 * 8) {
    adjacent = &orca::adjacent_matrix;
    adj_matrix = (int *)calloc((n * n) / adj_chunk + 1, sizeof(int));
    for (int i = 0; i < this->m; i++) {
      int a = edges[i].a, b = edges[i].b;
      adj_matrix[(a * n + b) / adj_chunk] |= (1 << ((a * n + b) % adj_chunk));
      adj_matrix[(b * n + a) / adj_chunk] |= (1 << ((b * n + a) % adj_chunk));
    }
  } else {
    adjacent = &orca::adjacent_list;
  }
  // set up adjacency, incidence lists
  inc.resize(n);
  for (int i = 0; i < this->m; ++i) {
    int a = edges[i].a, b = edges[i].b;
    inc[a].push_back(PII(b, i));
    inc[b].push_back(PII(a, i));
  }
  // initialize orbit counts
  orbit = (int64 **)malloc(n * sizeof(int64 *));
  for (int i = 0; i < n; i++)
    orbit[i] = (int64 *)calloc(73, sizeof(int64));
  return 1;
}

vector<long long> orca::countGraphlet(int g) {
  vector<long long> result;
  if (g == 5) {
    count5();
    result.resize(21);
    std::fill(result.begin(), result.end(), 0LL);
    for (int v = 0; v < this->n; ++v) {
      result[0] += orbit[v][17];
      result[1] += orbit[v][21];
      result[2] += orbit[v][23];
      result[3] += orbit[v][25];
      result[4] += orbit[v][30];
      result[5] += orbit[v][33];
      result[6] += orbit[v][34]; // five fold
      result[7] += orbit[v][38];
      result[8] += orbit[v][42];
      result[9] += orbit[v][44];
      result[10] += orbit[v][47];
      result[11] += orbit[v][50]; // two fold
      result[12] += orbit[v][52];
      result[13] += orbit[v][55]; // two fold
      result[14] += orbit[v][58];
      result[15] += orbit[v][61];
      result[16] += orbit[v][62];
      result[17] += orbit[v][65];
      result[18] += orbit[v][69];
      result[19] += orbit[v][70]; // two fold
      result[20] += orbit[v][72]; // five fold
    }
    result[6] /= 5;
    result[11] /= 2;
    result[13] /= 2;
    result[19] /= 2;
    result[20] /= 5;
  } else if (g == 4) {
    count4();
    result.resize(6);
    std::fill(result.begin(), result.end(), 0LL);
    for (int v = 0; v < this->n; ++v) {
      result[0] += orbit[v][4]; // two fold
      result[1] += orbit[v][7];
      result[2] += orbit[v][8]; // four fold
      result[3] += orbit[v][11];
      result[4] += orbit[v][12]; // two fold
      result[5] += orbit[v][14]; // four fold
    }
    result[0] /= 2;
    result[2] /= 4;
    result[4] /= 2;
    result[5] /= 4;
  }
  return result;
}
