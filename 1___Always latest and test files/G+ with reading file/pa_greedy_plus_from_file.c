// pa_greedy_plus_from_file.c (with progress logging)
// Deterministic Greedy+ (prefix-pruning) for (n, d)-Permutation Arrays
// Seeds from a file, extends deterministically (no randomness), and logs progress.
//
// Build:
//   gcc -O3 -march=native pa_greedy_plus_from_file.c -o pa_greedy_plus
//
// Run example:
//   ./pa_greedy_plus 13 12 pa_r_13_12.txt gplus_pa_13_12.txt
//
// Progress:
//   Prints a message to stderr every PROGRESS_STEP additions beyond the seed count.
//   You can change PROGRESS_STEP by editing the macro below.

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#ifndef MAX_N
#define MAX_N 28
#endif

#ifndef MAX_PA
#define MAX_PA 200000
#endif

#ifndef PROGRESS_STEP
#define PROGRESS_STEP 1
#endif

static int n_global, d_global;
static int T_global;                 // T = n - d/2
static int *pa;                      // MAX_PA x n
static int pa_count = 0;
static int seed_count = 0;           // how many came from the seed file
static int pi[MAX_N];
static int used[MAX_N+1];

static long long nodes = 0, pruned = 0, pushed = 0;

// ---------- helpers ----------
static inline int idx(int row, int col) { return row * n_global + col; }

static int same_perm(const int *a, const int *b) {
    for (int i = 0; i < n_global; ++i) if (a[i] != b[i]) return 0;
    return 1;
}

// LCS for prefix a[0..m-1] versus full b[0..n-1]
static int lcs_prefix(const int *a, int m, const int *b) {
    static int row0[MAX_N+1], row1[MAX_N+1];
    for (int j = 0; j <= n_global; ++j) row0[j] = 0;
    for (int i = 1; i <= m; ++i) {
        row1[0] = 0;
        for (int j = 1; j <= n_global; ++j) {
            if (a[i-1] == b[j-1]) row1[j] = row0[j-1] + 1;
            else row1[j] = (row0[j] > row1[j-1]) ? row0[j] : row1[j-1];
        }
        for (int j = 0; j <= n_global; ++j) row0[j] = row1[j];
    }
    return row0[n_global];
}

static int lcs_full(const int *a, const int *b) { return lcs_prefix(a, n_global, b); }

static inline int lev_distance_full(const int *a, const int *b) {
    int l = lcs_full(a, b);
    return 2 * (n_global - l);
}

static int prefix_ok(int m) {
    if (m <= 0) return 1;
    for (int r = 0; r < pa_count; ++r) {
        const int *P = &pa[idx(r,0)];
        int lcp = lcs_prefix(pi, m, P);
        if (lcp > T_global) return 0;  // prune
    }
    return 1;
}

static int full_ok() {
    for (int r = 0; r < pa_count; ++r) {
        const int *P = &pa[idx(r,0)];
        if (lev_distance_full(pi, P) < d_global) return 0;
    }
    return 1;
}

static void progress_maybe() {
    int added = pa_count - seed_count;
    if (added > 0 && (added % PROGRESS_STEP == 0)) {
        fprintf(stderr, "Found %d permutations so far (added %d after seed)...\n", pa_count, added);
    }
}

// deterministic backtracking with prefix pruning
static void backtrack(int l) {
    if (l == n_global) {
        if (full_ok()) {
            if (pa_count < MAX_PA) {
                memcpy(&pa[idx(pa_count,0)], pi, sizeof(int)*n_global);
                ++pa_count;
                progress_maybe();
            }
        }
        ++nodes;
        return;
    }
    for (int i = l; i < n_global; ++i) {
        int tmp = pi[l]; pi[l] = pi[i]; pi[i] = tmp;
        if (prefix_ok(l+1)) {
            backtrack(l+1);
            ++pushed;
        } else {
            ++pruned;
        }
        tmp = pi[l]; pi[l] = pi[i]; pi[i] = tmp;
    }
}

static int exists_in_pa(const int *perm) {
    for (int r = 0; r < pa_count; ++r)
        if (same_perm(perm, &pa[idx(r,0)])) return 1;
    return 0;
}

static void load_seed_file(const char *path) {
    if (!path || !*path) return;
    FILE *fp = fopen(path, "r");
    if (!fp) {
        fprintf(stderr, "WARN: Could not open seed file '%s' — starting from empty PA.\n", path);
        return;
    }
    int *buf = (int*)malloc(sizeof(int)*n_global);
    if (!buf) { fprintf(stderr, "ERROR: memory alloc for seed buf failed.\n"); fclose(fp); return; }
    while (1) {
        int ok = 1;
        for (int i = 0; i < n_global; ++i) {
            if (fscanf(fp, "%d", &buf[i]) != 1) { ok = 0; break; }
        }
        if (!ok) break;
        if (!exists_in_pa(buf)) {
            if (pa_count < MAX_PA) memcpy(&pa[idx(pa_count++,0)], buf, sizeof(int)*n_global);
            else { fprintf(stderr, "ERROR: MAX_PA capacity exceeded while loading seed file.\n"); break; }
        }
    }
    free(buf);
    fclose(fp);
}

static void write_pa(const char *out_path) {
    FILE *fp = fopen(out_path, "w");
    if (!fp) { fprintf(stderr, "ERROR: Could not open output file '%s'.\n", out_path); return; }
    for (int r = 0; r < pa_count; ++r) {
        const int *P = &pa[idx(r,0)];
        for (int j = 0; j < n_global; ++j) {
            if (j) fputc(' ', fp);
            fprintf(fp, "%d", P[j]);
        }
        fputc('\n', fp);
    }
    fclose(fp);
}

int main(int argc, char **argv) {
    if (argc < 3 || argc > 5) {
        fprintf(stderr, "Usage: %s <n> <d> [seed_file] [out_file]\n", argv[0]);
        fprintf(stderr, "Example: %s 13 12 pa_r_13_12.txt gplus_pa_13_12.txt\n", argv[0]);
        return 1;
    }
    n_global = atoi(argv[1]);
    d_global = atoi(argv[2]);
    if (n_global <= 0 || n_global > MAX_N || d_global < 0 || d_global > 2*n_global) {
        fprintf(stderr, "ERROR: invalid <n> or <d>. (n in 1..%d, d in 0..%d)\n", MAX_N, 2*MAX_N);
        return 1;
    }
    T_global = n_global - (d_global/2);

    const char *seed_file = (argc >= 4) ? argv[3] : NULL;
    char default_out[128];
    snprintf(default_out, sizeof(default_out), "gplus_pa_%d_%d.txt", n_global, d_global);
    const char *out_file  = (argc >= 5) ? argv[4] : default_out;

    pa = (int*)malloc((size_t)MAX_PA * (size_t)n_global * sizeof(int));
    if (!pa) { fprintf(stderr, "ERROR: memory alloc failed for PA buffer.\n"); return 2; }

    for (int i = 0; i < n_global; ++i) pi[i] = i+1;

    load_seed_file(seed_file);
    seed_count = pa_count;
    fprintf(stderr, "Loaded %d seed permutations from '%s'.\n", seed_count, seed_file ? seed_file : "(none)");

    clock_t start = clock();
    backtrack(0);
    clock_t end = clock();

    write_pa(out_file);

    double secs = (double)(end - start) / CLOCKS_PER_SEC;
    fprintf(stderr, "n=%d d=%d | PA size=%d (seed %d, added %d) | nodes=%lld, pruned=%lld, pushed=%lld | time=%.3fs\n",
            n_global, d_global, pa_count, seed_count, pa_count - seed_count, nodes, pruned, pushed, secs);
    fprintf(stdout, "Output written to: %s\n", out_file);

    free(pa);
    return 0;
}
