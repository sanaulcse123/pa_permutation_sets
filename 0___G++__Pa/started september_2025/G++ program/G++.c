
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define MAX_N 18
#define MAX_PA 362880  // practical PA sizes are far smaller; static cap

// -------------------------------
// Globals
// -------------------------------
static int pa[MAX_PA][MAX_N];                 // stored codewords
static int dp[MAX_PA][MAX_N + 1][MAX_N + 1];  // incremental LCS DP per stored codeword
static int pos_map[MAX_PA][MAX_N + 1];        // pos_map[i][val] = position (1..n) of 'val' in pa[i]

static int cand[MAX_N];   // current candidate (sigma)
static int n = 13, d = 18;
static int pa_count = 0;
static int LMAX;          // = n - d/2

// -------------------------------
// Utilities
// -------------------------------
static inline void swap_int(int *a, int *b) { int t = *a; *a = *b; *b = t; }

// Full LCS (O(n^2)) vs pa[idx]; also (re)initializes dp[idx] for this 'a'
static int lcs_full(int *a, int m, int idx) {
    for (int i = 0; i <= n; ++i) dp[idx][i][0] = 0;
    for (int j = 0; j <= m; ++j) dp[idx][0][j] = 0;
    for (int j = 1; j <= m; ++j) {
        for (int i = 1; i <= n; ++i) {
            if (a[j-1] == pa[idx][i-1]) dp[idx][i][j] = dp[idx][i-1][j-1] + 1;
            else {
                int up = dp[idx][i-1][j];
                int left = dp[idx][i][j-1];
                dp[idx][i][j] = (up > left) ? up : left;
            }
        }
    }
    return dp[idx][n][m];
}

// Incremental LCS update when we extend prefix to length m (uses dp[idx][*][m-1])
static int lcs_update_prefix(int *a, int m, int idx) {
    if (m == 1) {
        for (int i = 0; i <= n; ++i) dp[idx][i][0] = 0;
    }
    dp[idx][0][m] = 0;
    for (int i = 1; i <= n; ++i) {
        if (a[m-1] == pa[idx][i-1]) dp[idx][i][m] = dp[idx][i-1][m-1] + 1;
        else {
            int up = dp[idx][i-1][m];
            int left = dp[idx][i][m-1];
            dp[idx][i][m] = (up > left) ? up : left;
        }
    }
    return dp[idx][n][m];
}

// Compute end position t_L of a length-L increasing subsequence for the sequence of positions
// pos_map[idx][cand[0..m-1]] using patience sorting tails. If LIS < L, returns -1.
static int lis_endpos_for_length(int idx, int m, int L) {
    if (L <= 0) return -1;
    int tails[MAX_N]; // tails[k] = minimal ending position for length k+1 subseq
    int len = 0;
    for (int j = 0; j < m; ++j) {
        int p = pos_map[idx][ cand[j] ];
        // lower_bound on tails[0..len)
        int lo = 0, hi = len;
        while (lo < hi) {
            int mid = (lo + hi) >> 1;
            if (tails[mid] >= p) hi = mid; else lo = mid + 1;
        }
        tails[lo] = p;
        if (lo == len) ++len;
        if (len > L) {
            // G+ would already prune when LCS exceeds LMAX; continue anyway.
        }
    }
    if (len < L) return -1;
    return tails[L-1];
}

// After accepting pa[idx], fill pos_map for fast position queries.
static void compute_pos_map(int idx) {
    for (int i = 0; i < n; ++i) {
        int v = pa[idx][i];
        pos_map[idx][v] = i + 1; // 1-based position
    }
}

// -------------------------------
// Pruning checks
// -------------------------------

// Return 1 if current prefix cand[0..m-1] passes all checks; else 0 to prune.
static int is_valid_prefix(int m) {
    for (int idx = 0; idx < pa_count; ++idx) {
        // --- G+: incremental LCS ---
        int k = lcs_update_prefix(cand, m, idx);
        if (k > LMAX) return 0; // distance < d -> prune

        // --- G++: boundary "forced-(LMAX+1)" ---
        if (k == LMAX) {
            int tL = lis_endpos_for_length(idx, m, LMAX);
            if (tL != -1) {
                int maxPosRem = -1;
                for (int j = m; j < n; ++j) {
                    int v = cand[j];
                    int p = pos_map[idx][v];
                    if (p > maxPosRem) maxPosRem = p;
                }
                if (maxPosRem > tL) return 0; // an extra match is inevitable
            }
        }
    }
    return 1;
}

// At a leaf, verify full LCS against all stored codewords (defensive final check)
static int is_valid_full(void) {
    for (int idx = 0; idx < pa_count; ++idx) {
        int k = lcs_full(cand, n, idx);
        if (k > LMAX) return 0;
    }
    return 1;
}

// -------------------------------
// DFS search
// -------------------------------
static void dfs(int l) {
    if (l == n) {
        if (is_valid_full()) {
            memcpy(pa[pa_count], cand, n * sizeof(int));
            compute_pos_map(pa_count);
            // make sure dp boundary rows/cols are zeroed for this codeword
            for (int i = 0; i <= n; ++i) dp[pa_count][i][0] = 0;
            for (int j = 0; j <= n; ++j) dp[pa_count][0][j] = 0;
            ++pa_count;
        }
        return;
    }
    for (int i = l; i < n; ++i) {
        swap_int(&cand[l], &cand[i]);
        if (is_valid_prefix(l + 1)) dfs(l + 1);
        swap_int(&cand[l], &cand[i]);
    }
}

int main(int argc, char **argv) {
    if (argc >= 2) n = atoi(argv[1]);
    if (argc >= 3) d = atoi(argv[2]);
    LMAX = n - (d / 2);

    clock_t start = clock();

    // initialize candidate to identity
    for (int i = 0; i < n; ++i) cand[i] = i + 1;

    // run
    dfs(0);

    // report
    printf("n=%d, d=%d (LMAX=%d) -> PA size: %d\n", n, d, LMAX, pa_count);

    // write output (same naming as your G+ program)
    char filename[64];
    snprintf(filename, sizeof(filename), "gplus_pa_%d_%d.txt", n, d);
    FILE *fp = fopen(filename, "w");
    if (!fp) {
        fprintf(stderr, "Error: Cannot open output file %s\n", filename);
    } else {
        for (int i = 0; i < pa_count; ++i) {
            for (int j = 0; j < n; ++j) {
                if (j) fputc(' ', fp);
                fprintf(fp, "%d", pa[i][j]);
            }
            fputc('\n', fp);
        }
        fclose(fp);
    }

    double secs = (double)(clock() - start) / CLOCKS_PER_SEC;
    printf("Total running time: %.3f seconds\n", secs);

    return 0;
}
