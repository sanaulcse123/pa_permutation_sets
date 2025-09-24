// greedy_plus_plus.c — Greedy++ with lexicographic permutation generation
// Usage: ./greedy_plus_plus n d
//
// Change: enumerate σ′ lexicographically (v = 1..n, unused) at each depth,
// instead of swap-based DFS. This yields true dictionary order for all accepted
// permutations (and prints them as they are found).

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <time.h>

#define MAX_N    18
#define MAX_PA   20000   // practical cap

// Stored permutations (accepted π's)
static int pa[MAX_PA][MAX_N];
static unsigned char pos_idx[MAX_PA][MAX_N + 1]; // pos_idx[ib][v] in 1..n
static int pa_count = 0;

// LCS DP stacks per π and per depth (m = 0..n), each stack level stores a column over i=0..n
// dp_len[ib][m][i] = LCS(π[1..i], σ′[1..m])
// dp_pos[ib][m][i] = earliest last-π-position among all subsequences achieving dp_len[ib][m][i]
static unsigned char dp_len[MAX_PA][MAX_N + 1][MAX_N + 1];
static unsigned char dp_pos[MAX_PA][MAX_N + 1][MAX_N + 1];

static int n, d;
static int pi_cur[MAX_N];              // current candidate permutation under construction
static unsigned char used[MAX_N + 1];  // marks symbols in current prefix σ′
static int depth_m = 0;                // current prefix length

// --- print helper (prints the current candidate permutation) ---
static inline void print_current_perm(void) {
    for (int i = 0; i < n; ++i) {
        if (i) putchar(' ');
        printf("%d", pi_cur[i]);
    }
    putchar('\n');
}

// --- Extend DP one symbol for all stored π (σ′ grows by 'sym') ---
static inline int extend_all_dp_and_check_prune(int sym)
{
    // boundary for immediate prune: L >= n - d + 1
    const int boundary = n - d + 1;

    // Build a temporary list of unused symbols T (for β* check)
    int T_buf[MAX_N]; int T_sz = 0;
    for (int v = 1; v <= n; v++) if (!used[v]) T_buf[T_sz++] = v;

    // For each stored π, extend dp column from depth_m to depth_m+1
    for (int ib = 0; ib < pa_count; ib++) {
        unsigned char left_len = 0, left_pos = 255;

        dp_len[ib][depth_m + 1][0] = 0;
        dp_pos[ib][depth_m + 1][0] = 255;

        for (int i = 1; i <= n; i++) {
            unsigned char up_len = dp_len[ib][depth_m][i];
            unsigned char up_pos = dp_pos[ib][depth_m][i];

            unsigned char match_len = 0, match_pos = 255;
            if (pa[ib][i - 1] == sym) {
                match_len = (unsigned char)(dp_len[ib][depth_m][i - 1] + 1);
                match_pos = (unsigned char)i; // last symbol ends at i in π
            }

            // best among up, left, match; tie -> earliest last position
            unsigned char best_len = up_len, best_pos = up_pos;

            if (left_len > best_len) { best_len = left_len; best_pos = left_pos; }
            else if (left_len == best_len && left_pos < best_pos) { best_pos = left_pos; }

            if (match_len > best_len) { best_len = match_len; best_pos = match_pos; }
            else if (match_len == best_len && match_pos < best_pos) { best_pos = match_pos; }

            dp_len[ib][depth_m + 1][i] = best_len;
            dp_pos[ib][depth_m + 1][i] = best_pos;

            left_len = best_len; left_pos = best_pos;
        }

        // Check prune on the final cell i = n for this π
        int L = dp_len[ib][depth_m + 1][n];
        if (L >= boundary) return 1; // distance violated -> prune

        if (L == n - d) {
            int p_beta_star = (dp_pos[ib][depth_m + 1][n] == 255) ? INT_MAX
                                                                  : dp_pos[ib][depth_m + 1][n];
            // β* rule: if ∃ t ∈ T with posπ(t) > p_beta_star -> prune
            for (int k = 0; k < T_sz; k++) {
                int t = T_buf[k];
                if ((int)pos_idx[ib][t] > p_beta_star) return 1;
            }
        }
    }
    return 0; // not pruned
}

// --- Initialize DP stacks for a newly accepted π along current path ---
static inline void init_dp_for_new_pi_along_path(int ib_new)
{
    // Base column at m=0: all zeros
    for (int i = 0; i <= n; i++) {
        dp_len[ib_new][0][i] = 0;
        dp_pos[ib_new][0][i] = 255;
    }

    // Extend for each k = 1..depth_m using the actual path symbols pi_cur[0..k-1]
    for (int k = 1; k <= depth_m; k++) {
        int sym = pi_cur[k - 1];
        unsigned char left_len = 0, left_pos = 255;

        dp_len[ib_new][k][0] = 0;
        dp_pos[ib_new][k][0] = 255;

        for (int i = 1; i <= n; i++) {
            unsigned char up_len = dp_len[ib_new][k - 1][i];
            unsigned char up_pos = dp_pos[ib_new][k - 1][i];

            unsigned char match_len = 0, match_pos = 255;
            if (pa[ib_new][i - 1] == sym) {
                match_len = (unsigned char)(dp_len[ib_new][k - 1][i - 1] + 1);
                match_pos = (unsigned char)i;
            }

            unsigned char best_len = up_len, best_pos = up_pos;

            if (left_len > best_len) { best_len = left_len; best_pos = left_pos; }
            else if (left_len == best_len && left_pos < best_pos) { best_pos = left_pos; }

            if (match_len > best_len) { best_len = match_len; best_pos = match_pos; }
            else if (match_len == best_len && match_pos < best_pos) { best_pos = match_pos; }

            dp_len[ib_new][k][i] = best_len;
            dp_pos[ib_new][k][i] = best_pos;

            left_len = best_len; left_pos = best_pos;
        }
    }
}

// --- LEXICOGRAPHIC generator: choose next value v=1..n (unused), increasing ---
static void dfs_lex(int l)
{
    if (l == n) {
        // Reached a full permutation: print, accept, and prep DP for future checks
        print_current_perm();

        memcpy(pa[pa_count], pi_cur, n * sizeof(int));
        for (int i = 0; i < n; i++)
            pos_idx[pa_count][ pa[pa_count][i] ] = (unsigned char)(i + 1);

        init_dp_for_new_pi_along_path(pa_count);
        pa_count++;
        return;
    }

    for (int v = 1; v <= n; ++v) {
        if (used[v]) continue;          // maintain dictionary order via v = 1..n
        pi_cur[l] = v;
        used[v] = 1;

        int pruned = extend_all_dp_and_check_prune(v);
       if (!pruned) {
            depth_m++;
            dfs_lex(l + 1);
            depth_m--;
        }
        used[v] = 0;
    }
}

int main(int argc, char *argv[])
{
    clock_t t0 = clock();

    if (argc != 3) {
        fprintf(stderr, "Usage: %s n d\n", argv[0]);
        return 1;
    }
    n = atoi(argv[1]);
    d = atoi(argv[2]);

    if (n < 1 || n > MAX_N) {
        fprintf(stderr, "Error: n must be in [1, %d]\n", MAX_N);
        return 1;
    }
    if (d < 0 || d > n) {
        fprintf(stderr, "Error: d must be in [0, n]\n");
        return 1;
    }

    // init
    memset(used, 0, sizeof(used));
    depth_m = 0;
    pa_count = 0;

    // Begin lexicographic DFS
    dfs_lex(0);

    printf("n=%d, d=%d -> %d permutations\n", n, d, pa_count);

    // Write output
    char fname[64];
    snprintf(fname, sizeof(fname), "gpp_pa_%d_%d.txt", n, d);
    FILE *fp = fopen(fname, "w");
    if (fp) {
        for (int i = 0; i < pa_count; i++) {
            for (int j = 0; j < n; j++) {
                if (j) fputc(' ', fp);
                fprintf(fp, "%d", pa[i][j]);
            }
            fputc('\n', fp);
        }
        fclose(fp);
        printf("Saved: %s\n", fname);
    } else {
        fprintf(stderr, "Error: cannot open %s for writing\n", fname);
    }

    double secs = (double)(clock() - t0) / CLOCKS_PER_SEC;
    printf("Total running time: %.3f seconds\n", secs);
    return 0;
}
