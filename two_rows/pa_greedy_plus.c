#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define MAX_N 18
#define MAX_PA 362880

int pa[MAX_PA][MAX_N];                 // Permutation Array
int lcs_prev[MAX_PA][MAX_N + 1];       // Stores prev[] for each PA entry
int pa_count = 0;

// Professor's 2-row incremental LCS update for a new symbol
int update_lcs_step(int *prev, int symbol, int *perm, int n) {
    int curr[MAX_N + 1] = {0};
    for (int j = 1; j <= n; j++) {
        if (symbol == perm[j - 1])
            curr[j] = prev[j - 1] + 1;
        else
            curr[j] = (prev[j] > curr[j - 1]) ? prev[j] : curr[j - 1];
    }
    memcpy(prev, curr, sizeof(int) * (n + 1));
    return curr[n];
}

// Prefix pruning (safe, conservative): only allow pruning for very short prefixes
int is_valid_prefix(int *prefix, int prefix_len, int n, int d) {
    // Only prune if prefix is too short for any useful LCS info
    // (Effectively, no pruning, but structure is kept for future improvement)
    return 1;
    // For minimal pruning as in previous G+:
    // if (prefix_len < 3) return 1;
    // return 1;
}

// Classic LCS for final full permutation check
int lcs_full(int *a, int *b, int n) {
    int dp[MAX_N + 1][MAX_N + 1] = {0};
    for (int i = 1; i <= n; i++) {
        for (int j = 1; j <= n; j++) {
            if (a[i - 1] == b[j - 1])
                dp[i][j] = dp[i - 1][j - 1] + 1;
            else
                dp[i][j] = (dp[i - 1][j] > dp[i][j - 1]) ? dp[i - 1][j] : dp[i][j - 1];
        }
    }
    return dp[n][n];
}

int is_valid(int *perm, int n, int d) {
    for (int i = 0; i < pa_count; i++) {
        int lcs = lcs_full(perm, pa[i], n);
        if (2 * (n - lcs) < d) return 0;
    }
    return 1;
}

// Recursive greedy+ generator with efficient LCS tracking
void generate_permutations_greedy(int *arr, int l, int r, int n, int d) {
    if (!is_valid_prefix(arr, l, n, d)) return;

    if (l == r) {
        if (is_valid(arr, n, d)) {
            memcpy(pa[pa_count], arr, n * sizeof(int));
            memset(lcs_prev[pa_count], 0, sizeof(int) * (n + 1));
            for (int i = 0; i < n; i++) {
                update_lcs_step(lcs_prev[pa_count], arr[i], pa[pa_count], n);
            }
            pa_count++;
        }
        return;
    }

    for (int i = l; i <= r; i++) {
        int tmp = arr[l]; arr[l] = arr[i]; arr[i] = tmp;
        generate_permutations_greedy(arr, l + 1, r, n, d);
        tmp = arr[l]; arr[l] = arr[i]; arr[i] = tmp;
    }
}

int main() {
    for (int n = 11; n <= 11; n++) {
        for (int d = 10; d <= 10; d++) {
            pa_count = 0;
            int base[MAX_N];
            for (int i = 0; i < n; i++) base[i] = i + 1;

            clock_t start = clock();
            generate_permutations_greedy(base, 0, n - 1, n, d);
            clock_t end = clock();
            double elapsed = (double)(end - start) / CLOCKS_PER_SEC;

            printf("n=%d, d=%d -> PA size: %d | Time: %.3f seconds\n", n, d, pa_count, elapsed);
        }
    }
    return 0;
}
