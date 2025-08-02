#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define MAX_N 18
#define MAX_PA 362880

int pa[MAX_PA][MAX_N];
int pa_count = 0;

// Optimized LCS for partial arrays
int lcs_partial(int *a, int a_len, int *b, int b_len) {
    int dp[MAX_N + 1][MAX_N + 1] = {0};
    for (int i = 1; i <= a_len; i++) {
        for (int j = 1; j <= b_len; j++) {
            if (a[i - 1] == b[j - 1])
                dp[i][j] = dp[i - 1][j - 1] + 1;
            else
                dp[i][j] = (dp[i - 1][j] > dp[i][j - 1]) ? dp[i - 1][j] : dp[i][j - 1];
        }
    }
    return dp[a_len][b_len];
}

int lev_distance(int *a, int *b, int n) {
    return 2 * (n - lcs_partial(a, n, b, n));
}

int is_valid(int *perm, int n, int d) {
    for (int i = 0; i < pa_count; i++) {
        if (lev_distance(perm, pa[i], n) < d)
            return 0;
    }
    return 1;
}

// Prune only if prefix length >= 3 and PA is sufficiently populated
int is_valid_prefix(int *prefix, int prefix_len, int n, int d) {
    if (prefix_len < 3 || pa_count < 10) return 1;

    for (int i = 0; i < pa_count; i++) {
        int lcs_len = lcs_partial(prefix, prefix_len, pa[i], n);
        if (2 * (n - lcs_len) < d)
            return 0;
    }
    return 1;
}

// Optimized G+ permutation generation
void generate_permutations_gplus(int *arr, int l, int r, int n, int d) {
    if (!is_valid_prefix(arr, l, n, d)) return;

    if (l == r) {
        if (is_valid(arr, n, d)) {
            memcpy(pa[pa_count], arr, n * sizeof(int));
            pa_count++;
        }
        return;
    }

    for (int i = l; i <= r; i++) {
        int tmp = arr[l]; arr[l] = arr[i]; arr[i] = tmp;
        generate_permutations_gplus(arr, l + 1, r, n, d);
        tmp = arr[l]; arr[l] = arr[i]; arr[i] = tmp;
    }
}

int main() {
    for (int n = 11; n <= 11; n++) {
        for (int d = 10; d <= 10; d++) {
            pa_count = 0;
            int base[MAX_N];
            for (int i = 0; i < n; i++) base[i] = i + 1;

            clock_t start = clock();  // Start timer

            generate_permutations_gplus(base, 0, n - 1, n, d);

            clock_t end = clock();    // End timer
            double elapsed = (double)(end - start) / CLOCKS_PER_SEC;

            printf("n=%d, d=%d -> PA size: %d | Time: %.3f seconds\n", n, d, pa_count, elapsed);
        }
    }
    return 0;
}
