#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAX_N 7
#define MAX_PA 5040  // 7! = 5040

int perms[MAX_PA][MAX_N];
int pa[MAX_PA][MAX_N];
int perm_count = 0, pa_count = 0;

// Swap helper
void swap(int *a, int *b) {
    int t = *a; *a = *b; *b = t;
}

// Generate all permutations of [1..n]
void generate_permutations(int *arr, int l, int r, int n) {
    if (l == r) {
        memcpy(perms[perm_count], arr, n * sizeof(int));
        perm_count++;
        return;
    }
    for (int i = l; i <= r; i++) {
        swap(&arr[l], &arr[i]);  //  Corrected line
        generate_permutations(arr, l + 1, r, n);
        swap(&arr[l], &arr[i]);
    }
}

// LCS length between two arrays
int lcs(int *a, int *b, int n) {
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

// Levenshtein distance using theorem
int lev_distance(int *a, int *b, int n) {
    return 2 * (n - lcs(a, b, n));
}

// Check if permutation is valid to add to PA
int is_valid(int *perm, int n, int d) {
    for (int i = 0; i < pa_count; i++) {
        if (lev_distance(perm, pa[i], n) < d)
            return 0;
    }
    return 1;
}

int main() {
    for (int n = 3; n <= 7; n++) {
        for (int d = 2; d <= 10; d++) {
            perm_count = 0; pa_count = 0;
            int base[MAX_N];
            for (int i = 0; i < n; i++) base[i] = i + 1;

            generate_permutations(base, 0, n - 1, n);

            for (int i = 0; i < perm_count; i++) {
                if (is_valid(perms[i], n, d)) {
                    memcpy(pa[pa_count], perms[i], n * sizeof(int));
                    pa_count++;
                }
            }
            printf("n=%d, d=%d -> %d permutations\n", n, d, pa_count);
        }
    }
    return 0;
}
