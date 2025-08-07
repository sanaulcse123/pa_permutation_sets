#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>  // For measuring execution time

#define MAX_N 18
#define MAX_PA 362880  // 9! = 362880

int perms[MAX_PA][MAX_N];
int pa[MAX_PA][MAX_N];
int perm_count = 0, pa_count = 0;

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

// Swap helper
void swap(int *a, int *b) {
    int t = *a; *a = *b; *b = t;
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

// Generate all permutations of [1..n]
void generate_permutations(int *arr, int l, int r, int n, int d) {
    if (l == r) {
        if (is_valid(arr, n, d)) {
            memcpy(pa[pa_count], arr, n * sizeof(int));
            pa_count++;
        }
        return;
    }
    for (int i = l; i <= r; i++) {
        swap(&arr[l], &arr[i]);
        generate_permutations(arr, l + 1, r, n, d);
        swap(&arr[l], &arr[i]);
    }
}

int main() {
    clock_t start = clock();  // Start timing

    for (int n = 13; n <= 13; n++) {
        for (int d = 14; d <= 14; d++) {
            perm_count = 0; pa_count = 0;
            int base[MAX_N];
            for (int i = 0; i < n; i++) base[i] = i + 1;

            generate_permutations(base, 0, n - 1, n, d);
            printf("n=%d, d=%d -> %d permutations\n", n, d, pa_count);
        }
    }

    clock_t end = clock();  // End timing
    double total_time = (double)(end - start) / CLOCKS_PER_SEC;
    printf("Total execution time: %.2f seconds\n", total_time);

    return 0;
}
