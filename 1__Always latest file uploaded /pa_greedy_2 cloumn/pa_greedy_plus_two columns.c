#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAX_N 18
#define MAX_PA 362880  // 9! = 362880

int pa[MAX_PA][MAX_N];
int perm_count = 0, pa_count = 0;

// Two-row space-optimized LCS (professor's algorithm)
int lcs_tworow(int *a, int *b, int n) {
    int prev[MAX_N+1] = {0}, curr[MAX_N+1] = {0};
    for (int i = 1; i <= n; i++) {
        curr[0] = 0;
        for (int j = 1; j <= n; j++) {
            if (a[i-1] == b[j-1])
                curr[j] = prev[j-1] + 1;
            else
                curr[j] = (prev[j] > curr[j-1]) ? prev[j] : curr[j-1];
        }
        memcpy(prev, curr, sizeof(curr));
    }
    return curr[n];
}

// Levenshtein distance using theorem
int lev_distance(int *a, int *b, int n) {
    return 2 * (n - lcs_tworow(a, b, n));
}

// Check if permutation is valid to add to PA
int is_valid(int *perm, int n, int d) {
    for (int i = 0; i < pa_count; i++) {
        if (lev_distance(perm, pa[i], n) < d)
            return 0;
    }
    return 1;
}

// Swap helper
void swap(int *a, int *b) {
    int t = *a; *a = *b; *b = t;
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
    int n, d;
    // Example: n=11, d=10 (as in your table for L(11,10))
    for (n = 11; n <= 11; n++) {
        for (d = 12; d <= 12; d++) {
            pa_count = 0;
            int base[MAX_N];
            for (int i = 0; i < n; i++) base[i] = i + 1;

            generate_permutations(base, 0, n - 1, n, d);

            printf("n=%d, d=%d -> %d permutations\n", n, d, pa_count);
        }
    }
    return 0;
}
