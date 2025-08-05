#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAX_N 12
#define MAX_PA 5000000

int pa[MAX_PA][MAX_N];
int prev_row[MAX_PA][MAX_N + 1];  // Store previous row for each PA
int curr_row[MAX_PA][MAX_N + 1];  // For temporary computation
int pa_count = 0;

// Return max of two ints
int max(int a, int b) { return (a > b) ? a : b; }

// ----- Constant-time prefix pruning using 2-row LCS -----
int is_valid_prefix(int *prefix, int l, int n, int d) {
    for (int i = 0; i < pa_count; i++) {
    curr_row[i][0] = 0; // base case for LCS

    for (int j = 1; j <= n; j++) {
        if (prefix[l-1] == pa[i][j-1])
            curr_row[i][j] = prev_row[i][j-1] + 1;
        else
            curr_row[i][j] = (prev_row[i][j] > curr_row[i][j-1])
                            ? prev_row[i][j]
                            : curr_row[i][j-1];
    }
    int lcs = curr_row[i][n];
    int lev = 2 * (n - lcs);
    if (lev < d) return 0;
}

    return 1;
}

// ----- Swap helper -----
void swap(int *a, int *b) { int t = *a; *a = *b; *b = t; }

// ----- Recursive Greedy permutation generator -----
void generate_permutations(int *arr, int l, int r, int n, int d) {
    if (l > 0 && !is_valid_prefix(arr, l, n, d)) return;

    if (l == r) {
        if (pa_count >= MAX_PA) {
            printf("Max PA limit reached (%d), aborting.\n", MAX_PA);
            exit(1);
        }
        memcpy(pa[pa_count], arr, n * sizeof(int));
        // Store current row for future prefix extensions
        memcpy(prev_row[pa_count], curr_row[pa_count], (n + 1) * sizeof(int));
        pa_count++;
        return;
    }
    for (int i = l; i <= r; i++) {
        swap(&arr[l], &arr[i]);
        generate_permutations(arr, l + 1, r, n, d);
        swap(&arr[l], &arr[i]);
    }
}

// ----- Main -----
int main() {
    int n = 10, d = 12; // Change as needed for your experiment
    int base[MAX_N];
    for (int i = 0; i < n; i++) base[i] = i;

    generate_permutations(base, 0, n - 1, n, d);

    printf("n=%d, d=%d -> %d permutations (PA lower bound)\n", n, d, pa_count);
    return 0;
}
