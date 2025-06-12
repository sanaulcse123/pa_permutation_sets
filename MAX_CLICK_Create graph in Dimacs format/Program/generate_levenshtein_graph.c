#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAXN 9
#define MAX_PERMS 362880  // 9!

int factorial(int n) {
    int res = 1;
    for (int i = 2; i <= n; i++) res *= i;
    return res;
}

// Generate the next lexicographic permutation
int next_permutation(int *a, int n) {
    int i = n - 2;
    while (i >= 0 && a[i] >= a[i + 1]) i--;
    if (i < 0) return 0;

    int j = n - 1;
    while (a[j] <= a[i]) j--;
    int temp = a[i]; a[i] = a[j]; a[j] = temp;

    for (int k = i + 1, l = n - 1; k < l; k++, l--) {
        temp = a[k]; a[k] = a[l]; a[l] = temp;
    }
    return 1;
}

// Levenshtein distance between two permutations of equal length
int levenshtein(int *a, int *b, int n) {
    int dp[n+1][n+1];
    for (int i = 0; i <= n; i++) dp[i][0] = i;
    for (int j = 0; j <= n; j++) dp[0][j] = j;
    for (int i = 1; i <= n; i++)
        for (int j = 1; j <= n; j++)
            dp[i][j] = a[i-1] == b[j-1]
                     ? dp[i-1][j-1]
                     : 1 + (dp[i-1][j-1] < dp[i][j-1] ? (dp[i-1][j-1] < dp[i-1][j] ? dp[i-1][j-1] : dp[i-1][j]) : (dp[i][j-1] < dp[i-1][j] ? dp[i][j-1] : dp[i-1][j]));

    return dp[n][n];
}

void write_graph_dimacs(int n, int d, const char *filename) {
    int perm_count = factorial(n);
    int perms[perm_count][MAXN];

    // Generate all permutations
    int base[MAXN];
    for (int i = 0; i < n; i++) base[i] = i + 1;
    memcpy(perms[0], base, sizeof(int) * n);

    int count = 1;
    while (next_permutation(base, n)) {
        memcpy(perms[count], base, sizeof(int) * n);
        count++;
    }

    // Count edges
    int edge_count = 0;
    FILE *fp = fopen(filename, "w");
    if (!fp) {
        perror("fopen");
        return;
    }

    for (int i = 0; i < perm_count; i++) {
        for (int j = i + 1; j < perm_count; j++) {
            int dist = levenshtein(perms[i], perms[j], n);
            if (dist >= d) edge_count++;
        }
    }

    fprintf(fp, "p edge %d %d\n", perm_count, edge_count);

    for (int i = 0; i < perm_count; i++) {
        for (int j = i + 1; j < perm_count; j++) {
            int dist = levenshtein(perms[i], perms[j], n);
            if (dist >= d)
                fprintf(fp, "e %d %d\n", i + 1, j + 1); // 1-based index
        }
    }

    fclose(fp);
    printf("Graph written to %s (n=%d, d=%d, vertices=%d, edges=%d)\n",
           filename, n, d, perm_count, edge_count);
}

int main() {
    for (int n = 3; n <= 9; n++) {
        for (int d = 2; d <= n; d++) {
            char filename[256];
            snprintf(filename, sizeof(filename), "graph_n%d_d%d.dimacs", n, d);
            write_graph_dimacs(n, d, filename);
        }
    }
    return 0;
}
