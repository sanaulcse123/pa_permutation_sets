#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define MAX_N 18
#define MAX_PA 362880  // 9! = 362880

int pa[MAX_PA][MAX_N];
int dp[MAX_PA][MAX_N][MAX_N];
int perm_count = 0, pa_count = 0;
int n, d;
int pi[MAX_N];
int rejected_count = 0;  // New counter for rejected partial permutations

void print_perm(int *perm, int m) {
    for (int i = 0; i < m; i++) {
        printf("%d ", perm[i]);
    }
    printf("\n");
}

int lcs(int *a, int m, int ib) {
    for (int i = 0; i <= n; i++) dp[ib][i][0];
    for (int j = 0; j <= m; j++) {
        dp[ib][0][j] = 0;
        for (int i = 1; i <= n; i++) {
            if (a[j] == pa[ib][i-1])
                dp[ib][i][j+1] = dp[ib][i-1][j] + 1;
            else
                dp[ib][i][j+1] = (dp[ib][i-1][j+1] > dp[ib][i][j]) ? dp[ib][i-1][j+1] : dp[ib][i][j];
        }
    }
    return dp[ib][n][m+1];
}

int lcs_update(int *a, int m, int ib) {
    if (m == 0)
        for (int i = 0; i <= n; i++) dp[ib][i][0];
    dp[ib][0][m+1] = 0;
    for (int i = 1; i <= n; i++) {
        if (a[m] == pa[ib][i-1])
            dp[ib][i][m+1] = dp[ib][i-1][m] + 1;
        else
            dp[ib][i][m+1] = (dp[ib][i-1][m+1] > dp[ib][i][m]) ? dp[ib][i-1][m+1] : dp[ib][i][m];
    }
    return dp[ib][n][m+1];
}

int lev_distance(int *a, int m, int ib) {
    return 2 * (n - lcs_update(a, m, ib));
}

int is_valid(int m) {
    for (int i = 0; i < pa_count; i++) {
        if (lev_distance(pi, m, i) < d)
            return 0;
    }
    return 1;
}

void swap(int *a, int *b) {
    int t = *a; *a = *b; *b = t;
}

void generate_permutations(int l, int r) {
    if (l == r && is_valid(l)) {
        memcpy(pa[pa_count], pi, n * sizeof(int));
        lcs(pi, n, pa_count);
        pa_count++;
        return;
    }
    for (int i = l; i <= r; i++) {
        swap(&pi[l], &pi[i]);
        if (is_valid(l))
            generate_permutations(l + 1, r);
        else
            rejected_count++;  // Count rejected partial permutations here
        swap(&pi[l], &pi[i]);
    }
}

int main() {
    clock_t start, end;
    start = clock();

    // Example: n=13, d=14
    for (n = 13; n <= 13; n++) {
        for (d = 14; d <= 14; d++) {
            pa_count = 0;
            rejected_count = 0;
            for (int i = 0; i < n; i++) pi[i] = i + 1;
            generate_permutations(0, n - 1);
            printf("n=%d, d=%d -> %d permutations\n", n, d, pa_count);
            printf("Partial permutations rejected: %d\n", rejected_count);  // New output
            char filename[64];
            sprintf(filename, "gplus_pa_%d_%d.txt", n, d);
            FILE *fp = fopen(filename, "w");
            if (fp) {
                for (int i = 0; i < pa_count; i++) {
                    for (int j = 0; j < n; j++)
                        fprintf(fp, "%d ", pa[i][j]);
                    fprintf(fp, "\n");
                }
                fclose(fp);
            } else {
                printf("Error: Cannot open output file!\n");
            }
        }
    }

    end = clock();
    double total_time = (double)(end - start) / CLOCKS_PER_SEC;
    printf("Total running time: %.3f seconds\n", total_time);
    return 0;
}
