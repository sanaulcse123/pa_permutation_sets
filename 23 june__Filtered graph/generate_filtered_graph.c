#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAX_N 10
#define MAX_PERMS 4000000  // up to 9! = 362880 permutations

int** perm_list;
int perm_count = 0;
int** edge_list;
int edge_count = 0;

int computeLCS(int* a, int* b, int n) {
    int dp[MAX_N+1][MAX_N+1] = {{0}};
    for (int i = 1; i <= n; ++i)
        for (int j = 1; j <= n; ++j)
            dp[i][j] = (a[i-1] == b[j-1]) ? 1 + dp[i-1][j-1] : (dp[i-1][j] > dp[i][j-1] ? dp[i-1][j] : dp[i][j-1]);
    return dp[n][n];
}

int levenshtein(int* a, int* b, int n) {
    return 2 * (n - computeLCS(a, b, n));
}

void swap(int* a, int* b) {
    int t = *a;
    *a = *b;
    *b = t;
}

void generate(int* arr, int start, int n, int d, int* identity) {
    if (start == n) {
        if (levenshtein(arr, identity, n) >= d) {
            perm_list[perm_count] = malloc(n * sizeof(int));
            memcpy(perm_list[perm_count], arr, n * sizeof(int));
            perm_count++;
        }
        return;
    }
    for (int i = start; i < n; ++i) {
        swap(&arr[start], &arr[i]);
        generate(arr, start + 1, n, d, identity);
        swap(&arr[start], &arr[i]);
    }
}

int main() {
    int n, d;
    printf("Enter n (≤ 9): ");
    scanf("%d", &n);
    printf("Enter d (Levenshtein distance threshold): ");
    scanf("%d", &d);

    if (n > 9) {
        printf("n too large. Please use n ≤ 9 for this version.\n");
        return 1;
    }

    // Allocate memory
    perm_list = malloc(MAX_PERMS * sizeof(int*));
    edge_list = malloc(MAX_PERMS * MAX_N * sizeof(int*));

    for (int i = 0; i < MAX_PERMS * MAX_N; i++) {
        edge_list[i] = malloc(2 * sizeof(int));
    }

    int identity[MAX_N], perm[MAX_N];
    for (int i = 0; i < n; i++)
        identity[i] = perm[i] = i + 1;

    // Step 1: Generate filtered permutations
    generate(perm, 0, n, d, identity);

    printf("Filtered permutations count: %d\n", perm_count);

    // Step 2: Compute graph edges
    for (int i = 0; i < perm_count; i++) {
        for (int j = i + 1; j < perm_count; j++) {
            if (levenshtein(perm_list[i], perm_list[j], n) >= d) {
                edge_list[edge_count][0] = i + 1;
                edge_list[edge_count][1] = j + 1;
                edge_count++;
            }
        }
    }

    printf("Total edges in filtered graph: %d\n", edge_count);

    // Step 3: Write filtered permutations
    char perm_file[100];
    sprintf(perm_file, "filtered_permutations_n%d_d%d.txt", n, d);
    FILE* fperm = fopen(perm_file, "w");
    for (int i = 0; i < perm_count; i++) {
        for (int j = 0; j < n; j++)
            fprintf(fperm, "%d%c", perm_list[i][j], (j == n - 1) ? '\n' : ' ');
    }
    fclose(fperm);

    // Step 4: Write DIMACS file
    char graph_file[100];
    sprintf(graph_file, "graph_n%d_d%d_filtered.dimacs", n, d);
    FILE* fgraph = fopen(graph_file, "w");
    fprintf(fgraph, "p edge %d %d\n", perm_count, edge_count);
    for (int i = 0; i < edge_count; i++)
        fprintf(fgraph, "e %d %d\n", edge_list[i][0], edge_list[i][1]);
    fclose(fgraph);

    // Cleanup
    for (int i = 0; i < perm_count; i++)
        free(perm_list[i]);
    for (int i = 0; i < edge_count; i++)
        free(edge_list[i]);
    free(perm_list);
    free(edge_list);

    printf("Output complete. Files:\n  - %s\n  - %s\n", perm_file, graph_file);
    return 0;
}
