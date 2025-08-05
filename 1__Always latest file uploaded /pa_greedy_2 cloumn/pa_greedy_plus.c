#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAX_N 18
#define MAX_PA 362880  // 9! = 362880

int pa[MAX_PA][MAX_N];
int dp[MAX_PA][MAX_N][MAX_N];
int perm_count = 0, pa_count = 0;
int n, d;
int pi[MAX_N];

void print_perm( int *perm,int m) {
    for (int i = 0; i < m; i++) {
        printf("%d ", perm[i]);
    }
    printf("\n");
}

int lcs(int *a,int m, int ib) {    
    //a is new perm and 
    // a[m] is a new symbol in a
    // ib is index in pa[]
    for (int i = 0; i <= n; i++) dp[ib][i][0];
     for(int j=0; j<=m; j++){
       dp[ib][0][j] = 0;
       for (int i = 1; i <= n; i++) {
        if (a[j] == pa[ib][i-1])
          dp[ib][i][j+1] = dp[ib][i-1][j] + 1;
        else
          dp[ib][i][j+1] = (dp[ib][i-1][j+1] > dp[ib][i][j]) ? dp[ib][i-1][j+1] : dp[ib][i][j];  
       }
    }
 /*
    print_perm(a,m+1);
    print_perm(pa[ib],n);
    printf("lcs=%d m=%d \n", dp[ib][n][m+1],m);
    exit(0);
    */
    return dp[ib][n][m+1];
}

int lcs_update(int *a,int m, int ib) {    
    //a is new perm and 
    // a[m] is a new symbol in a
    // ib is index in pa[]
    if( m==0 )
       for (int i = 0; i <= n; i++) dp[ib][i][0];
    dp[ib][0][m+1] = 0;
    for (int i = 1; i <= n; i++) {
        if (a[m] == pa[ib][i-1])
          dp[ib][i][m+1] = dp[ib][i-1][m] + 1;
        else
          dp[ib][i][m+1] = (dp[ib][i-1][m+1] > dp[ib][i][m]) ? dp[ib][i-1][m+1] : dp[ib][i][m];  
    }

   
   // print_perm(a,m+1);
    //print_perm(pa[ib],n);
    //printf("lcs=%d m=%d \n", dp[ib][n][m+1],m);
   // exit(0);
    
    return dp[ib][n][m+1];
}
// Levenshtein distance using theorem
int lev_distance(int *a,int m, int ib) {
    return 2 * (n - lcs_update(a, m, ib));
}

// Check if permutation is valid to add to PA
int is_valid( int m) {
    for (int i = 0; i < pa_count; i++) {
        if (lev_distance(pi ,m, i) < d)
            return 0;
    }
    return 1;
}

// Swap helper
void swap(int *a, int *b) {
    int t = *a; *a = *b; *b = t;
}

// Generate all permutations of [1..n]
void generate_permutations( int l, int r) {
    if (l == r && is_valid(l)) {
        memcpy(pa[pa_count], pi, n * sizeof(int));
        lcs(pi,n,pa_count);
        pa_count++;
        printf("New perm ");
        print_perm(pi,n);
        return; 
    }
    for (int i = l; i <= r; i++) {
        swap(&pi[l], &pi[i]);
        if (is_valid(l))
            generate_permutations(l + 1, r);
        swap(&pi[l], &pi[i]);
    }
}

int main() {
    
    // Example: n=11, d=10 (as in your table for L(11,10))
    for (n = 13; n <= 13; n++) {
        for (d = 16; d <= 16; d++) {
            pa_count = 0;
            for (int i = 0; i < n; i++) pi[i] = i + 1;
            generate_permutations(0, n - 1);
            printf("n=%d, d=%d -> %d permutations\n", n, d, pa_count);
        }
    }
    return 0;
}

