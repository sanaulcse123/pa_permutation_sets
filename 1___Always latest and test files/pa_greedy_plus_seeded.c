// pa_greedy_plus_seeded.c
// Seed from a permutation file, then extend the PA for (n,d) by randomized greedy.
// Distance used: Levenshtein = 2 * (n - LCS), pairwise >= d required.

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define MAX_N 64
#define MAX_PA 200000   // raise if you expect very large sets

// ---- globals ----
static int n = 13, d = 12;
static int (*PA)[MAX_N];     // PA[k][i] = value at position i (0..n-1)
static int pa_count = 0;
static int max_pa_cap = 0;

// -------- utilities --------
static void *xmalloc(size_t sz) {
    void *p = malloc(sz);
    if (!p) { fprintf(stderr, "OOM (%zu bytes)\n", sz); exit(1); }
    return p;
}
static int cmp_int(const void *a, const void *b){ int x=*(const int*)a,y=*(const int*)b; return (x>y)-(x<y); }

// Fisher–Yates shuffle
static void shuffle(int *a, int m){
    for (int i = m-1; i > 0; --i){
        int j = rand() % (i+1);
        int t=a[i]; a[i]=a[j]; a[j]=t;
    }
}

// Compute LCS length of two permutations p and q in O(n log n)
// by mapping p through inverse of q and taking LIS length.
static int lcs_len_perm(const int *p, const int *q){
    static int inv[MAX_N], seq[MAX_N], tails[MAX_N];
    for (int i=0;i<n;++i) inv[q[i]] = i;             // value -> index in q
    for (int i=0;i<n;++i) seq[i] = inv[p[i]];        // mapped sequence
    // LIS on seq[]
    int len = 0;
    for (int i=0;i<n;++i){
        int x = seq[i];
        int lo=0, hi=len;
        while (lo < hi){
            int mid = (lo+hi)>>1;
            if (tails[mid] < x) lo = mid+1; else hi = mid;
        }
        tails[lo] = x;
        if (lo == len) ++len;
    }
    return len; // LCS length
}

static int lev_distance_perm(const int *a, const int *b){
    int l = lcs_len_perm(a,b);
    return 2 * (n - l);
}

// Check if candidate is valid vs all in current PA
static int is_valid_vs_set(const int *cand){
    for (int i=0;i<pa_count;++i){
        if (lev_distance_perm(cand, PA[i]) < d) return 0;
    }
    return 1;
}

// Append to PA
static void append_perm(const int *perm){
    memcpy(PA[pa_count], perm, n*sizeof(int));
    ++pa_count;
}

// Try to parse one permutation (n ints) from a line buffer.
// Returns 1 if exactly n ints 1..n were read, else 0.
static int parse_perm_line(const char *line, int *out){
    int cnt=0; const char *p=line; char *end;
    while (*p && cnt < n){
        long v = strtol(p, &end, 10);
        if (end==p) break;
        out[cnt++] = (int)v;
        p=end;
    }
    if (cnt != n) return 0;
    // basic sanity: values in 1..n (no full uniqueness check for speed)
    for (int i=0;i<n;++i) if (out[i] < 1 || out[i] > n) return 0;
    return 1;
}

// Load seed permutations from file; verifies pairwise validity.
static void load_seed(const char *path){
    FILE *fp = fopen(path, "r");
    if (!fp){ perror("open seed"); exit(1); }
    char *line = NULL; size_t cap=0;
    int tmp[MAX_N];
    while (getline(&line, &cap, fp) != -1){
        if (!parse_perm_line(line, tmp)) continue;
        // if first, accept; else ensure it’s valid vs all already loaded
        if (!is_valid_vs_set(tmp)){
            fprintf(stderr, "[WARN] Seed file contains a permutation that violates d=%d; skipping.\n", d);
            continue;
        }
        append_perm(tmp);
        if (pa_count >= max_pa_cap-1){
            max_pa_cap *= 2;
            PA = realloc(PA, (size_t)max_pa_cap * sizeof(*PA));
            if (!PA){ fprintf(stderr,"OOM while expanding PA\n"); exit(1); }
        }
    }
    free(line);
    fclose(fp);
    fprintf(stderr, "Loaded seed: %d permutations.\n", pa_count);
}

// Write current PA to file
static void dump_pa(const char *out_path){
    FILE *fp = fopen(out_path, "w");
    if (!fp){ perror("open output"); return; }
    for (int i=0;i<pa_count;++i){
        for (int j=0;j<n;++j){
            if (j) fputc(' ', fp);
            fprintf(fp, "%d", PA[i][j]);
        }
        fputc('\n', fp);
    }
    fclose(fp);
}

// Generate next valid permutation by random sampling.
// Returns 1 if a new perm was added, 0 if none found in 'tries'.
static int random_extend(int tries){
    int cand[MAX_N];
    // start from identity then shuffle each try
    for (int i=0;i<n;++i) cand[i]=i+1;
    for (int t=0;t<tries;++t){
        shuffle(cand, n);
        if (is_valid_vs_set(cand)){
            append_perm(cand);
            return 1;
        }
    }
    return 0;
}

static void usage(const char *argv0){
    fprintf(stderr,
        "Usage:\n"
        "  %s [-n 13] [-d 12] -i seed.txt [-o out.txt] [-T tries_per_round] [-R rounds] [-s seed]\n"
        "Default: n=13, d=12, out=gplus_pa_13_12.txt, tries_per_round=200000, rounds=200, seed=time\n",
        argv0);
}

int main(int argc, char **argv){
    const char *in_path = NULL;
    const char *out_path = "gplus_pa_13_12.txt";
    int tries_per_round = 200000; // increase if CPU is fast
    int rounds = 200;
    unsigned int rseed = (unsigned int)time(NULL);

    // parse CLI
    for (int i=1;i<argc;++i){
        if (!strcmp(argv[i],"-n") && i+1<argc) n = atoi(argv[++i]);
        else if (!strcmp(argv[i],"-d") && i+1<argc) d = atoi(argv[++i]);
        else if (!strcmp(argv[i],"-i") && i+1<argc) in_path = argv[++i];
        else if (!strcmp(argv[i],"-o") && i+1<argc) out_path = argv[++i];
        else if (!strcmp(argv[i],"-T") && i+1<argc) tries_per_round = atoi(argv[++i]);
        else if (!strcmp(argv[i],"-R") && i+1<argc) rounds = atoi(argv[++i]);
        else if (!strcmp(argv[i],"-s") && i+1<argc) rseed = (unsigned)strtoul(argv[++i],NULL,10);
        else { usage(argv[0]); return 1; }
    }
    if (!in_path){ usage(argv[0]); return 1; }
    if (n > MAX_N){ fprintf(stderr,"n too big for MAX_N=%d\n", MAX_N); return 1; }

    srand(rseed);

    // allocate PA storage
    max_pa_cap = 4096;
    PA = xmalloc((size_t)max_pa_cap * sizeof(*PA));

    // 1) load seed
    load_seed(in_path);
    if (pa_count == 0){
        fprintf(stderr, "No valid seed permutations were loaded. Exiting.\n");
        return 1;
    }

    // 2) verify the seed itself (pairwise) one more time
    for (int i=0;i<pa_count;i++){
        for (int j=i+1;j<pa_count;j++){
            int lev = lev_distance_perm(PA[i], PA[j]);
            if (lev < d){
                fprintf(stderr, "[ERROR] Seed pair (%d,%d) has distance %d < %d. Fix the seed file.\n",
                        i, j, lev, d);
                return 2;
            }
        }
    }

    // 3) greedy extension
    fprintf(stderr, "Extending with randomized greedy: n=%d d=%d seed=%u\n", n, d, rseed);
    clock_t t0 = clock();
    int added_total = 0;
    for (int r=0;r<rounds;++r){
        int before = pa_count;
        int added = 0;
        // try many random candidates this round
        for (int k=0;k<tries_per_round;++k){
            if (random_extend(1)) { ++added; }
        }
        added_total += added;

        // dump progress every round
        dump_pa(out_path);

        double elapsed = (double)(clock()-t0) / CLOCKS_PER_SEC;
        fprintf(stderr, "[round %d] added=%d  total=%d  elapsed=%.1fs  wrote=%s\n",
                r+1, added, pa_count, elapsed, out_path);

        // simple stopping condition: if a few rounds add nothing, break
        static int idle = 0;
        if (pa_count == before) { if (++idle >= 5) break; }
        else idle = 0;
    }

    // final write
    dump_pa(out_path);
    fprintf(stderr, "Done. Final PA size: %d\n", pa_count);
    return 0;
}
