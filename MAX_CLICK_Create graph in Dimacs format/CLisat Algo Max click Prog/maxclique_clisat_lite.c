/********************************************************************
 * CliSAT-lite  –  Exact Maximum-Clique solver (bit-parallel B&B)
 *   • Reads a DIMACS “p edge … / e u v” file (undirected, 1-based)
 *   • Uses greedy colouring for an upper bound
 *   • Bit-set branch-and-bound (BBMC style)
 *   • Good up to ~50 k vertices if graph is sparse; fine for n≤9!
 *******************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

#define WORD       uint64_t
#define POPCNT(x)  __builtin_popcountll((WORD)(x))

/* ----------------------------- globals -------------------------- */
static int   N   = 0;            /* #vertices                       */
static WORD *adj = NULL;         /* pointer to adjacency bitsets    */
static int   W   = 0;            /* #64-bit words per row           */

/* best clique found */
static int  best_size = 0;
static int *best      = NULL;
static int *stack_v   = NULL;

/* ------------------ helper: access adj matrix ------------------- */
static inline WORD *ROW(int v) {          /* row pointer v (0-based) */
    return adj + (size_t)v * W;
}

/* ------------------ greedy colouring upper bound ---------------- */
static int colour[1 << 20];               /* will resize dynamically */

static int greedy_colour(const WORD *cand, int c_sz)
{
    WORD left[W]; memcpy(left, cand, W * sizeof(WORD));
    int ub = 0, remaining = c_sz;

    while (remaining) {
        ub++;
        WORD forbid[W]; memset(forbid, 0, W * sizeof(WORD));

        for (int v = 0; v < N; ++v) {
            if (!(left[v / 64] & ((WORD)1 << (v % 64)))) continue;
            colour[v] = ub;               /* colour class index      */

            /* forbid neighbours */
            WORD *row = ROW(v);
            for (int w = 0; w < W; ++w) forbid[w] |= row[w];

            /* remove v from left */
            left[v / 64] &= ~((WORD)1 << (v % 64));
            --remaining;
        }
        /* keep only vertices NOT forbidden */
        for (int w = 0; w < W; ++w) left[w] &= ~forbid[w];
    }
    return ub;
}

/* ------------------ recursive branch-and-bound ------------------ */
static void expand(int depth, WORD *cand, int c_sz)
{
    if (c_sz == 0) {                       /* leaf */
        if (depth > best_size) {
            best_size = depth;
            memcpy(best, stack_v, depth * sizeof(int));
        }
        return;
    }

    /* colour bound */
    int ub = depth + greedy_colour(cand, c_sz);
    if (ub <= best_size) return;

    /* choose pivot = first set bit */
    int pivot = -1;
    for (int w = 0; w < W && pivot == -1; ++w)
        if (cand[w]) pivot = w * 64 + __builtin_ctzll(cand[w]);

    WORD ext[W]; memcpy(ext, cand, W * sizeof(WORD));

    while (c_sz) {
        /* select vertex with smallest colour to branch first */
        int v = -1, best_col = 1 << 30;

        for (int w = 0; w < W; ++w) {
            WORD bits = ext[w];
            while (bits) {
                int x = w * 64 + __builtin_ctzll(bits);
                if (colour[x] < best_col) { best_col = colour[x]; v = x; }
                bits &= bits - 1;
            }
        }
        if (v == -1) break;               /* should not happen */

        /* include v */
        stack_v[depth] = v;

        WORD newC[W]; int new_sz = 0;
        WORD *row = ROW(v);
        for (int w = 0; w < W; ++w) {
            newC[w] = cand[w] & row[w];
            new_sz += POPCNT(newC[w]);
        }
        if (depth + 1 + new_sz > best_size)
            expand(depth + 1, newC, new_sz);

        /* remove v from cand/ext */
        WORD mask = ~((WORD)1ULL << (v % 64));
        cand[v / 64] &= mask;
        ext [v / 64] &= mask;
        --c_sz;
    }
}

/* ------------------------- I/O & driver ------------------------- */
static void read_dimacs(const char *fname)
{
    FILE *fp = fopen(fname, "r");
    if (!fp) { perror("open"); exit(1); }

    char line[256];
    int m = 0;
    while (fgets(line, sizeof line, fp)) {
        if (line[0] == 'p')
            sscanf(line, "p edge %d %d", &N, &m);
    }
    rewind(fp);

    W = (N + 63) / 64;
    adj = calloc((size_t)N * W, sizeof(WORD));
    best     = malloc(N * sizeof(int));
    stack_v  = malloc(N * sizeof(int));

    while (fgets(line, sizeof line, fp)) {
        if (line[0] == 'e') {
            int u, v; sscanf(line, "e %d %d", &u, &v);
            --u; --v;
            ROW(u)[v / 64] |= (WORD)1 << (v % 64);
            ROW(v)[u / 64] |= (WORD)1 << (u % 64);
        }
    }
    fclose(fp);
}

int main(int argc, char **argv)
{
    if (argc != 2) {
        fprintf(stderr, "Usage: %s graph.dimacs\n", argv[0]);
        return 1;
    }
    read_dimacs(argv[1]);

    /* initial candidate = all vertices */
    WORD *cand = calloc(W, sizeof(WORD));
    for (int v = 0; v < N; ++v) cand[v / 64] |= (WORD)1 << (v % 64);

    expand(0, cand, N);

    printf("Maximum clique size = %d\nVertices:", best_size);
    for (int i = 0; i < best_size; ++i) printf(" %d", best[i] + 1);
    puts("");

    return 0;
}
