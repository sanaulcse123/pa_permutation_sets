import itertools

def lcs_length(a, b):
    n = len(a)
    dp = [[0]*(n+1) for _ in range(n+1)]
    for i in range(n):
        for j in range(n):
            if a[i] == b[j]:
                dp[i+1][j+1] = dp[i][j] + 1
            else:
                dp[i+1][j+1] = max(dp[i][j+1], dp[i+1][j])
    return dp[n][n]

def lev_distance(a, b):
    n = len(a)
    return 2 * (n - lcs_length(a, b))

def generate_permutations(n):
    return list(itertools.permutations(range(1, n+1)))

def save_graph(n, d):
    perms = generate_permutations(n)
    filename = f"graph_n{n}_d{d}.txt"
    with open(filename, "w") as f:
        for i, perm_i in enumerate(perms):
            neighbors = []
            for j, perm_j in enumerate(perms):
                if i != j and lev_distance(perm_i, perm_j) >= d:
                    neighbors.append(str(j))
            f.write(f"{i}: {' '.join(neighbors)}\n")
    print(f"Graph saved to {filename} for n={n}, d={d}")

if __name__ == "__main__":
    for n in range(3, 8):         # n = 3 to 7
        for d in range(2, 11):    # d = 2 to 10
            save_graph(n, d)
