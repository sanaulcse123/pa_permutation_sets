#include <iostream>
#include <vector>
#include <algorithm>
#include <fstream>

using namespace std;

// Compute LCS between two permutations
int computeLCS(const vector<int>& a, const vector<int>& b) {
    int n = a.size();
    vector<vector<int>> dp(n+1, vector<int>(n+1, 0));
    for(int i = 1; i <= n; ++i) {
        for(int j = 1; j <= n; ++j) {
            if(a[i-1] == b[j-1])
                dp[i][j] = 1 + dp[i-1][j-1];
            else
                dp[i][j] = max(dp[i-1][j], dp[i][j-1]);
        }
    }
    return dp[n][n];
}

// Compute Levenshtein distance using LCS-based relation
int levenshteinDistance(const vector<int>& a, const vector<int>& b) {
    int lcs = computeLCS(a, b);
    return 2 * (a.size() - lcs);
}

int main() {
    int n, d;
    cout << "Enter n (length of permutations): ";
    cin >> n;
    cout << "Enter d (Levenshtein distance threshold): ";
    cin >> d;

    vector<int> identity(n), perm(n);
    for (int i = 0; i < n; ++i)
        identity[i] = perm[i] = i + 1;

    vector<vector<int>> filtered; // Store filtered permutations

    // Step 1: Generate filtered permutations
    do {
        if (levenshteinDistance(perm, identity) >= d)
            filtered.push_back(perm);
    } while (next_permutation(perm.begin(), perm.end()));

    int V = filtered.size();
    vector<pair<int, int>> edges;

    // Step 2: Build edges based on pairwise Levenshtein distance
    for (int i = 0; i < V; ++i) {
        for (int j = i + 1; j < V; ++j) {
            if (levenshteinDistance(filtered[i], filtered[j]) >= d)
                edges.emplace_back(i + 1, j + 1); // 1-based indexing for DIMACS
        }
    }

    int E = edges.size();

    // Step 3: Write .clq graph file
    string filename = "graph_n" + to_string(n) + "_d" + to_string(d) + "_filtered.clq";
    ofstream fout(filename);
    fout << "c Generated Levenshtein permutation graph\n";
    fout << "p edge " << V << " " << E << "\n";
    for (auto& e : edges)
        fout << "e " << e.first << " " << e.second << "\n";
    fout.close();

    // Step 4: Also save filtered permutations for reference
    ofstream permout("filtered_permutations_n" + to_string(n) + "_d" + to_string(d) + ".txt");
    for (auto& p : filtered) {
        for (int i = 0; i < n; ++i)
            permout << p[i] << (i == n - 1 ? "\n" : " ");
    }
    permout.close();

    cout << "Graph written to " << filename << " with " << V << " vertices and " << E << " edges.\n";
    return 0;
}
