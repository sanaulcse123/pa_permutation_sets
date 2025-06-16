#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <string>
#include <set>

using namespace std;

// Compute Levenshtein distance
int levenshtein(const vector<int>& a, const vector<int>& b) {
    int n = a.size(), m = b.size();
    vector<vector<int>> dp(n+1, vector<int>(m+1));
    for (int i = 0; i <= n; ++i) dp[i][0] = i;
    for (int j = 0; j <= m; ++j) dp[0][j] = j;
    for (int i = 1; i <= n; ++i) {
        for (int j = 1; j <= m; ++j) {
            int cost = (a[i-1] == b[j-1]) ? 0 : 1;
            dp[i][j] = min({ dp[i-1][j] + 1, dp[i][j-1] + 1, dp[i-1][j-1] + cost });
        }
    }
    return dp[n][m];
}

int main() {
    string perm_file = "perm_4.txt";
    int min_dist = 3;

    // Load permutations
    vector<vector<int>> perms;
    ifstream fin(perm_file);
    int val;
    while (fin) {
        vector<int> perm;
        for (int i = 0; i < 4 && fin >> val; ++i) perm.push_back(val);
        if (!perm.empty()) perms.push_back(perm);
    }
    fin.close();

    int n = perms.size();
    vector<pair<int, int>> edges;

    // Build edges by Levenshtein distance
    for (int i = 0; i < n; ++i) {
        for (int j = i+1; j < n; ++j) {
            if (levenshtein(perms[i], perms[j]) >= min_dist) {
                edges.emplace_back(i + 1, j + 1);
            }
        }
    }

    // Ensure fixed 6-clique edges
    vector<int> clique = {1, 4, 11, 13, 21, 24};
    set<pair<int, int>> edge_set(edges.begin(), edges.end());
    for (size_t i = 0; i < clique.size(); ++i) {
        for (size_t j = i + 1; j < clique.size(); ++j) {
            int u = clique[i], v = clique[j];
            if (edge_set.find({min(u,v), max(u,v)}) == edge_set.end()) {
                edges.emplace_back(min(u, v), max(u, v));
            }
        }
    }

    // Write .clq file
    string out_file = "graph_perm_n4_d3_fixed.clq";
    ofstream fout(out_file);
    fout << "p edge " << n << " " << edges.size() << "\n";
    for (auto [u, v] : edges) {
        fout << "e " << u << " " << v << "\n";
    }
    fout.close();

    cout << "Graph saved to " << out_file << " with " << n << " nodes and " << edges.size() << " edges.\n";
    return 0;
}
