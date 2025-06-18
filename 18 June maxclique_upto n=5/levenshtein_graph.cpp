#include <iostream>
#include <vector>
#include <fstream>
#include <algorithm>
#include <string>

using namespace std;

const int MAXN = 5;

// Levenshtein distance via LCS for permutations: L(a,b) = 2 * (n - LCS)
int levenshtein_lcs(const vector<int>& a, const vector<int>& b, int n) {
    vector<vector<int>> dp(n + 1, vector<int>(n + 1, 0));
    for (int i = 1; i <= n; ++i) {
        for (int j = 1; j <= n; ++j) {
            if (a[i - 1] == b[j - 1])
                dp[i][j] = dp[i - 1][j - 1] + 1;
            else
                dp[i][j] = max(dp[i - 1][j], dp[i][j - 1]);
        }
    }
    return 2 * (n - dp[n][n]);
}

vector<vector<int>> generate_permutations(int n) {
    vector<vector<int>> result;
    vector<int> perm(n);
    for (int i = 0; i < n; ++i) perm[i] = i + 1;
    do {
        result.push_back(perm);
    } while (next_permutation(perm.begin(), perm.end()));
    return result;
}

void save_permutations(int n, const vector<vector<int>>& perms) {
    string filename = "permutations_n" + to_string(n) + ".txt";
    ofstream fout(filename);
    for (const auto& p : perms) {
        for (int val : p) fout << val << " ";
        fout << "\n";
    }
    fout.close();
    cout << " Saved: " << filename << "\n";
}

void generate_graph(int n, int d, const vector<vector<int>>& perms) {
    string filename = "graph_n" + to_string(n) + "_d" + to_string(d) + ".dimacs";
    ofstream fout(filename);

    int V = perms.size();
    int edge_count = 0;
    vector<pair<int, int>> edges;

    for (int i = 0; i < V; ++i) {
        for (int j = i + 1; j < V; ++j) {
            int dist = levenshtein_lcs(perms[i], perms[j], n);
            if (dist >= d) {
                edges.emplace_back(i + 1, j + 1);
                edge_count++;
            }
        }
    }

    fout << "p edge " << V << " " << edge_count << "\n";
    for (auto& e : edges) {
        fout << "e " << e.first << " " << e.second << "\n";
    }

    fout.close();
    cout << " Graph written: " << filename << " with " << edge_count << " edges\n";
}

int main() {
    for (int n = 3; n <= MAXN; ++n) {
        auto perms = generate_permutations(n);
        save_permutations(n, perms);
        for (int d = 2; d <= 10; ++d) {
            generate_graph(n, d, perms);
        }
    }
    return 0;
}
