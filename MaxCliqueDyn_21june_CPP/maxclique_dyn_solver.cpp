#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <algorithm>
#include <sstream>
using namespace std;

int n;
vector<vector<bool>> adj;
vector<vector<int>> permutations;

// Read DIMACS graph
void read_dimacs(const string& filename) {
    ifstream fin(filename);
    string line;
    int edges;
    while (getline(fin, line)) {
        if (line[0] == 'p') {
            sscanf(line.c_str(), "p edge %d %d", &n, &edges);
            adj.assign(n, vector<bool>(n, false));
        } else if (line[0] == 'e') {
            int u, v;
            sscanf(line.c_str(), "e %d %d", &u, &v);
            adj[u - 1][v - 1] = true;
            adj[v - 1][u - 1] = true;
        }
    }
}

// Read permutations
void read_permutations(const string& filename) {
    ifstream fin(filename);
    string line;
    while (getline(fin, line)) {
        istringstream iss(line);
        vector<int> p;
        int x;
        while (iss >> x) p.push_back(x);
        permutations.push_back(p);
    }
}

// Coloring upper bound heuristic
int greedy_coloring(const vector<int>& vertices) {
    int n = vertices.size();
    vector<int> color(n, -1);
    int max_color = 0;
    for (int i = 0; i < n; ++i) {
        bool used[100] = {};
        for (int j = 0; j < i; ++j)
            if (adj[vertices[i]][vertices[j]])
                used[color[j]] = true;
        int c = 0;
        while (used[c]) c++;
        color[i] = c;
        max_color = max(max_color, c + 1);
    }
    return max_color;
}

// MaxCliqueDyn recursive BnB
void maxclique_dyn(vector<int>& R, vector<int>& C, vector<int>& best, int& best_size) {
    if (C.empty()) {
        if ((int)R.size() > best_size) {
            best = R;
            best_size = R.size();
        }
        return;
    }
    if ((int)R.size() + greedy_coloring(C) <= best_size) return;

    for (size_t i = 0; i < C.size(); ++i) {
        int v = C[i];
        vector<int> R_new = R;
        R_new.push_back(v);

        vector<int> C_new;
        for (size_t j = i + 1; j < C.size(); ++j)
            if (adj[v][C[j]]) C_new.push_back(C[j]);

        maxclique_dyn(R_new, C_new, best, best_size);
    }
}

int main() {
    int perm_n, d;
    cout << "Enter value of n: ";
    cin >> perm_n;
    cout << "Enter value of d: ";
    cin >> d;

    string dimacs_file = "graph_n" + to_string(perm_n) + "_d" + to_string(d) + ".dimacs";
    string perm_file = "permutations_n" + to_string(perm_n) + ".txt";
    string clq_out = "clique_n" + to_string(perm_n) + "_d" + to_string(d) + ".clq";
    string perm_out = "clique_n" + to_string(perm_n) + "_d" + to_string(d) + ".perm";

    read_dimacs(dimacs_file);
    read_permutations(perm_file);

    vector<int> candidates(n), R, best;
    for (int i = 0; i < n; ++i) candidates[i] = i;
    int best_size = 0;

    maxclique_dyn(R, candidates, best, best_size);

    cout << "Max Clique Size: " << best_size << endl;
    cout << "Vertices in Max Clique:\n";
    for (int i : best) cout << i + 1 << " ";
    cout << "\nPermutations:\n";
    for (int i : best) {
        for (int x : permutations[i]) cout << x << " ";
        cout << endl;
    }

    ofstream clqfile(clq_out), permfile(perm_out);
    for (int i : best) clqfile << i + 1 << "\n";
    for (int i : best) {
        for (int x : permutations[i]) permfile << x << " ";
        permfile << "\n";
    }

    cout << "Saved clique indices to: " << clq_out << endl;
    cout << "Saved permutations to: " << perm_out << endl;
    return 0;
}
