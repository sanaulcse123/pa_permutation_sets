#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <set>
#include <algorithm>
#include <string>

using namespace std;

int max_clique_size = 0;
vector<int> best_clique;
vector<set<int>> graph;

void bron_kerbosch(set<int> R, set<int> P, set<int> X) {
    if (P.empty() && X.empty()) {
        if ((int)R.size() > max_clique_size) {
            max_clique_size = R.size();
            best_clique.assign(R.begin(), R.end());
        }
        return;
    }

    int u = *P.begin();
    set<int> P_without_neighbors;
    for (int v : P) {
        if (graph[u].find(v) == graph[u].end()) {
            P_without_neighbors.insert(v);
        }
    }

    for (int v : P_without_neighbors) {
        set<int> newR = R, newP, newX;
        newR.insert(v);
        for (int n : graph[v]) if (P.find(n) != P.end()) newP.insert(n);
        for (int n : graph[v]) if (X.find(n) != X.end()) newX.insert(n);
        bron_kerbosch(newR, newP, newX);
        P.erase(v);
        X.insert(v);
    }
}

bool read_graph(string filename, int &n) {
    ifstream in(filename);
    if (!in.is_open()) return false;

    string line;
    int u, v;
    n = 0;
    graph.clear();

    while (getline(in, line)) {
        if (line.empty()) continue;
        if (line[0] == 'p') {
            istringstream iss(line);
            string p, edge;
            int nodes, edges;
            iss >> p >> edge >> nodes >> edges;
            graph.resize(nodes + 1);
            n = nodes;
        } else if (line[0] == 'e') {
            istringstream iss(line);
            string e;
            iss >> e >> u >> v;
            graph[u].insert(v);
            graph[v].insert(u);
        }
    }
    return true;
}

void write_clique_file(const vector<int>& clique, const string& filename) {
    ofstream out(filename);
    for (int v : clique) out << v << " ";
    out << "\n";
    out.close();
}

void write_perm_file(const string& permFile, const vector<int>& clique, const string& outFile) {
    ifstream in(permFile);
    vector<string> lines;
    string line;
    while (getline(in, line)) lines.push_back(line);
    in.close();

    ofstream out(outFile);
    for (int idx : clique) {
        if (idx - 1 < lines.size())
            out << lines[idx - 1] << "\n";
    }
    out.close();
}

void run_for(int n, int d) {
    string dimacs_file = "graph_n" + to_string(n) + "_d" + to_string(d) + ".dimacs";
    string perm_file   = "permutations_n" + to_string(n) + ".txt";
    string clq_file    = "clique_n" + to_string(n) + "_d" + to_string(d) + ".clq";
    string perm_out    = "clique_n" + to_string(n) + "_d" + to_string(d) + ".perm";

    int nodes = 0;
    if (!read_graph(dimacs_file, nodes)) {
        cout << " Missing: " << dimacs_file << "\n";
        return;
    }

    set<int> R, P, X;
    for (int i = 1; i <= nodes; ++i) P.insert(i);
    max_clique_size = 0;
    best_clique.clear();
    bron_kerbosch(R, P, X);

    cout << " MaxClique size for n=" << n << ", d=" << d << " is: " << max_clique_size << "\n";
    write_clique_file(best_clique, clq_file);
    write_perm_file(perm_file, best_clique, perm_out);
}

int main() {
    for (int n = 3; n <= 5; ++n) {
        for (int d = 2; d <= 10; ++d) {
            run_for(n, d);
        }
    }
    return 0;
}
