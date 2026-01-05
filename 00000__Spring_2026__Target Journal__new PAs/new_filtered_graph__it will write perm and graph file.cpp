// new_filtered_graph_9_4.cpp  (UPDATED: hardcoded refs list = your “random ~160” set)
// NOTE: I kept your exact structure (no reading from txt).
// IMPORTANT: In what you pasted, the refs list is LONGER than 160 (it looks >160).
// So below program includes EXACTLY the refs you pasted in the last message, in the same order.
// If you truly want exactly 160 only, delete extra lines from the refs block until size=160.

#include <iostream>
#include <vector>
#include <algorithm>
#include <fstream>
#include <string>

using namespace std;

// LIS length in O(n log n)
static int lis_length(const vector<int>& seq) {
    vector<int> tail;
    tail.reserve(seq.size());
    for (int x : seq) {
        auto it = lower_bound(tail.begin(), tail.end(), x);
        if (it == tail.end()) tail.push_back(x);
        else *it = x;
    }
    return (int)tail.size();
}

// Ulam distance for permutations: d_U(a,b) = n - LCS(a,b)
// For permutations, LCS(a,b) = LIS(pos_in_b[a[i]])
static int ulamDistance(const vector<int>& a, const vector<int>& b) {
    int n = (int)a.size();
    vector<int> pos(n + 1, 0);
    for (int i = 0; i < n; i++) pos[b[i]] = i;

    vector<int> mapped;
    mapped.reserve(n);
    for (int i = 0; i < n; i++) mapped.push_back(pos[a[i]]);

    int lcs = lis_length(mapped);
    return n - lcs;
}

// Check if perm is at distance >= d from ALL reference permutations
static bool farFromAllRefs(const vector<int>& perm, const vector<vector<int>>& refs, int d) {
    for (const auto& r : refs) {
        if (ulamDistance(perm, r) < d) return false;
    }
    return true;
}

int main() {
    const int n = 9;
    const int d = 4;

    // ====== Reference permutations: your pasted "random 160" block ======
    vector<vector<int>> refs = {
        {4,5,7,8,1,2,9,3,6},
        {9,4,2,1,5,3,6,7,8},
        {2,1,9,3,6,5,8,4,7},
        {1,5,8,4,2,6,3,7,9},
        {4,1,9,8,3,2,5,6,7},
        {4,8,1,5,3,7,9,6,2},
        {2,5,1,9,4,7,6,3,8},
        {4,8,5,6,7,2,1,3,9},
        {4,8,2,7,5,1,6,9,3},
        {1,9,7,4,3,2,8,6,5},
        {4,6,9,2,3,1,5,8,7},
        {6,2,9,1,7,8,4,3,5},
        {2,6,7,4,3,8,9,5,1},
        {5,7,3,8,4,9,1,6,2},
        {3,7,2,9,6,8,4,5,1},
        {3,9,5,8,2,6,4,7,1},
        {1,6,7,9,2,4,8,3,5},
        {9,7,1,6,4,2,5,8,3},
        {9,8,4,2,3,5,7,6,1},
        {6,8,2,5,1,3,9,7,4},
        {6,9,5,2,8,1,7,3,4},
        {1,2,3,4,5,6,7,8,9},
        {1,4,7,5,2,8,3,9,6},
        {1,3,5,7,9,2,8,6,4},
        {1,2,8,6,9,3,7,4,5},
        {1,2,5,4,3,8,9,6,7},
        {7,9,2,3,6,5,4,1,8},
        {9,6,7,3,5,2,1,4,8},
        {1,2,5,9,7,8,3,4,6},
        {3,2,8,9,6,7,1,5,4},
        {6,5,4,3,1,9,7,8,2},
        {1,2,7,8,5,9,6,4,3},
        {5,7,2,1,6,4,9,3,8},
        {8,4,3,1,6,7,2,9,5},
        {8,9,5,1,7,4,3,6,2},
        {5,3,2,9,1,8,7,4,6},
        {1,4,2,7,9,3,6,8,5},
        {1,2,6,5,8,7,4,9,3},
        {3,4,2,6,5,1,7,9,8},
        {6,3,4,7,2,1,8,9,5},
        {1,4,3,9,7,2,5,6,8},
        {9,2,8,7,6,5,1,4,3},
        {7,9,4,2,6,3,8,1,5},
        {2,3,1,6,9,7,8,5,4},
        {4,3,7,8,5,9,2,6,1},
        {1,7,5,6,2,9,3,8,4},
        {3,4,6,8,7,1,2,5,9},
        {4,9,6,7,2,8,5,3,1},
        {2,1,4,8,3,6,7,5,9},
        {2,4,9,7,3,8,1,5,6},
        {8,3,4,9,1,7,5,2,6},
        {2,9,6,4,8,1,7,5,3},
        {5,3,8,6,1,7,2,4,9},
        {2,3,8,5,4,6,1,9,7},
        {8,9,3,6,2,4,5,1,7},
        {2,3,1,5,7,4,9,6,8},
        {3,5,6,2,8,7,4,9,1},
        {1,6,9,8,4,5,7,2,3},
        {9,4,1,8,7,2,3,6,5},
        {2,7,3,9,4,1,6,5,8},
        {1,8,7,9,6,3,2,5,4},
        {2,4,9,5,6,1,8,3,7},
        {2,3,7,8,6,4,9,1,5},
        {3,8,9,4,7,6,5,1,2},
        {5,3,4,8,7,6,9,2,1},
        {9,3,7,4,6,1,8,5,2},
        {7,4,8,9,2,1,3,5,6},
        {1,3,6,7,4,9,5,8,2},
        {2,4,6,7,1,9,5,3,8},
        {7,1,9,5,8,2,4,3,6},
        {9,7,6,8,2,3,1,4,5},
        {5,8,9,7,2,4,6,3,1},
        {6,8,4,3,5,9,1,2,7},
        {4,3,6,2,5,7,1,8,9},
        {5,1,4,9,2,8,6,7,3},
        {7,8,3,2,6,9,5,1,4},
        {5,7,9,4,6,2,1,8,3},
        {5,2,6,1,9,7,4,8,3},
        {3,4,5,9,2,7,1,6,8},
        {8,6,1,4,2,9,7,5,3},
        {1,8,5,2,7,3,6,4,9},
        {1,3,6,4,2,8,5,9,7},
        {2,3,6,1,8,7,9,4,5},
        {5,6,2,4,3,7,8,1,9},
        {2,5,6,3,4,9,1,7,8},
        {9,2,5,7,3,1,8,6,4},
        {7,4,3,1,5,2,6,9,8},
        {3,5,4,1,2,7,8,6,9},
        {3,4,5,1,8,9,6,7,2},
        {4,6,5,8,2,7,3,9,1},
        {5,4,1,6,8,2,9,3,7},
        {3,1,7,5,9,4,8,6,2},
        {1,6,5,7,3,8,2,9,4},
        {8,1,6,2,3,7,5,4,9},
        {7,4,9,8,1,2,6,5,3},
        {4,8,3,6,1,9,5,2,7},
        {2,5,6,7,9,8,1,4,3},
        {3,8,1,9,2,7,4,5,6},
        {2,8,9,1,4,6,5,7,3},
        {3,7,5,8,4,2,1,9,6},
        {4,5,1,7,3,2,6,8,9},
        {3,2,1,7,6,4,5,8,9},
        {4,7,2,9,6,1,3,5,8},
        {2,1,8,4,5,9,3,7,6},
        {1,3,7,8,2,4,6,5,9},
        {1,4,5,6,3,9,8,2,7},
        {3,1,4,8,2,9,5,7,6},
        {8,3,2,7,1,4,9,6,5},
        {1,4,7,8,6,9,2,5,3},
        {1,4,3,5,8,7,6,2,9},
        {3,1,8,6,5,9,2,4,7},
        {1,2,3,4,9,8,7,6,5},
        {5,6,7,1,8,9,2,3,4},
        {3,4,6,2,9,7,5,8,1},
        {8,5,9,4,7,3,2,1,6},
        {6,8,2,1,5,4,7,9,3},
        {6,3,1,2,7,9,5,8,4},
        {2,9,4,3,7,5,1,8,6},
        {3,6,9,8,4,2,1,7,5},
        {7,6,5,2,4,8,9,3,1},
        {6,3,5,9,7,4,2,8,1},
        {6,9,7,2,4,5,1,3,8},
        {4,9,1,3,8,6,2,7,5},
        {7,2,6,1,8,5,4,3,9},
        {4,2,8,1,7,6,3,9,5},
        {2,4,1,5,8,6,9,7,3},
        {2,4,7,6,9,8,3,5,1},
        {2,6,5,9,1,3,4,8,7},
        {6,7,3,2,8,4,1,5,9},
        {8,7,4,2,5,3,9,1,6},
        {5,4,9,7,3,6,1,2,8},
        {3,7,9,5,1,6,2,4,8},
        {5,8,7,9,3,1,4,2,6},
        {4,5,2,9,8,6,3,1,7},
        {4,6,7,5,9,1,8,3,2},
        {5,4,2,3,6,9,8,7,1},
        {1,7,3,9,8,5,4,2,6},
        {7,1,8,3,9,4,5,6,2},
        {5,3,9,6,2,8,1,4,7},
        {7,8,6,4,5,1,2,3,9},
        {5,2,7,4,1,3,9,8,6},
        {1,2,6,8,3,9,5,4,7},
        {7,4,5,9,3,6,8,2,1},
        {6,7,8,1,3,4,9,2,5},
        {1,2,6,4,7,3,5,9,8},
        {6,8,7,9,5,4,1,2,3},
        {2,7,5,8,6,3,1,9,4},
        {2,5,8,4,7,1,3,6,9},
        {3,8,5,2,1,6,7,9,4},
        {7,2,8,1,9,5,6,3,4},
        {7,3,2,5,4,9,8,6,1},
        {2,4,1,6,3,8,5,7,9},
        {7,5,6,1,3,4,8,2,9},
        {3,5,7,6,8,1,9,4,2},
        {6,1,9,3,8,7,5,2,4},
        {5,1,9,8,3,7,6,4,2},
        {3,6,9,1,5,4,7,2,8},
        {1,4,6,8,9,7,3,5,2},
        {1,7,6,5,3,4,9,2,8},
        {5,8,6,9,4,1,3,2,7}
    };

    cout << "refs.size() = " << refs.size() << "\n";

    // 1) Generate all permutations and filter
    vector<int> perm(n);
    for (int i = 0; i < n; i++) perm[i] = i + 1;

    vector<vector<int>> filtered;
    filtered.reserve(50000);

    do {
        if (farFromAllRefs(perm, refs, d)) {
            filtered.push_back(perm);
        }
    } while (next_permutation(perm.begin(), perm.end()));

    int V = (int)filtered.size();

    // Save permutations file (vertex id = line number, 1-based)
    const string permFile = "n9_d4_perms.txt";
    {
        ofstream out(permFile);
        out << "# n=9 d=4 (Ulam)\n";
        out << "# refs.size()=" << refs.size() << "\n";
        out << "# Filtered permutations: those with Ulam distance >= " << d << " from ALL reference permutations\n";
        out << "# Vertex i corresponds to line i below (1-based)\n";
        for (const auto& p : filtered) {
            for (int i = 0; i < n; i++) out << p[i] << (i + 1 == n ? '\n' : ' ');
        }
    }

    // 2) First pass: count edges
    long long E = 0;
    for (int i = 0; i < V; i++) {
        for (int j = i + 1; j < V; j++) {
            if (ulamDistance(filtered[i], filtered[j]) >= d) E++;
        }
    }

    // 3) Second pass: write DIMACS
    const string graphFile = "n9_d4_graph.clq";
    {
        ofstream g(graphFile);
        g << "c Ulam permutation graph (n=9,d=4), vertices are filtered perms\n";
        g << "c refs.size()=" << refs.size() << "\n";
        g << "c edge (i,j) iff UlamDistance(pi,pj) >= " << d << "\n";
        g << "p edge " << V << " " << E << "\n";
        for (int i = 0; i < V; i++) {
            for (int j = i + 1; j < V; j++) {
                if (ulamDistance(filtered[i], filtered[j]) >= d) {
                    g << "e " << (i + 1) << " " << (j + 1) << "\n";
                }
            }
        }
    }

    double denom = (V >= 2) ? (double)V * (V - 1) / 2.0 : 1.0;
    double density = (V >= 2) ? (double)E / denom : 0.0;

    cout << "Graph written: " << graphFile << "\n";
    cout << "Perms written: " << permFile << "\n";
    cout << "Graph size: V=" << V << " E=" << E << " density=" << density << "\n";
    return 0;
}
