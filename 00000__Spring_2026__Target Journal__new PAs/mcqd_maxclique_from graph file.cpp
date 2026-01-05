#include <fstream>
#include <iostream>
#include <set>
#include <string.h>
#include <map>
#include <assert.h>
#include <vector>
#include <string>

#include "mcqd.h"

using namespace std;

void read_dimacs(const string& name, bool** &conn, int &size) {
  ifstream f(name.c_str());
  string buffer;
  assert(f.is_open());

  set<int> v;
  multimap<int,int> e;
  int V_from_p = 0;

  while (getline(f, buffer)) {
    if (buffer.size() == 0) continue;

    if (buffer[0] == 'p') {
      // p edge V E
      char pch[16], edgech[16];
      int Vtmp = 0, Etmp = 0;
      if (sscanf(buffer.c_str(), "%15s %15s %d %d", pch, edgech, &Vtmp, &Etmp) == 4) {
        V_from_p = Vtmp;
      }
    }

    if (buffer[0] == 'e') {
      int vi, vj;
      sscanf(buffer.c_str(), "%*c %d %d", &vi, &vj);
      v.insert(vi);
      v.insert(vj);
      e.insert(make_pair(vi, vj));
    }
  }
  f.close();

  // allocate V+1 because vertices are 1..V
  if (V_from_p > 0) size = V_from_p + 1;
  else if (!v.empty()) size = *v.rbegin() + 1;
  else size = 1;

  conn = new bool*[size];
  for (int i = 0; i < size; i++) {
    conn[i] = new bool[size];
    memset(conn[i], 0, size * sizeof(bool));
  }

  for (auto it = e.begin(); it != e.end(); ++it) {
    conn[it->first][it->second] = true;
    conn[it->second][it->first] = true;
  }

  int V_print = (size > 0 ? size - 1 : 0);
  double denom = (V_print >= 2) ? (double)V_print * (V_print - 1) / 2.0 : 1.0;
  double density = (V_print >= 2) ? (double)e.size() / denom : 0.0;
  cout << "|V|=" << V_print << " |E|=" << e.size() << " density=" << density << "\n";
}

// Reads permutations file created by your generator.
// Returns lines[vertex_id] = permutation line (1-based). lines[0] is dummy.
vector<string> read_perm_file_1based(const string& permFile) {
  ifstream fin(permFile.c_str());
  if (!fin.is_open()) {
    cerr << "Error: cannot open permutation file: " << permFile << "\n";
    exit(1);
  }

  vector<string> lines;
  lines.push_back(""); // dummy for 1-based indexing

  string line;
  while (getline(fin, line)) {
    if (line.size() == 0) continue;
    if (line[0] == '#') continue;
    lines.push_back(line);
  }
  fin.close();
  return lines;
}

int main(int argc, char *argv[]) {
  // usage:
  //   ./mcqd_save_clique_perm graph.clq filtered_permutations.txt
  assert(argc == 3);

  string graphFile = argv[1];
  string permFile  = argv[2];

  cout << "graph: " << graphFile << "\n";
  cout << "perms: " << permFile  << "\n";

  bool **conn = nullptr;
  int size = 0;
  read_dimacs(graphFile, conn, size);

  vector<string> permLines = read_perm_file_1based(permFile);

  // Run dynamic max clique; checkpoint file best_clique_perms.txt will be updated on improvements
  Maxclique md(conn, size, 0.025);
  md.set_perm_lines(&permLines);

  int *qmax = nullptr;
  int qsize = 0;

  md.mcqdyn(qmax, qsize);

  cout << "DONE. Final best clique size = " << qsize << "\n";
  cout << "Best clique permutations saved (and updated during run) in: best_clique_perms.txt\n";

  delete [] qmax;

  for (int i = 0; i < size; i++) delete [] conn[i];
  delete [] conn;

  return 0;
}
