 How to Run the Project( for n=3,4 and 5)
==========================

1. First, run `levenshtein_graph.cpp`.
   This will generate:
   - `permutations_nX.txt` → contains all permutations of size `n`
   - `graph_nX_dY.dimacs` → DIMACS-format graph with edges where Levenshtein distance ≥ `d`

2. Next, run `maxclique_finder.cpp`.
   This will generate:
   - `clique_nX_dY.clq` → vertex indices forming the maximum clique
   - `clique_nX_dY.perm` → actual permutations corresponding to the clique

-----------------------------------------

 Description of Output Files
===============================

| File Name               | Description                                      |
|------------------------|--------------------------------------------------|
| permutations_nX.txt    | All permutations of length `n`                  |
| graph_nX_dY.dimacs     | Graph edges for pairs with Levenshtein ≥ `d`    |
| clique_nX_dY.clq       | Indices of permutations forming max clique      |
| clique_nX_dY.perm      | Actual permutations in the max clique           |
