import os
import networkx as nx
from networkx.algorithms.clique import find_cliques

def load_graph(filepath):
    G = nx.Graph()
    with open(filepath, 'r') as file:
        for line in file:
            if ":" not in line:
                continue
            node_str, neighbors_str = line.strip().split(":")
            node = int(node_str.strip())
            neighbors = list(map(int, neighbors_str.strip().split()))
            for neighbor in neighbors:
                G.add_edge(node, neighbor)
    return G

# Loop for n = 3 to 6 and d = 2 to 10
for n in range(3, 7):
    for d in range(2, 11):
        graph_file = f"graph_n{n}_d{d}.txt"
        output_file = f"max_clique_n{n}_d{d}.txt"

        if not os.path.exists(graph_file):
            print(f" Skipping missing file: {graph_file}")
            continue

        try:
            print(f" Processing {graph_file}")
            G = load_graph(graph_file)
            cliques = list(find_cliques(G))

            if not cliques:
                with open(output_file, "w") as f:
                    f.write("Max Clique Size: 0\nNo valid clique found.\n")
                continue

            max_size = max(len(c) for c in cliques)
            max_cliques = [c for c in cliques if len(c) == max_size]

            with open(output_file, "w") as f:
                f.write(f"Max Clique Size: {max_size}\n")
                f.write("Max Clique(s):\n")
                for clique in max_cliques:
                    f.write(" ".join(map(str, clique)) + "\n")

        except Exception as e:
            print(f" Error in {graph_file}: {e}")
            with open(output_file, "w") as f:
                f.write("Error processing graph file.\n")
