## ComputeDistances.py
## shortest path lengths between nodes of the Graph (20 rows x 50 columns)

import pandas as pd
import networkx as nx

# create grilled graph
def CreateGraph(num_rows, num_cols):
    all_nodes = num_cols * num_rows - 1
    G = nx.Graph()
    horizontal_edges = [(i, i + 1) for i in range(num_cols * num_rows) if (i + 1 - num_cols) % num_cols != 0]
    vertical_edges = [(i, i + num_cols) for i in range(num_cols * num_rows) if i + num_cols <= all_nodes]
    all_edges = horizontal_edges + vertical_edges
    G.add_edges_from(all_edges)
    return G


# Create matrix of shortest path lengths between nodes of the Graph (G_grill)
def spl(G_grill):
    distances = []
    for i in G_grill:
        dist = []
        for j in G_grill:
            dist.append(nx.shortest_path_length(G_grill, source=i, target=j))
        distances.append(dist)

    return distances


### MAIN ###

rows = 20
columns = 50

G_grill = CreateGraph(rows, columns)            # grid for rows and random arrangements

distances = spl(G_grill)                        # Matrix of short path lengths from node to node
distances_df = pd.DataFrame(distances)
distances_df.to_csv('distances_grill.csv', index=False, header=False)

