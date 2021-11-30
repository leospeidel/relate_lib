import tskit

ts = tskit.load("./example.trees")

#from https://github.com/tskit-dev/tskit/issues/186
edge_map = ts._tree_node_edges()

#for tree, edges in zip(ts.trees(), edge_map):
for tree in ts.trees():
    edges = next(edge_map)
    for n in tree.nodes():
        if edges[n] != -1: #e.g. root node
            e = ts.edge(edges[n])
            meta = [int(m) for m in e.metadata.split()]
            print("Equivalent child node in previous tree: ", meta[0], "; Equivalent child node in next tree: ", meta[1])

    break #break after first tree, remove if continue

