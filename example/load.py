import tskit

ts = tskit.load("./example.trees")

for tree in ts.trees():
    for n in tree.nodes():
        meta = [int(m) for m in ts.node(n).metadata.split()]
        print(n,meta)

    break #break after first tree, remove if continue

