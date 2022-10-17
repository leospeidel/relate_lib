#!/usr/bin/env python3

import argparse
import tskit
import matplotlib.pyplot as plt
import numpy as np
import itertools

# ----------------------------------------- #

parser = argparse.ArgumentParser(description="Compare branch statistics across windows from two tree sequences that span the same contig.")
parser.add_argument("--constrained-tree-sequence", type=str, help="A tree seqeuence with equivalent nodes collapsed across marginal trees and constrained node times.")
parser.add_argument("--uncompressed-tree-sequence", type=str, help="A tree sequence with separate nodes across marginal trees")
parser.add_argument("--check-mutations", action="store_true", help="Check that mutation tables ")
parser.add_argument("--output-prefix", type=str, help="Output prefix for plots")
parser.add_argument("--grid-size", type=int, default=251, help="Number of positions across which to compare tree sequences")
args = parser.parse_args()

# ----------------------------------------- #

def tree_tmrca(ts, grid):
    # TODO: should take averages over windows, not point evaluations
    out = np.empty(len(grid)-1)
    tree = ts.first()
    for i, coord in enumerate(grid[:-1]):
        tree.seek(coord)
        out[i] = tree.time(tree.root)
    return out


def pairwise_tmrca(ts, grid):
    out = np.zeros(int(ts.num_samples*(ts.num_samples-1)/2))
    tree = ts.first()
    for i, coord in enumerate(grid[:-1]):
        tree.seek(coord)
        tmp = []
        for j, k in itertools.combinations(np.arange(ts.num_samples), 2):
            tmp.append(tree.tmrca(j, k))
        out += np.array(tmp)
    return out
       

# ----------------------------------------- #

if __name__ == "__main__":
    ts1 = tskit.load(args.uncompressed_tree_sequence)
    ts2 = tskit.load(args.constrained_tree_sequence)
    
    assert ts1.sequence_length == ts2.sequence_length
    assert ts1.num_samples == ts2.num_samples

    if args.check_mutations:
        # TODO: could be a bit more comprehensive
        assert ts1.num_mutations == ts2.num_mutations
        assert ts1.num_sites == ts2.num_sites
        assert np.array_equal(ts1.tables.sites.position, ts2.tables.sites.position)
        assert np.array_equal(ts1.tables.mutations.site, ts2.tables.mutations.site)
    
    windows = np.linspace(0, ts1.sequence_length, args.grid_size)
    
    div1 = ts1.diversity(sample_sets=[[i for i in range(ts1.num_samples)]], windows=windows, mode='branch')
    div2 = ts2.diversity(sample_sets=[[i for i in range(ts2.num_samples)]], windows=windows, mode='branch')
    tajd1 = ts1.Tajimas_D(sample_sets=[[i for i in range(ts1.num_samples)]], windows=windows, mode='branch')
    tajd2 = ts2.Tajimas_D(sample_sets=[[i for i in range(ts2.num_samples)]], windows=windows, mode='branch')
    tree_tmrca1 = tree_tmrca(ts1, windows)
    tree_tmrca2 = tree_tmrca(ts2, windows)
    pair_tmrca1 = np.log10(pairwise_tmrca(ts1, windows))
    pair_tmrca2 = np.log10(pairwise_tmrca(ts2, windows))
    
    fig, ax = plt.subplots(nrows=2, ncols=2, constrained_layout=True)
    for title, pair, coord in zip(
        ["Diversity (branch)", "Tajima's D (branch)", "log10 mean pair TMRCA", "Root time"],
        [(div1,div2),(tajd1,tajd2),(pair_tmrca1,pair_tmrca2),(tree_tmrca1,tree_tmrca2)],
        [(0,0),(0,1),(1,0),(1,1)],
        ):
        i, j = coord
        a, b = pair
        ax[i,j].scatter(a, b, s=10)
        ax[i,j].set_title(title)
        ax[i,j].axline((np.mean(a),np.mean(b)), slope=1, color='black', ls='dashed')
    fig.supxlabel("Values in uncompressed tree sequence")
    fig.supylabel("Values in compressed+constrained tree sequence")
    plt.savefig(args.output_prefix + ".png")

