#!/usr/bin/env python3

import scipy.sparse
import numpy as np
import cvxpy as cp
import array

import tskit
#import tsdate

import argparse

# ----------------------------------------- #

parser = argparse.ArgumentParser(
    description=
        "Given a tree sequence with a set of node times that result in negative "
        "branch lengths, find a set of valid node times such that the relative "
        "squared error is minimal. Original (invalid) node times should be stored "
        "as double precision in the node metadata of the input tree sequence."
)
parser.add_argument("--tree-sequence", type=str, help="Tree sequence with unconstrained node times in node metadata (e.g. node times such that some branch lengths are negative)")
parser.add_argument("--output", type=str, help="Output tree sequence with constrained node times (or table collection with approximate node times if constraint is not met)") 
parser.add_argument("--skip-optimization", action="store_true", help="Skip minimization of relative square error")
parser.add_argument("--clear-mutations", action="store_true", help="Clean mutations and sites from tree sequence")
parser.add_argument("--plot-constrained-vs-unconstrained", action="store_true", help="Save plot of log10 constrained/unconstrained node times to <output>.png")
parser.add_argument("--minimum-branch-length", type=float, default=1e-3, help="Minimum allowed branch length")
parser.add_argument("--eps", type=float, default=1e-10, help="Optimization convergence tolerance")
parser.add_argument("--max-iters", type=int, default=10000, help="Maximum iterations for convex solver")
parser.add_argument("--verbose", action="store_true", help="Print output from convex solver")
parser.add_argument("--starting-state", type=str, help="Table collection with node times to use as starting point")
args = parser.parse_args()

# ----------------------------------------- #

def node_age_to_branch_length(ts):
    """
    Sparse matrix mapping node age to branch lengths.
    """
    edge_id = np.arange(ts.num_edges, dtype=np.int32)
    row_idx = np.tile(edge_id, 2)
    col_idx = np.concatenate([
        ts.tables.edges.parent - ts.num_samples,
        ts.tables.edges.child - ts.num_samples,
    ])
    value = np.repeat([1, -1], ts.num_edges)
    out_of_bounds = col_idx < 0
    row_idx = row_idx[~out_of_bounds]
    col_idx = col_idx[~out_of_bounds]
    value = value[~out_of_bounds]
    return scipy.sparse.coo_array(
        (value, (row_idx, col_idx)),
        shape=(ts.num_edges, ts.num_nodes-ts.num_samples),
    ).tocsr()


def constrain_ages_quadprog(ts, node_age, min_len, starting_state=None, **kwargs):
    """
    Use quadratic programming to find the solution to,

    min \sum_i (t_i - \hat{t}_i)^2 / \hat{t}_i
    s.t. A t_i >= eps, t_i >= eps

    where t_i is the age of internal node i, \hat{t}_i is an unconstrained node age (e.g. posterior mean),
    A is a sparse matrix mapping internal node ages to branch lengths, and eps is the minimum possible
    branch length.

    NB: Would have to include an offset for ancient samples
    NB: OSQP is super fast but uses sparse matrix factorization to solve subproblem. 
        SCS with an indirect solver would scale better.
    """
    if starting_state is not None:
        tables = tskit.TableCollection.load(starting_state)
        assert tables.nodes.num_rows == ts.num_nodes
        assert sum(tables.nodes.flags == tskit.NODE_IS_SAMPLE) == ts.num_samples
        xhat = tables.nodes.time[ts.num_samples:]
    else:
        xhat = node_age[ts.num_samples:]
    x = cp.Variable(len(xhat))
    x.value = xhat.copy()
    eps_edge = np.full(ts.num_edges, min_len)
    eps_node = np.full(ts.num_nodes - ts.num_samples, min_len)
    A = node_age_to_branch_length(ts)
    P = scipy.sparse.spdiags(1/np.sqrt(xhat), 0, len(xhat), len(xhat))
    obj = cp.Minimize(cp.sum_squares(P @ x - P @ xhat))
    constr = [A @ x >= eps_edge, x >= eps_node]
    prob = cp.Problem(obj, constr)
    result = prob.solve(**kwargs)
    min_edge = np.min(A @ x.value)
    return np.concatenate([node_age[:ts.num_samples], x.value]), min_edge


# ----------------------------------------- #

if __name__ == "__main__":
    ts = tskit.load(args.tree_sequence)
    node_times = np.array(array.array('d', ts.tables.nodes.metadata.tobytes()))
    if args.verbose:
        print(f"Loaded {args.tree_sequence} ...")
    # TODO: should also be able to load in from unsorted TableCollection with unconstrained node times in tables.nodes.time

    if not np.all(node_times[:ts.num_samples] == 0.0):
        # TODO: This should be easy enough to support, just add a constant vector to cvxpy constraints
        raise ValueError("Currently ancient samples are not supported")

    node_times_constrained = node_times.copy()
    if not args.skip_optimization:
        if args.verbose:
            print(f"Finding constrained branch lengths ...")
        solver_args = {"solver":"OSQP", "eps_rel":args.eps, "eps_abs":args.eps, "max_iter":args.max_iters, "verbose":args.verbose}
        node_times_constrained, min_edge_length = constrain_ages_quadprog(ts, node_times_constrained, args.minimum_branch_length, args.starting_state, **solver_args)

    # TODO: this doesn't work to "brute force" the constraint -- why?
    #node_times_constrained = tsdate.core.constrain_ages_topo(ts, node_times_constrained, eps=args.eps)

    if args.plot_constrained_vs_unconstrained:
        if args.verbose:
            print(f"Saving plot of constrained vs unconstrained times to {args.output + '.png'} ...")
        import matplotlib.pyplot as plt
        x, y = np.log10(node_times[ts.num_samples:]), np.log10(node_times_constrained[ts.num_samples:])
        plt.scatter(x, y, s=10)
        plt.axline((np.mean(x), np.mean(y)), slope=1, color='black', ls='dashed')
        plt.xlabel("log10 node age (unconstrained)")
        plt.ylabel("log10 node age (constrained)")
        plt.savefig(args.output + ".png")

    tables = ts.dump_tables()
    if args.clear_mutations:
        if args.verbose:
            print(f"Removing mutations and sites ... ")
       tables.mutations.clear()
       tables.sites.clear()
    else:
        tables.mutations.time.fill(tskit.UNKNOWN_TIME)
    tables.nodes.time = node_times_constrained

    if min_edge_length < 0.0:
        # Likely hasn't converged, or args.eps is too high
        print(
            f"Saving tables to {args.output} because node times do not meet constraint.\n"
            f"These may be used as a starting point for further optimization, e.g.\n "
            f"  ./add-constrained-node-ages.py --tree-sequence {args.tree_sequence} "
            f"--max-iters 20000 --verbose --starting-state {args.output}"
        )
        tables.dump(args.output)
    else:
        tables.sort()
        tables.tree_sequence().dump(args.output)
        if args.verbose:
            print(f"Saved output tree sequence to {args.output}.")
