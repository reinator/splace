"""Command line for SPLACE: single-threshold GMYC (Pons et al. 2006) with
pGMYC (Zhang et al. 2013), in its Python 3 port by iTaxoTools (GMYC.py next
to this file, GPL-3.0, commit 48445d1, unchanged).

The command-line block of that GMYC.py cannot run (its gmyc() writes files
named with `self.unique` outside a class), so this script calls the same
functions in the same order as gmyc() and prints the result as JSON. Only
the computation is used: the ete3 drawing classes (NodeStyle, which need
PyQt) and matplotlib are replaced by stand-ins when they are missing, since
SPLACE draws the species on its own tree viewer.

    python splace_gmyc.py ULTRAMETRIC_TREE [--pvalue 0.01] [--out result.json]
"""
import argparse
import json
import os
import sys
import types


def _stub_missing():
    try:
        import matplotlib  # noqa: F401
        matplotlib.use("Agg")
        import matplotlib.pyplot  # noqa: F401
    except Exception:
        mpl = types.ModuleType("matplotlib")
        mpl.use = lambda *a, **k: None
        plt = types.ModuleType("matplotlib.pyplot")
        mpl.pyplot = plt
        sys.modules["matplotlib"] = mpl
        sys.modules["matplotlib.pyplot"] = plt
    try:
        import nexus  # noqa: F401
    except Exception:
        nx = types.ModuleType("nexus")
        nx.NexusReader = None
        sys.modules["nexus"] = nx
    import ete3
    if not hasattr(ete3, "NodeStyle"):
        ete3.NodeStyle = dict


def run(tree_text, pvalue):
    _stub_missing()
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    import GMYC as g
    from scipy.optimize import fmin_l_bfgs_b

    # As GMYC.py's gmyc(): every internal node as the threshold, the best
    # likelihood, the null model and the likelihood-ratio test (2 df).
    min_change, max_iters = 0.1, 100
    best_llh, best_num_spe, best_node = float("-inf"), -1, None
    utree = g.um_tree(tree_text)
    for tnode in utree.nodes:
        wt_list, num_spe = utree.get_waiting_times(threshold_node=tnode)
        tt = g.tree_time(wt_list, num_spe)
        last_llh, change, cnt = float("-inf"), float("inf"), 0
        while change > min_change and cnt < max_iters:
            cnt += 1
            para, nn, cc = fmin_l_bfgs_b(g.tar_fun, [1, 1], args=tuple([tt]), bounds=[[0, 10], [0, 10]], approx_grad=True)
            tt.update(para[0], para[1])
            logl = tt.sum_llh()
            change = abs(logl - last_llh)
            last_llh = logl
        final_llh = tt.sum_llh()
        if final_llh > best_llh:
            best_llh, best_num_spe, best_node = final_llh, num_spe, tnode
    null_logl = g.optimize_null_model(utree)
    utree.get_waiting_times(threshold_node=best_node)
    species = [[leaf.name for leaf in sp] for sp in utree.species_list]
    p = g.lh_ratio_test(null_llh=null_logl, llh=best_llh, df=2).get_p_value()
    if p >= pvalue:
        species = [[leaf.name for leaf in utree.tree.get_leaves()]]
    # Plain Python numbers (scipy gives numpy types, which json refuses).
    return {
        "method": "GMYC (single threshold), pGMYC",
        "best_logl": float(best_llh),
        "null_logl": float(null_logl),
        "pvalue": float(p),
        "threshold_species": int(best_num_spe),
        "lrt_passed": bool(p < pvalue),
        "species": species,
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("tree", help="rooted ultrametric tree (Newick file)")
    ap.add_argument("--pvalue", type=float, default=0.01, help="p-value of the likelihood-ratio test (default 0.01)")
    ap.add_argument("--out", help="JSON file for the result (default: standard output)")
    ap.add_argument("--version", action="version", version="splace-gmyc 1.0 (pGMYC 1.1, iTaxoTools 48445d1)")
    a = ap.parse_args()
    with open(a.tree, encoding="utf-8") as f:
        text = f.read().strip()
    result = run(text, a.pvalue)
    out = json.dumps(result, indent=1)
    if a.out:
        with open(a.out, "w", encoding="utf-8") as f:
            f.write(out + "\n")
    print(out)
    print("Final number of estimated species by GMYC: %d" % len(result["species"]), file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
