"""
Demo / CLI runner for the monotone tree-drawing algorithms.

A rooted tree is described by a balanced-parenthesis string (see the thesis /
README). Give one with --paren, or generate a random tree with --random N.
Then pick one or more algorithms to draw it.

Examples
--------
    # Random 8-node tree, run every algorithm, show the plots
    python demo.py --random 8 --algo all --show

    # Specific tree from a parenthesis string, save the He & He drawing
    python demo.py --paren "(()())(())" --algo hehe --save out.png
"""

import argparse

import Tree_Operations as tops
import Angelini as ang
import He_He as hh


def build_tree(args):
    if args.paren is not None:
        if not tops.isvalid(args.paren):
            raise SystemExit(f"Not a valid balanced-parenthesis string: {args.paren!r}")
        paren = args.paren
    else:
        paren = tops.random_tree(args.random)
    print(f"Tree (parenthesis): {paren}")
    return tops.paren_to_nxgraph(paren)


def run_algo(name, g, show, save):
    """Run one algorithm and print its grid size. Returns [area, maxx, maxy, pos]."""
    filename = save or f"{name}.png"
    do_save = save is not None
    if name == "bfs":
        res = ang.getGridArea(g, algo="bfs", display=show, save=do_save, filename=filename)
    elif name == "dfs":
        res = ang.getGridArea(g, algo="dfs", display=show, save=do_save, filename=filename)
    elif name == "hehe":
        res = hh.PathDrawAlgorithm(g, display=show, save=do_save, filename=filename)
    else:
        raise ValueError(name)
    area, maxx, maxy = res[0], res[1], res[2]
    print(f"  {name:5s} -> grid {maxx} x {maxy}  (bounding area {area})")
    if do_save:
        print(f"          saved to {filename}")
    return res


def main():
    p = argparse.ArgumentParser(description="Draw a rooted tree with monotone-drawing algorithms.")
    src = p.add_mutually_exclusive_group(required=True)
    src.add_argument("--paren", help="balanced-parenthesis string describing the tree")
    src.add_argument("--random", type=int, metavar="N", help="generate a random tree with N nodes")
    p.add_argument("--algo", choices=["bfs", "dfs", "hehe", "all"], default="all",
                   help="which algorithm(s) to run (default: all)")
    p.add_argument("--show", action="store_true", help="display the drawing(s) in a window")
    p.add_argument("--save", metavar="FILE",
                   help="save the drawing to FILE (a per-algorithm name is used when --algo all)")
    args = p.parse_args()

    g = build_tree(args)
    print(f"Nodes: {g.number_of_nodes()}, edges: {g.number_of_edges()}")

    algos = ["bfs", "dfs", "hehe"] if args.algo == "all" else [args.algo]
    print("Grid sizes:")
    for name in algos:
        if args.save is None:
            save = None
        elif len(algos) > 1:
            # Multiple algos can't share one filename; use a per-algorithm name.
            save = f"{name}.png"
        else:
            save = args.save
        run_algo(name, g, args.show, save)


if __name__ == "__main__":
    main()
