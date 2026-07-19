"""Demo / CLI runner for the monotone tree-drawing algorithms.

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

from __future__ import annotations

import argparse

import networkx as nx

import Angelini as ang
import He_He as hh
import Tree_Operations as tops

# Maps an algorithm name to a callable(graph, **draw_kwargs) -> [area, maxx, maxy, pos].
ALGORITHMS = {
    "bfs": lambda g, **kw: ang.getGridArea(g, algo="bfs", **kw),
    "dfs": lambda g, **kw: ang.getGridArea(g, algo="dfs", **kw),
    "hehe": lambda g, **kw: hh.PathDrawAlgorithm(g, **kw),
}


def build_tree(args: argparse.Namespace) -> nx.Graph:
    """Build the tree from the CLI arguments and echo its parenthesis string."""
    paren = args.paren if args.paren is not None else tops.random_tree(args.random)
    print(f"Tree (parenthesis): {paren}")
    return tops.paren_to_nxgraph(paren)


def run_algo(name: str, g: nx.Graph, show: bool, save: str | None) -> list:
    """Run one algorithm, print its grid size, and return its result list."""
    res = ALGORITHMS[name](g, display=show, save=save is not None, filename=save or f"{name}.png")
    area, maxx, maxy = res[0], res[1], res[2]
    print(f"  {name:5s} -> grid {maxx} x {maxy}  (bounding area {area})")
    if save is not None:
        print(f"          saved to {save}")
    return res


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Draw a rooted tree with monotone-drawing algorithms.")
    src = p.add_mutually_exclusive_group(required=True)
    src.add_argument("--paren", help="balanced-parenthesis string describing the tree")
    src.add_argument("--random", type=int, metavar="N", help="generate a random tree with N nodes")
    p.add_argument(
        "--algo",
        choices=[*ALGORITHMS, "all"],
        default="all",
        help="which algorithm(s) to run (default: all)",
    )
    p.add_argument("--show", action="store_true", help="display the drawing(s) in a window")
    p.add_argument(
        "--save",
        metavar="FILE",
        help="save the drawing to FILE (a per-algorithm name is used when --algo all)",
    )
    return p.parse_args(argv)


def main(argv: list[str] | None = None) -> None:
    args = parse_args(argv)
    algos = list(ALGORITHMS) if args.algo == "all" else [args.algo]

    try:
        g = build_tree(args)
        print(f"Nodes: {g.number_of_nodes()}, edges: {g.number_of_edges()}")

        print("Grid sizes:")
        for name in algos:
            if args.save is None:
                save = None
            elif len(algos) > 1:
                # Multiple algorithms can't share one file; use a per-algorithm name.
                save = f"{name}.png"
            else:
                save = args.save
            run_algo(name, g, args.show, save)
    except (TypeError, ValueError, OSError) as exc:
        raise SystemExit(f"error: {exc}") from exc


if __name__ == "__main__":
    main()
