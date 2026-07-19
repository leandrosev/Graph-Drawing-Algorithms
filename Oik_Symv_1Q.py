"""Oikonomou & Symvonis one-quadrant monotone tree-drawing algorithm.

The single-quadrant variant: the tree is rooted and every node is given an
angular wedge (proportional to its subtree size) inside the first quadrant
``[0, pi/2]``. Edge directions are primitive integer vectors chosen within each
wedge, so the whole drawing lives in one quadrant on a compact ``~n x n`` grid.

This is the same geometric idea as :mod:`Oik_Symv_Balanced` (the four-quadrant
balanced version); here the wedge spans a quarter-turn and every coordinate is
non-negative, which makes it convenient for composing sub-drawings from a chosen
root and origin.

Reference:
    Anargyros Oikonomou, Antonios Symvonis,
    "Simple Compact Monotone Tree Drawings".
    https://arxiv.org/pdf/1708.09653.pdf

"""

from __future__ import annotations

import math

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np

Coord = tuple[int, int]


# ---------------------------------------------------------------------------
# Validation
# ---------------------------------------------------------------------------
def _validate_tree(tree: object, root: object) -> None:
    """Validate that *tree* is a non-empty undirected tree containing *root*.

    Raises
    ------
    TypeError
        If *tree* is not an undirected simple ``networkx.Graph``.
    ValueError
        If *tree* is empty, is not a tree, or does not contain *root*.
    """
    if not isinstance(tree, nx.Graph):
        raise TypeError(f"tree must be a networkx.Graph, got {type(tree).__name__!r}")
    if tree.is_directed():
        raise TypeError("tree must be undirected; got a directed graph")
    if tree.is_multigraph():
        raise TypeError("tree must be a simple graph; multigraphs are not supported")
    if tree.number_of_nodes() == 0:
        raise ValueError("tree is empty; at least one node is required")
    if not nx.is_tree(tree):
        raise ValueError("input graph is not a tree (a tree must be connected and acyclic)")
    if root not in tree:
        raise ValueError(f"root {root!r} is not a node of the tree")


# ---------------------------------------------------------------------------
# Root finding
# ---------------------------------------------------------------------------
def gravity_root_finder(tree: nx.Graph):
    """Return the centroid of *tree* — a node whose every subtree has <= n/2 nodes.

    Optional helper: pass its result as ``root`` for the most balanced drawing.
    """
    g = tree

    def bfs_parents(r):
        parent = {r: r}
        queue = [r]
        while queue:
            u = queue.pop(0)
            for v in g.neighbors(u):
                if v not in parent:
                    parent[v] = u
                    queue.append(v)
        return parent

    def count_subtrees(u, parent, child_count):
        total = 1
        for v in g.neighbors(u):
            if parent[v] == u:
                count_subtrees(v, parent, child_count)
                total += child_count[v]
        child_count[u] = total

    n = g.number_of_nodes()
    for r in list(g.nodes()):
        degree = g.degree(r)
        parent = bfs_parents(r)
        child_count = {}
        count_subtrees(r, parent, child_count)
        balanced = sum(1 for v in g.neighbors(r) if child_count[v] <= n / 2)
        if balanced == degree:
            return r
    return None  # unreachable for a valid tree


# ---------------------------------------------------------------------------
# Main algorithm
# ---------------------------------------------------------------------------
def oik_symv_1q(
    g: nx.Graph,
    root: object = "A",
    roots_origin: Coord = (0, 0),
    display: bool = False,
    save: bool = False,
    w_labels: bool = False,
    n_size: float = 25,
    n_color: str = "black",
    l_color: str = "black",
    filename: str = "oik_symv_1Q.png",
) -> list:
    """Compute the one-quadrant balanced monotone drawing of tree *g*.

    Parameters
    ----------
    g : networkx.Graph
        An undirected tree.
    root : hashable, default ``"A"``
        The node to root the drawing at. Pass ``gravity_root_finder(g)`` for the
        most balanced result.
    roots_origin : (int, int)
        Grid position placed at the root (useful when composing sub-drawings).
    display, save : bool
        Whether to show the drawing and/or save it to *filename*.
    w_labels, n_size, n_color, l_color :
        Passed through to ``networkx.draw_networkx``.
    filename : str
        Output path used when ``save=True``.

    Returns
    -------
    list
        ``[graph_area, maxx, maxy, position]`` where ``position`` maps each node
        to its ``(x, y)`` grid coordinate and ``graph_area = maxx * maxy``.

    Raises
    ------
    TypeError, ValueError
        If *g* is not a valid tree or *root* is not one of its nodes.
    """
    _validate_tree(g, root)
    if save and not filename:
        raise ValueError("filename must be a non-empty path when save=True")

    G = g.copy()

    def bfs_parents(r):
        parent = {r: r}
        queue = [r]
        while queue:
            u = queue.pop(0)
            for v in G.neighbors(u):
                if v not in parent:
                    parent[v] = u
                    queue.append(v)
        return parent

    def count_subtrees(u, parent, child_count):
        total = 1
        for v in G.neighbors(u):
            if parent[v] == u:
                count_subtrees(v, parent, child_count)
                total += child_count[v]
        child_count[u] = total

    def assign_angles(u, parent, child_count, angle):
        (th1, th2) = angle[u]
        Tu = child_count[u]
        a2 = th1
        for v in G.neighbors(u):
            if parent[v] == u:
                Tv = child_count[v]
                a1 = a2
                a2 = a1 + (th2 - th1) * Tv / (Tu - 1)
                angle[v] = (a1, a2)
                assign_angles(v, parent, child_count, angle)

    def assign_coords(u, parent, angle, position):
        p = parent[u]
        (xp, yp) = position[p]
        (th1, th2) = angle[u]
        th = th2 - th1

        # Geometry preserved verbatim from the original: the first `if` is
        # superseded by the `if/elif/else` block below (which always assigns),
        # kept as-is so the drawing is byte-identical to the reference version.
        if th1 < math.atan(0.5) and math.atan(0.5) < th2 and th > math.atan(0.5):
            (x, y) = (2, 1)
        if th1 < math.pi / 4 and th2 > math.pi / 4:
            (x, y) = (1, 1)
        elif th1 < math.atan(2) and math.atan(2) < th2 and th > math.atan(0.5):
            (x, y) = (1, 2)
        else:
            d = math.ceil(1 / (th2 - th1))
            if th2 <= math.pi / 4:
                (x, y) = (d, math.floor(math.tan(th1) * d + 1))
            elif th2 > math.pi / 2:
                (x, y) = (1, d)
            else:
                (x, y) = (math.floor(math.tan(math.pi / 2 - th2) * d + 1), d)

        if p == u:
            position[u] = roots_origin
        else:
            position[u] = (x + xp, y + yp)

        for v in G.neighbors(u):
            if parent[v] == u:
                assign_coords(v, parent, angle, position)

    # --- run the algorithm --------------------------------------------------
    r = root
    parent = bfs_parents(r)
    n = len(parent)
    child_count: dict[object, int] = {}
    count_subtrees(r, parent, child_count)
    position: dict[object, Coord] = {r: roots_origin}
    angle: dict[object, tuple[float, float]] = {r: (0, math.pi / 2)}
    assign_angles(r, parent, child_count, angle)
    assign_coords(r, parent, angle, position)

    xs = [p[0] for p in position.values()]
    ys = [p[1] for p in position.values()]
    maxx = max(xs) + 1
    maxy = max(ys) + 1
    graph_area = maxx * maxy

    if display or save:
        _draw(
            G,
            position,
            roots_origin,
            maxx,
            maxy,
            n,
            w_labels,
            n_size,
            n_color,
            l_color,
            display,
            save,
            filename,
        )

    return [graph_area, maxx, maxy, position]


def _draw(
    G,
    position,
    roots_origin,
    maxx,
    maxy,
    n,
    w_labels,
    n_size,
    n_color,
    l_color,
    display,
    save,
    filename,
):
    """Render the drawing, always closing the figure afterwards."""
    fig = plt.figure(figsize=(7, 7))
    try:
        nx.draw_networkx(
            G,
            pos=position,
            with_labels=w_labels,
            node_size=n_size,
            node_color=n_color,
            font_color=l_color,
            width=1.2,
        )
        plt.grid(color="gray", alpha=0.5)
        ax = plt.gca()
        ax.set_ylim([-1, maxy + 1])
        ax.set_xlim([-1, maxx + 1])
        plt.title(f"Oikonomou-Symvonis 1-Quadrant\n\nGrid Size: {maxx} x {maxy} ({n} nodes)")
        if roots_origin == (0, 0):
            plt.xticks(np.arange(0, maxx + 1, 1))
            plt.yticks(np.arange(0, maxy + 1, 1))
        ax.set_aspect(1)
        ax.set_xticklabels([])
        ax.set_yticklabels([])
        if save:
            try:
                plt.savefig(filename, bbox_inches="tight")
            except OSError as exc:
                raise OSError(f"could not save drawing to {filename!r}: {exc}") from exc
        if display:
            plt.show()
    finally:
        plt.close(fig)
