"""Oikonomou & Symvonis fine-tuned balanced monotone tree-drawing algorithm.

A "balanced" monotone drawing: the tree is rooted at its centroid (the
gravity/balance node) and each node is given an angular wedge proportional to
the size of its subtree. Edge directions are then chosen as primitive integer
vectors that fall inside each wedge, keeping the drawing compact and monotone.

Unlike the Angelini and He & He algorithms, this one relies on simple geometry
rather than number theory (Stern-Brocot / Farey), which makes it the easiest of
the family to follow.

Reference:
    Anargyros Oikonomou, Antonios Symvonis,
    "Simple Compact Monotone Tree Drawings".
    https://arxiv.org/pdf/1708.09653.pdf
"""

from __future__ import annotations

import math

import matplotlib.pyplot as plt
import networkx as nx

Coord = tuple[int, int]


# ---------------------------------------------------------------------------
# Validation
# ---------------------------------------------------------------------------
def _validate_tree(tree: object) -> None:
    """Validate that *tree* is a non-empty undirected simple ``networkx`` tree.

    Raises
    ------
    TypeError
        If *tree* is not an undirected simple ``networkx.Graph``.
    ValueError
        If *tree* is empty or is not a tree (connected and acyclic).
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


# ---------------------------------------------------------------------------
# Root finding
# ---------------------------------------------------------------------------
def gravity_root_finder(tree: nx.Graph):
    """Return the centroid of *tree* — a node whose every subtree has <= n/2 nodes.

    Such a node always exists for a tree and yields the most balanced rooting.
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
# Geometry: map an angular wedge (th1, th2) to a primitive integer vector
# ---------------------------------------------------------------------------
def get_xy(th1: float, th2: float) -> Coord:
    """Return a primitive integer vector whose direction lies within (th1, th2).

    Assumes ``0 <= th1 < th2 <= pi/2`` (the caller mirrors the other quadrant).
    """
    th = th2 - th1
    if th > math.pi / 4:
        (x, y) = (1, 1)

    elif th <= math.pi / 4 and th > math.atan(0.5):
        if th1 >= math.pi / 4:
            (x, y) = (1, 2)
        elif th1 < math.pi / 4 and th1 >= math.atan(0.5):
            (x, y) = (1, 1)
        elif th1 < math.atan(0.5):
            (x, y) = (2, 1)

    if th <= math.atan(0.5):
        d = math.ceil(1 / (th2 - th1))
        if th2 <= math.pi / 4 and th2 > th1 and th1 >= 0:
            (x, y) = (d, math.floor(math.tan(th1) * d + 1))
        elif th2 > math.pi / 4 and th1 < math.pi / 4 and th2 > th1:
            (x, y) = (1, 1)
        elif th2 > th1 and th1 >= math.pi / 4:
            (x, y) = (math.floor(math.tan(math.pi / 2 - th2) * d + 1), d)

    return (x, y)


# ---------------------------------------------------------------------------
# Main algorithm
# ---------------------------------------------------------------------------
# The leftmost wedge ends at pi; float accumulation leaves it slightly short of
# math.pi, so a generous tolerance is used to recognise "the wedge that reaches
# pi" and point it straight left. (Preserved from the original implementation.)
_NEAR_PI = 3.14


def new_balanced_tree(
    g: nx.Graph,
    display: bool = False,
    save: bool = False,
    roots_path: bool = False,  # reserved / unused; kept for backward compatibility
    filename: str = "tree3.png",
) -> list:
    """Compute the balanced monotone drawing of tree *g* and return its metrics.

    Parameters
    ----------
    g : networkx.Graph
        An undirected tree. Its root is found automatically (the centroid).
    display, save : bool
        Whether to show the drawing and/or save it to *filename*.
    roots_path : bool
        Reserved for a root-path variant; currently unused.
    filename : str
        Output path used when ``save=True``.

    Returns
    -------
    list
        ``[graph_area, maxx, maxy, position]`` where ``position`` maps each node
        to its ``(x, y)`` grid coordinate, ``maxx``/``maxy`` are the horizontal
        and vertical spans, and ``graph_area = maxx * maxy``.

    Raises
    ------
    TypeError, ValueError
        If *g* is not a valid non-empty undirected tree (see ``_validate_tree``).
    """
    _validate_tree(g)
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

        if th1 < math.pi / 2 and th2 > math.pi / 2 and th2 < math.pi:
            (x, y) = (0, 1)
        elif th2 <= math.pi / 2 and th2 > th1 and th2 < math.pi:
            (x, y) = get_xy(th1, th2)
        elif th1 >= math.pi / 2 and th2 > th1 and th2 < math.pi:
            (x1, y1) = get_xy(math.pi - th2, math.pi - th1)
            (x, y) = (-x1, y1)
        if th2 >= _NEAR_PI:
            (x, y) = (-1, 0)

        if p == u:
            position[u] = (0, 0)
        else:
            position[u] = (x + xp, y + yp)

        for v in G.neighbors(u):
            if parent[v] == u:
                assign_coords(v, parent, angle, position)

    # --- run the algorithm --------------------------------------------------
    r = gravity_root_finder(g)
    parent = bfs_parents(r)
    n = len(parent)
    child_count: dict[object, int] = {}
    count_subtrees(r, parent, child_count)
    position: dict[object, Coord] = {r: (0, 0)}
    angle: dict[object, tuple[float, float]] = {r: (0, math.pi)}
    assign_angles(r, parent, child_count, angle)
    assign_coords(r, parent, angle, position)

    xs = [p[0] for p in position.values()]
    ys = [p[1] for p in position.values()]
    maxx = abs(max(xs) - min(xs))
    maxy = abs(max(ys) - min(ys))
    grid_area = n**2
    graph_area = maxx * maxy

    if display or save:
        _draw(G, position, g, graph_area, grid_area, maxx, maxy, n, display, save, filename)

    return [graph_area, maxx, maxy, position]


def _draw(G, position, g, graph_area, grid_area, maxx, maxy, n, display, save, filename):
    """Render the drawing, always closing the figure afterwards."""
    fig = plt.figure(figsize=(7, 7))
    try:
        nx.draw_networkx(
            G, pos=position, with_labels=True, node_size=30, arrows=False, node_color="red"
        )
        plt.grid(color="gray")
        plt.suptitle(
            "Oikonomou-Symvonis Fine Tuned Balanced Algorithm",
            fontsize=16,
            fontweight="bold",
        )
        ax = plt.gca()
        cover = round(graph_area / grid_area, 4) * 100 if grid_area else 0
        xlabel = (
            f"\nNo nodes: {g.number_of_nodes()}, Graph area: {graph_area}, "
            f"Grid area cover (nxn): {cover} %"
            f"\nMax x: {maxx} = {round(maxx / n, 2)}*n,"
            f"\nMax y: {maxy} = {round(maxy / n, 2)}*n\n"
        )
        plt.xlabel(xlabel, horizontalalignment="center")
        ax.set_aspect(1)
        plt.xticks([p[0] for p in position.values()])
        plt.yticks([p[1] for p in position.values()])
        plt.tight_layout(pad=3.5)
        if save:
            try:
                plt.savefig(filename)
            except OSError as exc:
                raise OSError(f"could not save drawing to {filename!r}: {exc}") from exc
        if display:
            plt.show()
    finally:
        plt.close(fig)
