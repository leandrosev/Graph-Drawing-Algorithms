"""Angelini et al. monotone tree-drawing algorithms (BFS-based and DFS-based).

A monotone drawing of a rooted tree assigns each edge a slope drawn from a
carefully chosen set of primitive vectors so that, for every pair of nodes,
the path between them is monotone with respect to some direction.

This module implements the two variants of Angelini et al.:

* ``"bfs"`` -- slopes come from the Stern-Brocot / ``fusc`` sequence
  (OEIS A002487), giving a balanced grid.
* ``"dfs"`` -- slopes are ``k/1`` for ``k = 1 .. n-1``, giving a tall, narrow
  grid.

Reference:
    Angelini, Colasante, Di Battista, Frati, Patrignani,
    "Monotone Drawings of Graphs", Journal of Graph Algorithms and
    Applications, 2012.
    https://www.emis.de/journals/JGAA/accepted/2012/Angelini+2012.16.1.pdf

Assumptions
-----------
The input must be an undirected tree (connected, acyclic) whose node ordering
starts at ``root_node`` and lists parents before their children -- exactly the
ordering produced by :func:`Tree_Operations.paren_to_nxgraph`. These
assumptions are validated and, when violated, raise informative exceptions
rather than producing a silently wrong drawing.

@author: Evangelidakis Leandros
@School of Applied Mathematics and Physical Sciences
"""

from __future__ import annotations

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np

# Stern-Brocot ratios can involve division by (temporarily) zero denominators;
# the results are never used, so silence the warnings rather than pollute output.
np.seterr(divide="ignore", invalid="ignore")

VALID_ALGORITHMS: tuple[str, ...] = ("bfs", "dfs")

# A node's coordinate is a (y, x) pair; positions handed to matplotlib are (x, y).
Coord = tuple[float, float]


# ---------------------------------------------------------------------------
# Validation
# ---------------------------------------------------------------------------
def _validate_tree(tree: object, root_node: object) -> None:
    """Validate that *tree* is a rooted undirected tree containing *root_node*.

    Raises
    ------
    TypeError
        If *tree* is not an undirected simple ``networkx.Graph``.
    ValueError
        If *tree* is empty, is not a tree, or does not contain *root_node*.
    """
    if not isinstance(tree, nx.Graph):
        raise TypeError(f"tree must be a networkx.Graph, got {type(tree).__name__!r}")
    if tree.is_directed():
        raise TypeError("tree must be undirected; got a directed graph")
    if tree.is_multigraph():
        raise TypeError("tree must be a simple graph; multigraphs are not supported")
    if tree.number_of_nodes() == 0:
        raise ValueError("tree is empty; at least the root node is required")
    if not nx.is_tree(tree):
        raise ValueError("input graph is not a tree (a tree must be connected and acyclic)")
    if root_node not in tree:
        raise ValueError(f"root_node {root_node!r} is not a node of the tree")


# ---------------------------------------------------------------------------
# Slope generation
# ---------------------------------------------------------------------------
def _slope_map(n: int, algo: str) -> dict[float, Coord]:
    """Return ``ratio -> (y, x)`` for the ``n - 1`` edge slopes, sorted by ratio.

    For ``"bfs"`` the numerators/denominators are consecutive terms of the
    Stern-Brocot (``fusc``) sequence; for ``"dfs"`` the slopes are ``k/1``.
    An empty mapping is returned for a single-node tree (``n == 1``).
    """
    if algo == "bfs":
        v = np.zeros(2 * n)
        v[0] = 1
        for i in range(n):
            v[2 * i] = v[i]
            v[2 * i + 1] = v[i] + v[i + 1]
        x = v[1:n]
        y = v[2 : n + 1]
    else:  # "dfs"; algo is already validated by the caller
        y = np.arange(1, n, dtype=float)
        x = np.ones(n - 1)

    # ratio -> (y, x); the ratios of a Stern-Brocot run (and of k/1) are unique,
    # so no two edges collapse onto the same dictionary key.
    ratio_to_coord = dict(zip(list(np.divide(y, x)), list(zip(y, x))))
    return dict(sorted(ratio_to_coord.items()))


# ---------------------------------------------------------------------------
# Drawing
# ---------------------------------------------------------------------------
def _draw(
    g: nx.Graph,
    pos: dict[object, Coord],
    algo: str,
    maxx: int,
    maxy: int,
    n: int,
    w_labels: bool,
    n_size: float,
    n_color: str,
    l_color: str,
    display: bool,
    save: bool,
    filename: str,
) -> None:
    """Render the drawing, optionally saving and/or displaying it.

    The matplotlib figure is always closed, even if saving or showing fails,
    so repeated calls do not leak figures.
    """
    titles = {
        "bfs": "Angelini et al. BFS-based",
        "dfs": "Angelini et al. DFS-based",
    }
    fig = plt.figure(figsize=(7, 7))
    try:
        nx.draw_networkx(
            g,
            pos=pos,
            with_labels=w_labels,
            node_size=n_size,
            arrows=False,
            node_color=n_color,
            font_color=l_color,
        )
        plt.grid(color="gray")
        plt.title(f"{titles[algo]}\n\nGrid Size: {maxx} x {maxy} ({n} nodes)")

        ax = plt.gca()
        ax.set_ylim([-1, maxy + 1])
        ax.set_xlim([-1, maxx + 1])
        plt.xticks(np.arange(0, maxx + 1, 1))
        plt.yticks(np.arange(0, maxy + 1, 1))
        ax.set_yticklabels([])
        ax.set_xticklabels([])

        if save:
            try:
                plt.savefig(filename, bbox_inches="tight")
            except OSError as exc:
                raise OSError(f"could not save drawing to {filename!r}: {exc}") from exc
        if display:
            plt.show()
    finally:
        plt.close(fig)


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------
def getGridArea(
    tree: nx.Graph,
    root_node: object = "A",
    display: bool = False,
    save: bool = False,
    algo: str = "bfs",
    view: object = None,  # kept for backward compatibility; unused
    filename: str = "tree1.png",
    w_labels: bool = False,
    n_size: float = 25,
    n_color: str = "black",
    l_color: str = "black",
) -> list:
    """Compute a monotone drawing of *tree* and return its grid metrics.

    Parameters
    ----------
    tree : networkx.Graph
        An undirected rooted tree, ordered from ``root_node`` with parents
        before children (as produced by ``Tree_Operations.paren_to_nxgraph``).
    root_node : hashable, default ``"A"``
        The root; must be the first node in the tree's ordering.
    display, save : bool
        Whether to show the drawing in a window and/or save it to ``filename``.
    algo : {"bfs", "dfs"}
        Which slope-assignment variant to use (case-insensitive).
    view : unused
        Accepted for backward compatibility; ignored.
    filename : str
        Output path used when ``save=True``.
    w_labels, n_size, n_color, l_color :
        Passed through to ``networkx.draw_networkx``.

    Returns
    -------
    list
        ``[graph_area, maxx, maxy, gridpos]`` where ``gridpos`` maps each node
        to its integer ``(x, y)`` grid position and ``graph_area = maxx * maxy``.

    Raises
    ------
    TypeError
        If *tree* is not an undirected simple ``networkx.Graph``.
    ValueError
        If *tree* is empty/not a tree, ``root_node`` is missing or not the
        first node, ``algo`` is unknown, or ``save=True`` with an empty
        ``filename``.
    RuntimeError
        If the tree does not satisfy the rooted/ordered layout assumptions and
        coordinates cannot be resolved.
    """
    # --- validate arguments -------------------------------------------------
    algo = str(algo).lower()
    if algo not in VALID_ALGORITHMS:
        raise ValueError(f"algo must be one of {VALID_ALGORITHMS}, got {algo!r}")
    _validate_tree(tree, root_node)
    if save and not filename:
        raise ValueError("filename must be a non-empty path when save=True")

    g = tree.copy()
    n = g.number_of_nodes()

    # Adjacency list as node -> [children...] in the tree's (CCW) node ordering.
    graph: dict[object, list] = {}
    for line in nx.generate_adjlist(g):
        parts = line.split(" ")
        graph[parts[0]] = parts[1:]

    # The layout relies on the ordering starting at the root.
    first_node = next(iter(graph))
    if first_node != root_node:
        raise ValueError(
            f"root_node {root_node!r} must be the first node in the tree's "
            f"ordering, but the ordering starts at {first_node!r}. Build the "
            "tree so the root is inserted first (e.g. via "
            "Tree_Operations.paren_to_nxgraph)."
        )

    try:
        gridpos = _compute_positions(g, graph, n, algo, root_node)
    except (KeyError, IndexError) as exc:
        raise RuntimeError(
            "failed to resolve node coordinates; ensure the tree is rooted at "
            f"{root_node!r} with nodes ordered from the root "
            "(as produced by Tree_Operations.paren_to_nxgraph)"
        ) from exc

    xs = [pos[0] for pos in gridpos.values()]
    ys = [pos[1] for pos in gridpos.values()]
    maxx = int(max(xs)) + 1
    maxy = int(max(ys)) + 1
    graph_area = maxx * maxy

    if display or save:
        _draw(
            g,
            gridpos,
            algo,
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

    return [graph_area, maxx, maxy, gridpos]


def _compute_positions(
    g: nx.Graph,
    graph: dict[object, list],
    n: int,
    algo: str,
    root_node: object,
) -> dict[object, tuple[int, int]]:
    """Assign each node an integer ``(x, y)`` grid position.

    This is the core of the Angelini et al. algorithm: decompose the tree into
    per-node subtrees, assign each subtree a contiguous run of sorted slopes,
    turn slopes into per-node offsets, and accumulate offsets from the root.
    """
    # Map each node to the list of its child subtrees (each a list of nodes).
    subtree_map: dict[object, list] = {}
    for k in list(g.nodes()):
        children = graph[k]
        if not children:
            subtree_map[k] = [k]
            continue
        newg = g.copy()
        newg.remove_node(k)
        subtree_map[k] = []
        for child in children:
            component = list(nx.bfs_tree(newg, child).nodes())
            subtree_map[k].append(component)
            newg.remove_node(child)

    # Sorted slope ratios and the mapping ratio -> (y, x) coordinate.
    slope_map = _slope_map(n, algo)
    sorted_ratios = list(slope_map.keys())

    # For each node, split its slope run among its subtrees and hand each child
    # the last (steepest) slope of its subtree as the child's own offset.
    #   subtrees_seq: tuple(subtree nodes) -> slope run
    #   coordmap:     node -> (y, x) offset relative to its parent
    #   T_u:          node -> slope run available to its own subtree
    subtrees_seq: dict[tuple, list] = {}
    coordmap: dict[object, Coord] = {root_node: (0, 0)}
    T_u: dict[object, list] = {root_node: list(sorted_ratios)}
    for node in graph:
        children = graph[node]
        if not children:
            continue
        lengths = [0]
        for subtree in subtree_map[node]:
            lengths.append(len(subtree))
        for j, subtree in enumerate(subtree_map[node]):
            start = 1 + sum(lengths[: j + 1])
            end = sum(lengths[: j + 2])
            subtrees_seq[tuple(subtree)] = T_u[node][start - 1 : end]
        for child in children:
            for subtree in subtree_map[node]:
                if child in subtree:
                    seq = subtrees_seq[tuple(subtree)]
                    coordmap[child] = slope_map[seq[-1]]
                    T_u[child] = seq
                    break

    # Accumulate offsets from the root to get absolute (y, x) positions.
    grid_map: dict[object, Coord] = {root_node: (0, 0)}
    for parent in graph:
        py, px = grid_map[parent]
        for child in graph[parent]:
            if child not in grid_map:
                cy, cx = coordmap[child]
                grid_map[child] = (cy + py, cx + px)

    # Expose positions to matplotlib as integer (x, y) tuples.
    return {node: (int(grid_map[node][1]), int(grid_map[node][0])) for node in g.nodes()}
