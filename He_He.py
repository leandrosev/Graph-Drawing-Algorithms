"""He & He optimal monotone tree-drawing algorithm.

Implementation of the path-based monotone drawing of Dayu He and Xin He,
"Optimal Monotone Drawings of Trees" (https://arxiv.org/pdf/1604.03921v1.pdf).

The idea: decompose the tree into a *length-decreasing path decomposition*
(LDPD), group those paths into ``c``-partition levels, build a matching set of
primitive slope vectors from Farey sequences, assign each path a vector, and
accumulate the vectors from the root to obtain grid coordinates.

The public entry point is :func:`path_draw_algorithm`. The remaining functions
(``farey_seq``, ``ldpd``, ``c_partition``, ``construct_prim_vectors`` ...) are
the building blocks and are kept individually importable.

Assumptions
-----------
The input must be an undirected tree (connected, acyclic) whose node ordering
starts at ``root`` and lists parents before their children -- exactly the
ordering produced by :func:`Tree_Operations.paren_to_nxgraph`. Violations are
validated and raise informative exceptions instead of a silently wrong drawing.

@author: Evangelidakis Leandros
@School of Applied Mathematics and Physical Sciences
"""

from __future__ import annotations

import math
from string import ascii_uppercase

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np

Vector = tuple[int, int]


# ---------------------------------------------------------------------------
# Validation
# ---------------------------------------------------------------------------
def _validate_tree(tree: object, root: object) -> None:
    """Validate that *tree* is a rooted undirected tree containing *root*.

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
        raise ValueError("tree is empty; at least the root node is required")
    if not nx.is_tree(tree):
        raise ValueError("input graph is not a tree (a tree must be connected and acyclic)")
    if root not in tree:
        raise ValueError(f"root {root!r} is not a node of the tree")


# ---------------------------------------------------------------------------
# Farey sequences (primitive slope vectors)
# ---------------------------------------------------------------------------
def invert(tup: Vector) -> Vector:
    """Swap the two components of a 2-tuple: ``(a, b) -> (b, a)``."""
    return (tup[1], tup[0])


def farey_seq(n: int, descending: bool = False) -> list[Vector]:
    """Return the ``n``-th Farey sequence closed under coordinate inversion.

    Each fraction ``a/b`` is represented as the tuple ``(a, b)``; the result is
    sorted by slope ``b / a`` and contains, for every ``(a, b)``, its inverse
    ``(b, a)`` as well.
    """
    result = []
    a, b, c, d = 0, 1, 1, n
    if descending:
        a, c = 1, n - 1
    result.append((a, b))

    while (c <= n and not descending) or (a > 0 and descending):
        k = int((n + b) / d)
        a, b, c, d = c, d, (k * c - a), (k * d - b)
        result.append((a, b))

    result.pop(0)
    for i in range(len(result)):
        if invert(result[i]) not in result:
            result.append(invert(result[i]))

    return sorted(result, key=lambda x: x[1] / float(x[0]))


def farey_seq2(n: int, descending: bool = False) -> list[Vector]:
    """Like :func:`farey_seq`, but bracketed by the axis vectors.

    The vertical vector ``(1, 0)`` is prepended and the horizontal vector
    ``(0, 1)`` appended, giving a sequence spanning the full first quadrant.
    """
    result = []
    a, b, c, d = 0, 1, 1, n
    if descending:
        a, c = 1, n - 1
    result.append((a, b))
    while (c <= n and not descending) or (a > 0 and descending):
        k = int((n + b) / d)
        a, b, c, d = c, d, (k * c - a), (k * d - b)
        result.append((a, b))

    result.pop(0)
    for i in range(len(result)):
        if invert(result[i]) not in result:
            result.append(invert(result[i]))

    return [(1, 0)] + sorted(result, key=lambda x: x[1] / float(x[0])) + [(0, 1)]


# ---------------------------------------------------------------------------
# Path decompositions
# ---------------------------------------------------------------------------
def _adjacency(tree: nx.Graph) -> dict[object, list]:
    """Return ``node -> [children...]`` in the tree's node ordering."""
    graph: dict[object, list] = {}
    for line in nx.generate_adjlist(tree):
        parts = line.split(" ")
        graph[parts[0]] = parts[1:]
    return graph


def path_decomposition(tree: nx.Graph, root: object = "A") -> list[list]:
    """Decompose *tree* into leaf-to-junction paths (first leaf reaches root)."""
    graph = _adjacency(tree)
    leafs = [v for v in graph if graph[v] == []]

    b = list(nx.shortest_path(tree, leafs[0], root))
    seen = set(b)
    bset = [b]
    for v in leafs[1:]:
        path_to_root = list(nx.shortest_path(tree, v, root))
        for u in path_to_root:
            if u in seen:
                bset.append(list(nx.shortest_path(tree, v, u)))
                seen.update(path_to_root)
                break

    for path in bset:
        del path[-1]

    return bset


def ldpd(tree: nx.Graph, root: object = "A") -> list[list]:
    """Length-Decreasing Path Decomposition of *tree*.

    Returns a list of paths, ordered by non-increasing length, that partition
    the tree's edges: the longest leaf-to-root path first, then each remaining
    leaf's path up to the first already-covered node.
    """
    graph = _adjacency(tree)
    leafs = [v for v in graph if graph[v] == []]

    leafs_paths = [list(nx.shortest_path(tree, u, root)) for u in leafs]

    maxlen_path = max(leafs_paths, key=len)
    covered = set(maxlen_path)
    paths = [maxlen_path]
    used_leafs = [maxlen_path[0]]

    while len(leafs) != len(used_leafs):
        candidates = []
        for v in leafs:
            if v not in used_leafs:
                path_to_root = list(nx.shortest_path(tree, v, root))
                for u in path_to_root:
                    if u in covered:
                        candidates.append(list(nx.shortest_path(tree, v, u)))
                        break
        longest = max(candidates, key=len)
        paths.append(longest)
        covered.update(longest)
        # The leaf that owns this path is its first node. (The original code
        # used ``longest[0][0]``, which is only the first *character* of the
        # label and breaks for multi-character labels, i.e. trees with > 26
        # nodes; ``longest[0]`` is the full leaf label.)
        used_leafs.append(longest[0])

    return paths


def c_partition(tree: nx.Graph, ldpd: list[list], c: int) -> list[list[list]]:
    """Partition an LDPD into ``ceil(log_c(n))`` levels by path length."""
    n = tree.number_of_nodes()
    K = math.ceil(math.log(n, c))
    D = []

    level0 = [b for b in ldpd if (n - 1) / c <= len(b) - 1 <= (n - 1)]
    D.append(level0)

    for j in range(2, K + 1):
        level = [b for b in ldpd if (n - 1) / (c**j) <= len(b) - 1 < (n - 1) / (c ** (j - 1))]
        D.append(level)

    return D


def construct_prim_vectors(f: int, d: int, n: int) -> list:
    """Build the primitive slope-vector set and its per-level breakdown.

    Returns ``[union_of_R, All_Rs]`` where ``union_of_R`` is the full sorted
    vector list and ``All_Rs`` holds the vectors introduced at each level.
    """
    c = f + 1
    K = math.ceil(math.log(n, c))

    def find_elements_between(source, start, end, howmany):
        begin = source.index(start)
        end = source.index(end)
        between = source[begin + 1 : end]
        return between[0:howmany]

    P = farey_seq(d)

    P1 = P.index((1, 1))
    S1 = P[0:f]
    S2 = P[P1 + 1 : P1 + f + 1]
    R1 = list(S1) + [(1, 1)] + list(S2)

    Pd2 = farey_seq2(d**2)
    R2 = []
    newR1 = [(1, 0)] + R1 + [(0, 1)]

    for first, second in zip(newR1, newR1[1:]):
        R2 += find_elements_between(Pd2, first, second, f)

    All_Rs = [R1, R2]
    union_of_R = [*R1, *R2]
    union_of_R = sorted(union_of_R, key=lambda x: x[1] / float(x[0]))

    for j in range(3, K + 1):
        if (1, 0) not in union_of_R and (0, 1) not in union_of_R:
            union_of_R.insert(0, (1, 0))
            union_of_R.append((0, 1))

        R = []
        jsource = farey_seq2(d**j)
        for first, second in zip(union_of_R, union_of_R[1:]):
            R += find_elements_between(jsource, first, second, f)

        union_of_R.remove((0, 1))
        union_of_R.remove((1, 0))
        union_of_R += R
        union_of_R = sorted(union_of_R, key=lambda x: x[1] / float(x[0]))
        All_Rs.append(R)

    return [union_of_R, All_Rs]


def get_leafs(tree: dict[object, list], start: object = None) -> list:
    """Return the leaves (childless nodes) of an adjacency-list ``dict``."""
    return [node for node in tree if tree[node] == []]


# ---------------------------------------------------------------------------
# Drawing
# ---------------------------------------------------------------------------
def _draw_hehe(
    tree: nx.Graph,
    pos: dict[object, Vector],
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
    """Render the drawing, always closing the figure afterwards."""
    fig = plt.figure(figsize=(7, 5))
    try:
        nx.draw_networkx(
            tree,
            pos=pos,
            with_labels=w_labels,
            node_size=n_size,
            arrows=False,
            node_color=n_color,
            font_color=l_color,
        )
        plt.grid(color="gray")
        ax = plt.gca()
        ax.set_ylim([-1, maxy + 1])
        ax.set_xlim([-1, maxx + 1])
        ax.set_xticklabels([])
        ax.set_yticklabels([])
        plt.title(f"He & He Optimal Algorithm\n\nGrid Size: {maxx + 1} x {maxy + 1} ({n} nodes)")
        ax.set_aspect(1)
        plt.xticks(np.arange(0, maxx + 1, 1))
        plt.yticks(np.arange(0, maxy + 1, 1))

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
def path_draw_algorithm(
    g: nx.Graph,
    root: object = "A",
    display: bool = False,
    save: bool = False,
    filename: str = "He_He.png",
    f: int = 3,
    d: int = 3,
    w_labels: bool = False,
    n_size: float = 25,
    n_color: str = "black",
    l_color: str = "black",
) -> list:
    """Compute the He & He monotone drawing of *g* and return its grid metrics.

    Parameters
    ----------
    g : networkx.Graph
        An undirected rooted tree, ordered from *root* with parents before
        children (as produced by ``Tree_Operations.paren_to_nxgraph``).
    root : hashable, default ``"A"``
        The root; must be the first node in the tree's ordering.
    display, save : bool
        Whether to show and/or save the drawing to *filename*.
    filename : str
        Output path used when ``save=True``.
    f, d : int
        Algorithm parameters controlling the primitive-vector construction
        (``f`` vectors per gap, Farey order ``d``). Both must be >= 1.
    w_labels, n_size, n_color, l_color :
        Passed through to ``networkx.draw_networkx``.

    Returns
    -------
    list
        ``[graph_area, maxx, maxy, gridpos]`` where ``gridpos`` maps each node
        to its ``(x, y)`` grid position and ``graph_area = maxx * maxy``.
        ``maxx``/``maxy`` are the maximum coordinates (the drawn grid spans
        ``maxx + 1`` by ``maxy + 1`` points).

    Raises
    ------
    TypeError
        If *g* is not an undirected simple ``networkx.Graph``.
    ValueError
        If *g* is empty/not a tree, *root* is missing or not the first node,
        ``f``/``d`` are not positive integers, or ``save=True`` with an empty
        ``filename``.
    RuntimeError
        If the tree does not satisfy the rooted/ordered layout assumptions and
        the drawing cannot be constructed.
    """
    # --- validate arguments -------------------------------------------------
    _validate_tree(g, root)
    if not (isinstance(f, int) and f >= 1):
        raise ValueError(f"f must be an integer >= 1, got {f!r}")
    if not (isinstance(d, int) and d >= 1):
        raise ValueError(f"d must be an integer >= 1, got {d!r}")
    if save and not filename:
        raise ValueError("filename must be a non-empty path when save=True")

    tree = g.copy()
    graph = _adjacency(tree)

    first_node = next(iter(graph))
    if first_node != root:
        raise ValueError(
            f"root {root!r} must be the first node in the tree's ordering, but "
            f"the ordering starts at {first_node!r}. Build the tree so the root "
            "is inserted first (e.g. via Tree_Operations.paren_to_nxgraph)."
        )

    try:
        gridpos = _compute_positions(tree, graph, root, f, d)
    except (KeyError, IndexError, ValueError, ZeroDivisionError) as exc:
        raise RuntimeError(
            "failed to construct the He & He drawing; ensure the tree is rooted "
            f"at {root!r} with nodes ordered from the root, and that f/d are "
            "large enough for the tree size"
        ) from exc

    n = tree.number_of_nodes()
    xs = [pos[0] for pos in gridpos.values()]
    ys = [pos[1] for pos in gridpos.values()]
    maxx = int(max(xs))
    maxy = int(max(ys))
    graph_area = maxx * maxy

    if display or save:
        _draw_hehe(
            tree,
            gridpos,
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
    tree: nx.Graph,
    graph: dict[object, list],
    root: object,
    f: int,
    d: int,
) -> dict[object, Vector]:
    """Assign each node an ``(x, y)`` grid position via the He & He method."""
    leafs = get_leafs(graph, root)
    t = len(leafs)
    n = tree.number_of_nodes()

    union_of_R, all_Rs = construct_prim_vectors(f, d, n)

    # Order leaves/paths by the integer rank of their leaf label.
    letters = list(ascii_uppercase) + [c1 + c2 for c1 in ascii_uppercase for c2 in ascii_uppercase]
    for_relabel = dict(zip(letters, range(n + 1)))

    # NOTE: B and D deliberately share the same path-list objects; truncating
    # the paths in B below is observed through D when building b_levels.
    B = ldpd(tree, root)
    B = sorted(B, key=lambda path: for_relabel[path[0]])
    D = c_partition(tree, B, f + 1)

    # Level of each primitive vector (its index in all_Rs).
    R_levels = {}
    for i, level in enumerate(all_Rs):
        for vec in level:
            R_levels[tuple(vec)] = i

    # Drop the shared endpoint of each path (it belongs to the parent path).
    for path in B:
        path.remove(path[-1])

    # Level of each (now truncated) path.
    b_levels = {}
    for i, level in enumerate(D):
        for path in level:
            b_levels[tuple(path)] = i

    # Assign to each path the next unused vector whose level fits the path.
    assigned_edges = {}
    last_index = 0
    for li in range(t):
        bl = B[li]
        bl_level = b_levels[tuple(bl)]
        for vector in union_of_R:
            if R_levels[vector] <= bl_level and union_of_R.index(vector) > last_index:
                assigned_edges[tuple(bl)] = vector
                last_index = union_of_R.index(vector)
                break

    # Every node inherits the vector of the path it belongs to.
    abs_coords = {}
    for path in B:
        for node in path:
            abs_coords[node] = assigned_edges[tuple(path)]

    # Accumulate vectors from the root to obtain absolute coordinates.
    grid_points: dict[object, Vector] = {root: (0, 0)}
    for node in graph:
        for child in graph[node]:
            if child not in grid_points:
                (cx, cy) = abs_coords[child]
                (px, py) = grid_points[node]
                grid_points[child] = (cx + px, cy + py)

    return grid_points
