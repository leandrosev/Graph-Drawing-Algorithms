"""Tree generation and conversion utilities.

Every rooted, ordered tree corresponds to a balanced string of parentheses:
starting from the root and reading counter-clockwise, ``'('`` descends to a new
child and ``')'`` returns to the parent. This module generates such strings
(exhaustively or at random) and converts them into ``networkx`` graphs suitable
for the monotone-drawing algorithms in :mod:`Angelini` and :mod:`He_He`.

Node labels follow the scheme ``A, B, ..., Z, AA, AB, ...`` with ``'A'`` as the
root, which is the ordering the drawing algorithms expect.

Input conventions differ between the two exhaustive generators (kept for
backward compatibility):

* :func:`gen_paren_lazy` takes the number of parenthesis *pairs* (non-root nodes).
* :func:`gen_paren_fast` takes the total number of *nodes*.

So ``gen_paren_lazy(3)`` and ``gen_paren_fast(4)`` both enumerate the 4-node trees.

@author: Evangelidakis Leandros
@School of Applied Mathematics and Physical Sciences
"""

from __future__ import annotations

import random
from itertools import permutations
from string import ascii_uppercase

import matplotlib.pyplot as plt
import networkx as nx

# Available node labels: A..Z then AA..ZZ (702 total). 'A' is always the root.
_LABELS = list(ascii_uppercase) + [c1 + c2 for c1 in ascii_uppercase for c2 in ascii_uppercase]
_MAX_NODES = len(_LABELS)


def is_valid(parenth: str) -> bool:
    """Return ``True`` if *parenth* is a balanced string of parentheses.

    The empty string is considered valid (it represents the root-only tree).
    Any character other than ``'('`` or ``')'`` makes the string invalid.

    Raises
    ------
    TypeError
        If *parenth* is not a string.
    """
    if not isinstance(parenth, str):
        raise TypeError(f"expected a string, got {type(parenth).__name__!r}")
    if parenth == "":
        return True
    if parenth[0] == ")":
        return False

    left = right = 0
    for char in parenth:
        if char == "(":
            left += 1
        elif char == ")":
            right += 1
        else:
            return False  # not a parenthesis character
        if right > left:
            return False
    return right == left


def gen_paren_lazy(n: int) -> list[str]:
    """Enumerate every balanced string with *n* pairs of parentheses (brute force).

    *n* is the number of pairs, i.e. the number of non-root nodes; the result
    describes all rooted trees with ``n + 1`` nodes.

    Complexity: ``O((2n)!)`` permutations filtered by :func:`is_valid` -- use
    :func:`gen_paren_fast` for anything but tiny ``n``.

    Raises
    ------
    ValueError
        If *n* is not a non-negative integer.
    """
    if not isinstance(n, int) or n < 0:
        raise ValueError(f"n must be a non-negative integer, got {n!r}")
    perms = {"".join(p) for p in permutations("(" * n + ")" * n)}
    return [s for s in perms if is_valid(s)]


def gen_paren_fast(n: int) -> list[str]:
    """Enumerate every balanced string describing an *n*-node rooted tree.

    *n* is the total number of nodes (including the root). Builds each string
    incrementally, only ever extending valid prefixes.

    Complexity: ``O(4^m / m^1.5)`` for ``m = n - 1`` pairs (the Catalan number).

    Raises
    ------
    ValueError
        If *n* is not a positive integer.
    """
    if not isinstance(n, int) or n < 1:
        raise ValueError(f"n must be a positive integer, got {n!r}")

    pairs = n - 1
    if pairs == 0:
        return [""]  # single-node (root-only) tree

    res: list[str] = []

    def process(string: str = "", opened: int = 0, closed: int = 0) -> None:
        if len(string) == 2 * pairs:
            res.append(string)
            return
        if opened < pairs:
            process(string + "(", opened + 1, closed)
        if closed < opened:
            process(string + ")", opened, closed + 1)

    process()
    return res


def random_tree(n: int) -> str:
    """Return the parenthesis string of a random rooted tree with *n* nodes.

    Uses rejection sampling: shuffle ``n - 1`` pairs of parentheses until a
    balanced arrangement is found (expected ``O(n^1.5)``).

    Raises
    ------
    ValueError
        If *n* is not a positive integer.
    """
    if not isinstance(n, int) or n < 1:
        raise ValueError(f"n must be a positive integer, got {n!r}")

    chars = ["(", ")"] * (n - 1)
    while True:
        random.shuffle(chars)
        parenth = "".join(chars)
        if is_valid(parenth):
            return parenth


def paren_to_nxgraph(parenthesis: str) -> nx.Graph:
    """Convert a balanced parenthesis string into a rooted ``networkx.Graph``.

    Nodes are labelled ``A`` (root), ``B``, ``C``, ... in the order they are
    created, giving the ordered rooted tree the drawing algorithms expect.

    Raises
    ------
    TypeError
        If *parenthesis* is not a string.
    ValueError
        If *parenthesis* is not balanced, or the tree would exceed the
        702-node label scheme.
    """
    if not isinstance(parenthesis, str):
        raise TypeError(f"expected a string, got {type(parenthesis).__name__!r}")
    if not is_valid(parenthesis):
        raise ValueError(f"not a valid balanced-parenthesis string: {parenthesis!r}")

    labels = _LABELS.copy()
    g = nx.Graph()
    g.add_node("A")
    parent = labels.pop(0)  # 'A'
    child = "A"
    index = 0
    nodes = ["A"]
    opened = closed = 0
    parent_of = {"A": "A"}

    for paren in parenthesis:
        if paren == "(":
            if not labels:
                raise ValueError(
                    f"tree too large: the label scheme supports at most {_MAX_NODES} nodes"
                )
            index += 1
            child = labels.pop(0)
            g.add_node(child)
            nodes.append(child)
            g.add_edge(parent, child)
            parent_of[child] = parent
            parent = nodes[index]
            opened += 1
        else:
            parent = parent_of[child]
            child = parent_of[child]
            closed += 1
        if opened == closed:
            opened = closed = 0
            parent = "A"

    return g


def draw(g: nx.Graph) -> None:
    """Display *g* with networkx's spectral layout (a quick, non-monotone view).

    Raises
    ------
    TypeError
        If *g* is not a ``networkx.Graph``.
    """
    if not isinstance(g, nx.Graph):
        raise TypeError(f"g must be a networkx.Graph, got {type(g).__name__!r}")
    fig = plt.figure(figsize=(7, 7))
    try:
        nx.draw_spectral(g, node_size=20)
        plt.title("NetworkX drawing functions")
        plt.show()
    finally:
        plt.close(fig)
