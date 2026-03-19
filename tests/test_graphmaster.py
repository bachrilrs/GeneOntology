#!/usr/bin/env python3
# AUTHOR: Bachri Laroussi
# DATE: 2025-12-19
# M1 BBS - S7 Graphs

import math
import src.graphmaster as gm


def test_create_graph():
    g = gm.create_graph()
    assert g["nodes"] == {}
    assert g["edges"] == {}
    assert g["directed"] is True


def test_add_node_and_edge():
    g = gm.create_graph()
    gm.add_node(g, "A", {"type": "X"})
    assert "A" in g["nodes"]
    assert g["nodes"]["A"]["type"] == "X"

    gm.add_edge(g, "A", "B")
    assert gm.edge_exists(g, "A", "B")
    assert gm.node_exists(g, "A")
    assert gm.node_exists(g, "B")


def test_nb_nodes_nb_edges_directed():
    g = gm.create_graph()
    gm.add_edge(g, "A", "B")
    gm.add_edge(g, "A", "C")

    assert gm.nb_nodes(g) == 3
    assert gm.nb_edges(g) == 2


def test_nb_edges_neighbors_undirected():
    g = gm.create_graph(directed=False)
    gm.add_edge(g, "A", "B")
    gm.add_edge(g, "B", "C")

    assert gm.nb_edges(g) == 2
    assert gm.neighbors(g, "A") == ["B"]
    assert gm.neighbors(g, "B") == ["A", "C"]


def test_neighbors_directed():
    g = gm.create_graph()
    gm.add_edge(g, "A", "B")
    gm.add_edge(g, "A", "C")
    gm.add_edge(g, "C", "D")

    assert gm.neighbors(g, "A") == ["B", "C"]
    assert gm.neighbors(g, "C") == ["D"]
    assert gm.neighbors(g, "B") == []


def test_read_delim_degrees():
    g = gm.read_delim("data/dressing.tsv", column_separator="\t")

    assert gm.edge_exists(g, "socks", "shoes")
    assert gm.out_degree(g, "underwear") == 2
    assert gm.out_degree(g, "tie") == 1
    assert gm.in_degree(g, "jacket") == 2
    assert gm.in_degree(g, "shoes") == 3


def test_bfs_basic():
    g = gm.create_graph()
    gm.add_edge(g, "A", "B")
    gm.add_edge(g, "A", "C")
    gm.add_edge(g, "B", "D")
    gm.add_edge(g, "C", "D")

    state, dist, pred = gm.bfs(g, "A")

    assert dist["A"] == 0
    assert dist["B"] == 1
    assert dist["C"] == 1
    assert dist["D"] == 2
    assert pred["D"] in ("B", "C")


def test_is_acyclic_true():
    g = gm.create_graph(directed=True)
    gm.add_edge(g, "A", "B")
    gm.add_edge(g, "A", "C")
    gm.add_edge(g, "B", "D")
    gm.add_edge(g, "C", "D")

    parcours = gm.dfs(g)

    assert gm.is_acyclic(g) is True
    assert "BACK EDGE" not in parcours["classification"].values()


def test_is_acyclic_false():
    g = gm.create_graph()
    gm.add_edge(g, "X", "Y")
    gm.add_edge(g, "Y", "Z")
    gm.add_edge(g, "Z", "X")

    assert gm.is_acyclic(g) is False


def test_topological_sort():
    g = gm.read_delim("data/dressing.tsv", column_separator="\t")
    ordre = gm.topological_sort(g)
    pos = {n: i for i, n in enumerate(ordre)}

    assert pos["underwear"] < pos["pants"]
    assert pos["pants"] < pos["belt"]
    assert pos["belt"] < pos["jacket"]
    assert pos["underwear"] < pos["shoes"]
    assert pos["socks"] < pos["shoes"]


def test_path():
    g = gm.create_graph()
    gm.add_edge(g, "A", "B")
    gm.add_edge(g, "C", "D")

    res = gm.bfs(g, "A")

    assert gm.path(res, "A", "B") is not None
    assert gm.path(res, "A", "D") is None


def test_bfs_unreachable_nodes():
    g = gm.create_graph()
    gm.add_edge(g, "A", "B")
    gm.add_edge(g, "C", "D")

    state, dist, pred = gm.bfs(g, "A")

    assert dist["A"] == 0
    assert dist["B"] == 1
    assert dist["C"] == -math.inf
    assert dist["D"] == -math.inf
    assert pred["C"] is None
    assert pred["D"] is None


def test_select_nodes():
    g = gm.create_graph()
    gm.add_node(g, "GO:1", {"namespace": "biological_process"})
    gm.add_node(g, "GO:2", {"namespace": "cellular_component"})

    assert gm.select_nodes(g, "namespace", "biological_process") == ["GO:1"]


def test_filter_edges():
    g = gm.create_graph()
    gm.add_edge(g, "A", "B", {"relationship": "is a"})
    gm.add_edge(g, "B", "C", {"relationship": "part of"})

    assert gm.filter_edges(g, "relationship", "is a") == [("A", "B")]