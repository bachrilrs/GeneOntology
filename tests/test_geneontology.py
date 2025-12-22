#!/usr/bin/env python3
# AUTHOR: Bachri Laroussi
# DATE: 2025-12-19
# M1 BBS - S7 Graphs

"""
Mini test-suite for src.geneontology on a tiny GO graph.

Graph orientation (as in your loader):
  child GO term  -->  parent GO term   via 'is a' / 'part of'
  GeneProduct    -->  GO term          via 'annotation'
"""

from src.geneontology import *


def make_mini_go():
    """Build a tiny GO DAG with 3 terms and 2 gene products."""
    go = create_graph(directed=True, weighted=False)
    go["alt_id"] = {}  # keep compatibility with real GO graphs

    # GO terms (biological_process)
    add_node(go, "GO:1", {"type": "GOTerm", "namespace": "biological_process"})  # root (most general)
    add_node(go, "GO:2", {"type": "GOTerm", "namespace": "biological_process"})
    add_node(go, "GO:3", {"type": "GOTerm", "namespace": "biological_process"})  # leaf (most specific)

    # Hierarchy: child -> parent
    add_edge(go, "GO:2", "GO:1", {"relationship": "is a"})
    add_edge(go, "GO:3", "GO:2", {"relationship": "is a"})

    # Gene products
    add_node(go, "P1", {"type": "GeneProduct", "id": "P1"})
    add_node(go, "P2", {"type": "GeneProduct", "id": "P2"})

    # Annotations: GP -> GO term
    add_edge(go, "P1", "GO:2", {"relationship": "annotation", "evidence-codes": ["EXP"]})
    add_edge(go, "P2", "GO:3", {"relationship": "annotation", "evidence-codes": ["IDA"]})

    return go


# Test merging results from different branches
def make_branching_go():
    """
    GO:1 (BP root)
      ├─ GO:2 (BP)
      │    └─ GO:3 (BP)
      └─ GO:4 (MF)

    Annotations:
      P1 -> GO:3 (EXP)
      P2 -> GO:4 (IDA)
    """
    go = create_graph(directed=True, weighted=False)
    go["alt_id"] = {}

    # GO terms (2 namespaces)
    add_node(go, "GO:1", {"type": "GOTerm", "namespace": "biological_process"})
    add_node(go, "GO:2", {"type": "GOTerm", "namespace": "biological_process"})
    add_node(go, "GO:3", {"type": "GOTerm", "namespace": "biological_process"})
    add_node(go, "GO:4", {"type": "GOTerm", "namespace": "molecular_function"})

    # Hierarchy: child -> parent
    add_edge(go, "GO:2", "GO:1", {"relationship": "is a"})
    add_edge(go, "GO:3", "GO:2", {"relationship": "is a"})
    add_edge(go, "GO:4", "GO:1", {"relationship": "is a"})

    # Gene products
    add_node(go, "P1", {"type": "GeneProduct", "id": "P1"})
    add_node(go, "P2", {"type": "GeneProduct", "id": "P2"})

    # Annotations
    add_edge(go, "P1", "GO:3", {"relationship": "annotation", "evidence-codes": ["EXP"]})
    add_edge(go, "P2", "GO:4", {"relationship": "annotation", "evidence-codes": ["IDA"]})

    return go

def test_geneproducts_union_recursive():
    go = make_branching_go()

    # GO:1 doit récupérer P1 (via GO:3) et P2 (via GO:4)
    res = set(GeneProducts(go, "GO:1", recursive=True))
    assert res == {"P1", "P2"}

    # GO:2 doit récupérer P1 (via GO:3) mais pas P2
    res2 = set(GeneProducts(go, "GO:2", recursive=True))
    assert res2 == {"P1"}

    # GO:4 est une feuille -> seulement P2
    assert GeneProducts(go, "GO:4", recursive=True) == ["P2"]



def test_node_types():
    go = make_mini_go()
    assert is_goterm(go, "GO:1")
    assert is_goterm(go, "GO:2")
    assert is_goterm(go, "GO:3")
    assert is_geneproduct(go, "P1")
    assert is_geneproduct(go, "P2")


def test_counts():
    go = make_mini_go()
    assert count_goterm(go) == 3
    assert count_geneproducts(go) == 2
    assert count_annotations(go) == 2

def test_count_goterm_by_namespace():
    go = make_branching_go()
    d = count_goterm(go, by_namespace=True)

    assert d["biological_process"] == 3   # GO:1, GO:2, GO:3
    assert d["molecular_function"] == 1   # GO:4
    assert d["cellular_component"] == 0

def test_evidence_code_distribution():
    go = make_branching_go()
    dist = evidence_code_distribution(go)

    assert dist["EXP"] == 1
    assert dist["IDA"] == 1

def tests_summary_statistics():
    go = make_branching_go()
    stats = summary(go)
    
    assert stats["total_goterms"] == 4
    assert stats["total_gp"] == 2
    assert stats["annotations"] == 2
    assert stats["avg_ann_per_GP"] == 1.0

def test_goterms_direct_and_recursive():
    go = make_mini_go()

    assert GOTerms(go, "P1", recursive=False) == ["GO:2"]
    assert GOTerms(go, "P1", recursive=True) == ["GO:1", "GO:2"]

    assert GOTerms(go, "P2", recursive=False) == ["GO:3"]
    assert GOTerms(go, "P2", recursive=True) == ["GO:1", "GO:2", "GO:3"]



def test_geneproducts_direct_and_recursive():
    go = make_mini_go()

    # Testing Precomputing indexes for efficiency
    idx = go_to_gps_index(go)
    children = build_children_index(go)
    assert 	GeneProducts(go, "GO:1", True) == GeneProducts(go, "GO:1", True, idx, children)

    # Also test without precomputed indexes
    # Direct annotations
    assert GeneProducts(go, "GO:1", recursive=False) == []
    assert GeneProducts(go, "GO:2", recursive=False) == ["P1"]
    assert GeneProducts(go, "GO:3", recursive=False) == ["P2"]

    # Recursive: include descendants (more specific terms)
    assert set(GeneProducts(go, "GO:1", recursive=True)) == {"P1", "P2"}
    assert set(GeneProducts(go, "GO:2", recursive=True)) == {"P1", "P2"}
    assert GeneProducts(go, "GO:3", recursive=True) == ["P2"]

    
    

def test_max_depth_bp():
    go = make_mini_go()

    # Subgraph restricted to GO terms of biological_process
    dag = induced_goterm_subgraph(go, namespace="biological_process")

    # Deepest path should end at the leaf GO:3 (depth 2)
    max_d, deepest, path = max_depth_go(
        dag,
        namespace=None,          # already filtered by induced_goterm_subgraph
        return_path=True,
        reverse_path=True        # path from most specific -> most general
    )

    assert max_d == 2
    assert deepest == "GO:3"
    assert path == ["GO:3", "GO:2", "GO:1"]


def run_tests():
    test_node_types()
    test_counts()
    test_goterms_direct_and_recursive()
    test_geneproducts_direct_and_recursive()
    test_max_depth_bp()
    test_geneproducts_union_recursive()
    test_count_goterm_by_namespace()
    test_evidence_code_distribution()
    tests_summary_statistics()

    print("All geneontology tests passed.")


if __name__ == "__main__":
    run_tests()