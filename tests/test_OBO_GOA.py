#!/usr/bin/env python3
# AUTHOR: Bachri Laroussi
# DATE: 2025-12-21
# M1 BBS - S7 Graphs

"""
Mini test-suite for src.geneontology on real GO graph and GOA annotations.
Graph orientation (as in your loader):
  child GO term  -->  parent GO term   via 'is a' / 'part of'
  GeneProduct    -->  GO term          via 'annotation'
"""

from pathlib import Path
import pytest
from src.geneontology import *

OBO_FILE = Path("data/go-basic.obo")
GOA_FILE = Path("data/25.H_sapiens.goa")

pytestmark = [
    pytest.mark.slow,
    pytest.mark.skipif(
        not OBO_FILE.exists() or not GOA_FILE.exists(),
        reason="Fichiers GOA/OBO absents"
    )
]


@pytest.fixture(scope="module")
def go_graph():
    g = load_OBO(str(OBO_FILE))
    load_GOA(g, str(GOA_FILE))
    return g


def test_geneproducts_direct(go_graph):
    assert "P04637" in GeneProducts(go_graph, "GO:0003677")


def test_geneproducts_recursive(go_graph):
    assert "P04637" in GeneProducts(go_graph, "GO:0005634", recursive=True)


def test_goterms_direct(go_graph):
    assert "GO:0003677" in GOTerms(go_graph, "P04637")


def test_goterms_recursive(go_graph):
    assert "GO:0005634" in GOTerms(go_graph, "P04637", recursive=True)


def test_max_depth_go_original(go_graph):
    assert max_depth_go(go_graph, "biological_process") == 17
    assert max_depth_go(go_graph, "molecular_function") == 12
    assert max_depth_go(go_graph, "cellular_component") == 14


def test_max_depth_go_transposed(go_graph):
    gt = transpose_graph(go_graph)
    assert max_depth_go(gt, "biological_process") == 17
    assert max_depth_go(gt, "molecular_function") == 12
    assert max_depth_go(gt, "cellular_component") == 14