#!/usr/bin/env python3
# AUTHOR: Bachri Laroussi
# DATE: 2025-12-21
# M1 BBS - S7 Graphs

"""
Docstring for test_OBO_GOA
Mini test-suite for src.geneontology on real GO graph and GOA annotations.
Graph orientation (as in your loader):
  child GO term  -->  parent GO term   via 'is a' / 'part of'
  GeneProduct    -->  GO term          via 'annotation'
""" 
from src.geneontology import *
from pprint import pprint   
import os

if not os.path.exists('data/go-basic.obo'):
    print("Please download go-basic.obo from https://current.geneontology.org/ontology/go-basic.obo")
    exit(1)
if not os.path.exists('data/25.H_sapiens.goa'):
    print("Please download 25.H_sapiens.goa from https://ftp.ebi.ac.uk/pub/databases/GO/goa/proteomes/25.H_sapiens.goa.")
    exit(1)

print('Gene Ontology functions test.\nNote that we used the versions that are precised in Projet_GeneOntology.pdf\n')

g = load_OBO('data/go-basic.obo')
load_GOA(g, 'data/25.H_sapiens.goa')
print('GO graph loaded from OBO and GOA files.\n')

assert ('P04637' in GeneProducts(g, 'GO:0003677')) == True # molecular_function
assert ('P04637' in GeneProducts(g, 'GO:0005634',recursive=True)) == True # recursive cellular_component

assert ('GO:0003677' in GOTerms(g, 'P04637')) == True  # list of GO terms for TP53
assert ('GO:0005634' in GOTerms(g, 'P04637',recursive=True)) == True  # list of GO terms for TP53


assert max_depth_go(g,'biological_process') == 17  # max depth of the GO graph
assert max_depth_go(g,'molecular_function') == 12  
assert max_depth_go(g,'cellular_component') == 14  


gt = transpose_graph(g) # in both cases, we should obtain the same results

assert max_depth_go(gt,'biological_process') == 17   # max depth of the transposed GO graph
assert max_depth_go(gt,'molecular_function') == 12  
assert max_depth_go(gt,'cellular_component') == 14  

print('All tests passed successfully.')