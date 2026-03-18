#!/usr/bin/env python3
# AUTHOR: Bachri Laroussi
# DATE: 2025-12-19
# M1 BBS - S7 Graphs
from pprint import pprint
import src.graphmaster as gm
import math

 
print('Graphmaster functions test.')
print("We used data/dressing.tsv for some tests.")
print('...')

g = gm.create_graph()
assert g['nodes'] == {}
assert g['edges'] == {}
assert g['directed'] == True


# Test add_node / add edge

gm.add_node(g, 'A', {'type': 'X'})
assert 'A' in g['nodes']
assert g['nodes']['A']['type'] == 'X'

gm.add_edge(g , 'A' , 'B')
assert gm.edge_exists(g, 'A', 'B')
assert gm.node_exists(g, 'A')
assert gm.node_exists(g, 'B')

# Test nb_nodes / nb_edges
g = gm.create_graph()
gm.add_edge(g, 'A', 'B')
gm.add_edge(g, 'A', 'C')

assert gm.nb_nodes(g) == 3
assert gm.nb_edges(g) == 2  # directed graph

# Test for a non directed graph
g = gm.create_graph(directed=False)
gm.add_edge(g, 'A', 'B')
gm.add_edge(g, 'B', 'C')


assert gm.nb_edges(g) == 2   #
assert gm.neighbors(g, 'A') == ['B']
assert gm.neighbors(g, 'B') == ['A', 'C']

# Test neighbors
g = gm.create_graph()
gm.add_edge(g, 'A', 'B')
gm.add_edge(g, 'A', 'C')
gm.add_edge(g, 'C', 'D')

assert gm.neighbors(g, 'A') == ['B', 'C']
assert gm.neighbors(g, 'C') == ['D']
assert gm.neighbors(g, 'B') == []

g = gm.read_delim('data/dressing.tsv' , column_separator='\t')

assert gm.edge_exists(g, 'socks', 'shoes')

assert gm.out_degree(g , 'underwear') == 2
assert gm.out_degree(g , 'tie') == 1
assert gm.in_degree(g , 'jacket') == 2
assert gm.in_degree(g , 'shoes') == 3

g = gm.create_graph()
gm.add_edge(g, 'A', 'B')
gm.add_edge(g, 'A', 'C')
gm.add_edge(g, 'B', 'D')
gm.add_edge(g, 'C', 'D')

# Test BFS
state, dist, pred = gm.bfs(g, 'A')

assert dist['A'] == 0
assert dist['B'] == 1
assert dist['C'] == 1
assert dist['D'] == 2
assert pred['D'] in ('B','C')

# Test is_acyclic function 
g = gm.create_graph(directed=True)
gm.add_edge(g, 'A', 'B')
gm.add_edge(g, 'A', 'C')
gm.add_edge(g, 'B', 'D')
gm.add_edge(g, 'C', 'D')

parcours = gm.dfs(g)

assert gm.is_acyclic(g) == True
assert 'BACK EDGE' not in parcours['classification'].values()

g = gm.create_graph()
gm.add_edge(g, 'X', 'Y')
gm.add_edge(g, 'Y', 'Z')
gm.add_edge(g, 'Z', 'X')

assert gm.is_acyclic(g) == False

g = gm.read_delim('data/dressing.tsv' , column_separator='\t')
res = gm.bfs(g , 'underwear')


# Test de topological_sort()
ordre = gm.topological_sort(g)
pos = {n: i for i, n in enumerate(ordre)}

assert pos['underwear'] < pos['pants']
assert pos['pants'] < pos['belt']
assert pos['belt'] < pos['jacket']
assert pos['underwear'] < pos['shoes']
assert pos['socks'] < pos['shoes']



# Test path() 
g = gm.create_graph()
gm.add_edge(g, 'A', 'B')
gm.add_edge(g, 'C', 'D')  # isolated component

res = gm.bfs(g, 'A')

#  A → B exists
assert gm.path(res, 'A', 'B') is not None

#  A → D doesn't exists
assert gm.path(res, 'A', 'D') is None

# Test BFS
g = gm.create_graph()
gm.add_edge(g, 'A', 'B')
gm.add_edge(g, 'C', 'D')

state, dist, pred = gm.bfs(g, 'A')

assert dist['A'] == 0
assert dist['B'] == 1

assert dist['C'] == -math.inf
assert dist['D'] == -math.inf
assert pred['C'] is None
assert pred['D'] is None

# Test for select_nodes and filter_edges
g = gm.create_graph()
gm.add_node(g, 'GO:1', {'namespace':'biological_process'})
gm.add_node(g, 'GO:2', {'namespace':'cellular_component'})
assert gm.select_nodes(g, 'namespace', 'biological_process') == ['GO:1']

g = gm.create_graph()
gm.add_edge(g, 'A','B', {'relationship':'is a'})
gm.add_edge(g, 'B','C', {'relationship':'part of'})
assert gm.filter_edges(g, 'relationship', 'is a') == [('A','B')]

print("Everything's good!")
