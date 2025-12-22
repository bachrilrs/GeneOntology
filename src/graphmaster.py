#!/usr/bin/env python3
# AUTHOR: Bachri Laroussi
# DATE: 2025-11-26
# M1 BBS - S7 Graphs
from pprint import pprint
import math


def create_graph(directed = True, weighted = False, weight_attribute = None):
    """Create and return an empty graph structure.

    Parameters
    ----------
    directed : bool, optional
        If True, edges are directed (default True).
    weighted : bool, optional
        If True, edges may carry a numeric weight.
    weight_attribute : str or None, optional
        Name of the edge attribute used to store weights.

    Returns
    -------
    dict
        A graph represented as a dictionary with 'nodes' and 'edges' mappings.
    """
    g = { 'nodes': {}, 'edges': {}, 'directed': directed, 'weighted': weighted, 'weight_attribute': weight_attribute }
    return g
        

def add_node(g, node_id, attributes = None): 
    """Add a node to the graph and return its attribute dictionary.

    Parameters
    ----------
    g : dict
        Graph object.
    node_id :
        Identifier of the node (string, int, ...).
    attributes : dict, optional
        Initial attributes for the node.

    Returns
    -------
    dict
        The attribute dictionary associated with this node.
    """
    if node_id not in g['nodes']: # ensure node does not already exist
        if attributes is None: # create empty attributes if not provided
            attributes = {}
        g['nodes'][node_id] = attributes
        g['edges'][node_id] = {} # init outgoing edges
    return g['nodes'][node_id] # return node attributes
        

def add_edge(g, node_id1, node_id2, attributes = None): 
    """Add a (directed) edge from node_id1 to node_id2 and return its attributes.

    If the nodes do not exist yet, they are created automatically.

    Parameters
    ----------
    g : dict
        Graph object.
    node_id1 : hashable
        Source node identifier.
    node_id2 : hashable
        Target node identifier.
    attributes : dict, optional
        Initial attributes for the edge.

    Returns
    -------
    dict
        The attribute dictionary associated with this edge.
    """
    if not node_exists(g , node_id1): add_node(g, node_id1) # ensure node1 exists
    if not node_exists(g , node_id2): add_node(g, node_id2) # ensure node2 exists
   
    # add edge(s) only if they do not exist
    if not edge_exists(g,node_id1,node_id2):
        if attributes is None: # create empty attributes if not provided
            attributes = {}
        g['edges'][node_id1][node_id2] = attributes
        if not g['directed']:
            g['edges'][node_id2][node_id1] = g['edges'][node_id1][node_id2] # share the same attributes as n1->n2
    return g['edges'][node_id1][node_id2] # return edge attributes

def node_exists(g, n):
    """Return True if the given node identifier exists in the graph."""
    return n in g['nodes']
	
def edge_exists(g, n1 , n2):
    """Return True if an edge from n1 to n2 exists in the graph."""
    return node_exists(g,n1) and n2 in g['edges'][n1] 

def nodes(g):
    """Return the list of node identifiers of the graph, sorted."""
    return sorted(g['nodes'].keys())
        
        

def read_delim(filename, column_separator='\t', directed=True, weighted=False, weight_attribute=None): 
    """Read an edge list from a delimited text file and build a graph.

    The file format is:
        node_id1 <sep> node_id2 <sep> att1 <sep> att2 ...

    Parameters
    ----------
    filename : str
        Path to the input file.
    column_separator : str, optional
        Column separator used in the file (default '\\t').
    directed : bool, optional
        If True, build a directed graph.
    weighted : bool, optional
        If True, edges may carry a numeric weight.
    weight_attribute : str or None, optional
        Name of the edge attribute used to store weights.

    Returns
    -------
    dict
        A graph built from the file.
    """
    g = create_graph(directed, weighted, weight_attribute)
    with open(filename) as f: 
        # GET COLUMNS NAMES
        tmp = f.readline().rstrip()
        attNames= tmp.split(column_separator)
        # REMOVES FIRST TWO COLUMNS WHICH CORRESPONDS TO THE LABELS OF THE CONNECTED VERTICES
        attNames.pop(0)  # remove first column name (source node not to be in attribute names)
        attNames.pop(0)  # remove second column (target node ...)
        # PROCESS THE REMAINING LINES
        row = f.readline().rstrip()
        while row:
            vals = row.split(column_separator)
            u = vals.pop(0)
            v = vals.pop(0)
            att = {}
            for i in range(len(attNames)):
                att[ attNames[i] ] = vals[i]
            add_edge(g, u, v, att)
            row = f.readline().rstrip() # NEXT LINE
        return g

def nb_nodes(g):
    """Return the number of nodes in the graph."""
    return len(g['nodes'])

def nb_edges(g):
    """Return the number of edges in the graph.
    For undirected graphs, each edge is counted once.
    """
    total = sum(len(adj) for adj in g['edges'].values())
    return total if g['directed'] else total //2

   
def neighbors(g, node_id):
    """Return the sorted list of neighbors (successors) of a given node."""

    return sorted(list(g['edges'].get(node_id, {}).keys()))


def enqueue(queue,element):
    """Append an element at the end of a queue (FIFO list)."""
    queue.append(element)


def dequeue(queue):
    """Remove and return the first element of a queue (FIFO list).

    Raises
    ------
    IndexError
        If the queue is empty.
    """
    if len(queue)==0:
        raise IndexError("Empty queue")
    return queue.pop(0) # first in first out, so pop(0)
        

def out_degree(g, node_id):
    """Return the out-degree of a node (number of outgoing edges)."""
    return len(g['edges'][node_id])

def in_degree(g,node_id):
    """Return the in-degree of a node (number of incoming edges).
    Visits all nodes to count how many have an edge to node_id."""
    cpt= 0
    for n in nodes(g):
        if node_id in neighbors(g,n):
            cpt +=1
    return cpt

def degree(g,node_id):
    """Return the total degree of a node (in-degree + out-degree)."""
    return (in_degree(g,node_id) + out_degree(g,node_id))
        

def select_nodes(g,attribut,attribut_value):
    """Return a list of node identifiers having a given attribute value.

    Parameters
    ----------
    g : dict
        Graph object.
    attribut : str
        Attribute name to test on nodes.
    attribut_value : any
        Value the attribute must have.

    Returns
    -------
    list of hashable
        Node identifiers matching the condition.
        Important point for the other functions that use Select nodes.
    """
    res = list()
    for n in nodes(g):
        if attribut in g['nodes'][n] and g['nodes'][n][attribut] == attribut_value:
            res.append(n)
    return res 


def filter_edges(g,attribut , attribut_value):
    """Return a list of edges (u, v) whose attribute matches a given value.

    Parameters
    ----------
    g : dict
        Graph object.
    attribut : str
        Edge attribute name.
    attribut_value : any
        Value the attribute must have.

    Returns
    -------
    list of tuple
        List of (u, v) edge pairs.
    """
    res = list()
    for u in nodes(g):
        for v in neighbors(g,u):
            if attribut in g['edges'][u][v] and g['edges'][u][v][attribut] == attribut_value:
                res.append((u,v)) # adding edge in res in form of tuple
    return res


def relationships(g):
    """Print the number of edges for each relationship type in the graph."""
    list = ['is a' , 'part of' , 'annotation']
    res = {}
    for el in list:
        res[el] = len(filter_edges(g,'relationship',el))	
    return res
	

def bfs(g,s):
    """Perform a breadth-first search (BFS) from source node s.

    Parameters
    ----------
    g : dict
        Graph object.
    s : hashable
        Source node identifier.

    Returns
    -------
    tuple
        (state, distance, predecessor) dictionaries:
        - state[u] in {'Unexplored', 'Discovered', 'Finished'}
        - distance[u] is the BFS distance from s (or -inf if unreachable)
        - predecessor[u] is the parent of u in the BFS tree (or None).
    """
    state = {}
    distance = {}
    predecessor = {}

    for node in nodes(g):
        state[node] = 'Unexplored'
        distance[node] = -math.inf # float('-inf') works too
        predecessor[node] = None
        
    state[s]= 'Discovered'   
    distance[s] = 0
    queue = []
    enqueue(queue,s)
    while queue :	
        u = dequeue(queue)
        for v in neighbors(g,u):
            if state[v]=='Unexplored':
                state[v] = 'Discovered'
                distance[v]= distance[u] + 1
                predecessor[v] = u
                enqueue(queue,v)
                
        state[u]='Finished'
        
    return state , distance , predecessor


def path(bfs_res, source, target):
    """Reconstruct the shortest path from source to target using BFS results.

    Parameters
    ----------
    bfs_res : tuple
        Result of bfs(g, source): (state, distance, predecessor).
    source : hashable
        Source node identifier.
    target : hashable
        Target node identifier.

    Returns
    -------
    str or None
        A string describing the path and its distance, or None if no path exists.
    """
    state, distance, predecessor = bfs_res
        
    if source == target:
        return [source]
    route = []
    cur = target
    while cur is not None and cur != source:
        route.append(cur)
        cur = predecessor[cur]

    if cur is None:
        return None

    route.append(source)
    route.reverse()
    return f"{route} \n La distance entre {source} et {target} est de {distance[target]}."

def dfs_visit(g,u,parcours):
    """Recursive DFS visit from node u, updating the traversal structure.

    Parameters
    ----------
    g : dict
        Graph object.
    u : hashable
        Starting node for this DFS visit.
    parcours : dict
        DFS traversal structure containing state, times and edge classifications.
    """
    parcours['state'][u] = 'Discovered' 
    parcours['time'] += 1
    parcours['discovered'][u]= parcours['time']

    for v in neighbors(g,u):
        if parcours['state'][v] == 'Unexplored':
            parcours['predecessor'][v] = u
            parcours['classification'][(u,v)] = 'TREE EDGE'
            dfs_visit(g,v,parcours)
        elif parcours['state'][v] == 'Discovered':
            parcours['classification'][(u,v)] = 'BACK EDGE'
        elif parcours['discovered'][u] > parcours['discovered'][v]:
            parcours['classification'][(u,v) ] = "CROSS EDGE"
        else:
            parcours['classification'][(u,v)] = "FORWARD EDGE"
    parcours['state'][u] = 'Finished'
    parcours['time'] += 1
    parcours['finished'][u] = parcours['time']
		
def dfs(g):
    """Perform a depth-first search (DFS) on the whole graph.

    Returns
    -------
    dict
        Traversal information with keys:
        'time', 'state', 'discovered', 'finished',
        'classification', 'predecessor'.
    """
    parcours = {'time' : 0,
        'state' : {},
        'discovered' : {},
        'finished' : {},
        'classification' : {},
        'predecessor' : {},
    }
   
    for u in nodes(g):
        parcours['state'][u] = 'Unexplored'
        parcours['predecessor'][u] = None


    for v in nodes(g):
        if parcours['state'][v] == 'Unexplored':
            dfs_visit(g,v,parcours)
    return parcours


		
def is_acyclic(g):
    """Check whether the graph is acyclic using DFS edge classification.

    Returns
    -------
    str
        A textual message indicating if the graph is a DAG or contains a cycle.
    """

    parcours = dfs(g)
    for el in parcours['classification'].values():
        if el == 'BACK EDGE':
            return False
    return True

	

def topological_sort(g):
    """Compute a topological ordering of the nodes of a DAG.

    The graph must be acyclic; otherwise a ValueError is raised.

    Returns
    -------
    list of hashable
        List of node identifiers ordered topologically.
    """
    # Check graph is acyclic.
    if not is_acyclic(g):
        raise ValueError("impossible to sort: DAG.")
    if nb_nodes(g) == 0:
        return False
        
    parcours = dfs(g)

    # Gets nodes and their 'finished' time.
    nodes = list(parcours['finished'].keys())
    fins = parcours['finished']

    # Sort nodes by finished time 
    # for i in range(len(nodes)):
    #    for j in range(i+1, len(nodes)):
    #       if fins[nodes[i]] < fins[nodes[j]]:
    #          som = nodes[i]
    #          nodes[i] = nodes[j]
    #          nodes[j] = som

        
    sorted_nodes = list(fins.keys())
    sorted_nodes.sort(key=lambda n: fins[n], reverse=True)

        
    return sorted_nodes

def transpose_graph(g):
    """Return the transpose of a directed graph.

    All nodes are copied, and for each edge u→v an edge v→u is created.

    Parameters
    ----------
    g : dict
        Original graph.

    Returns
    -------
    dict
        New graph representing the transpose of g.
    """
    gt = create_graph(
        directed = g['directed'],
        weighted = g['weighted'],
        weight_attribute = g['weight_attribute']
    )
        
    # Add nodes
    for n in g['nodes']:
        add_node(gt, n, g['nodes'][n].copy())

    # Reverse edges
    for u in g['edges']:
        for v in g['edges'][u]:
            add_edge(gt, v, u, g['edges'][u][v].copy())

    return gt

def induced_subgraph(g, nodes):
    """Build the subgraph induced by a given set of nodes.

    Parameters
    ----------
    g : dict
        Original graph.
    nodes : iterable
        Iterable of node identifiers to keep.

    Returns
    -------
    dict
        New graph containing only these nodes and edges between them.
    """
    sg = create_graph(g['directed'], g['weighted'], weight_attribute=g['weight_attribute'])

    # add nodes (with their attributs)
    for t in nodes:
        add_node(sg, t, g['nodes'][t].copy())

    # add edges if they're in nodes
    node_set = set(nodes)
    for t in nodes:
        for p in neighbors(g, t):  
            if p in node_set:
                add_edge(sg, t, p, g['edges'][t][p].copy())
    return sg
		

if __name__ == "__main__":
	print("# Graphmaster module")

