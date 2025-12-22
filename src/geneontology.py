#!/usr/bin/env python3
# AUTHOR: Bachri Laroussi
# DATE: 2025-12-10
# M1 BBS - S7 Graphs

import re
from collections import Counter
from src.graphmaster import *

file_obo = 'data/go-basic.obo'
def load_OBO(filename='go-basic.obo'):
	"""
	parse the OBO file and returns the graph
	obsolete terms are discarded
	only is_a and part_of relationships are loaded

	Extract of a file to be parsed:
	[Term]
	id: GO:0000028
	name: ribosomal small subunit assembly
	namespace: biological_process
	def: "The aggregation, arrangement and bonding together of constituent RNAs and proteins to form the small ribosomal subunit." [GOC:jl]
	subset: gosubset_prok
	synonym: "30S ribosomal subunit assembly" NARROW [GOC:mah]
	synonym: "40S ribosomal subunit assembly" NARROW [GOC:mah]
	is_a: GO:0022618 ! ribonucleoprotein complex assembly
	relationship: part_of GO:0042255 ! ribosome assembly
	relationship: part_of GO:0042274 ! ribosomal small subunit biogenesis
	"""
	def parseTerm(lines):
		# search for obsolete
		for l in lines:
			if l.startswith('is_obsolete: true'):
				return
		# otherwise create node
		go_id = re_go_id.match(lines.pop(0)).group(1)
		go_attr = add_node(go_graph, go_id) # add node to graph and get the node attribute dict
		go_attr['type'] = 'GOTerm'
		for line in lines:
			if re_go_name.match(line): go_attr['name'] = re_go_name.match(line).group(1)
			elif re_go_namespace.match(line): go_attr['namespace'] = re_go_namespace.match(line).group(1)
			elif re_go_def.match(line): go_attr['def'] = re_go_def.match(line).group(1)
			elif re_go_alt_id.match(line): go_graph['alt_id'][ re_go_alt_id.match(line).group(1) ] = go_id  # alt_id → go_id
			elif re_go_is_a.match(line): 
				parent_id = re_go_is_a.match(line).group(1)
				add_edge(go_graph, go_id, parent_id, { 'relationship': 'is a' })
			elif re_go_part_of.match(line): 
				parent_id = re_go_part_of.match(line).group(1)
				add_edge(go_graph, go_id, parent_id, { 'relationship': 'part of' })
	# method main
	go_graph          = create_graph(directed=True, weighted=False)
	go_graph['alt_id'] = {} # alternate GO ids
	# regexp to parse term lines
	re_go_id          = re.compile(r'^id:\s+(GO:\d+)\s*$')
	re_go_name        = re.compile(r'^name:\s+(.+)\s*$')
	re_go_namespace   = re.compile(r'^namespace:\s+(.+)\s*$')
	re_go_def         = re.compile(r'^def:\s+(.+)\s*$')
	re_go_alt_id      = re.compile(r'^alt_id:\s+(GO:\d+)\s*$')
	re_go_is_a        = re.compile(r'^is_a:\s+(GO:\d+)\s')
	re_go_xref        = re.compile(r'^xref:\s+(\S+)\s*$')
	re_go_part_of      = re.compile(r'^relationship:\s+part_of\s+(GO:\d+)\s')
	# buffer each term lines, then parse lines to create GOTerm node
	with open(filename) as f:
		line = f.readline().rstrip()
		# skip header until first [Term] is reached
		while not line.startswith('[Term]'): 
			line = f.readline().rstrip()
		buff = []  
		line = f.readline()
		stop = False
		while line and not stop:
			line = line.rstrip()
			# new Term
			if line.startswith('[Term]'):
				parseTerm(buff)
				buff=[]
			# last Term
			elif line.startswith('[Typedef]'):
				parseTerm(buff)
				stop=True
			# or append to buffer
			else:
				buff.append(line)
			line = f.readline()
	return go_graph

def load_GOA(go, filename, warnings=True):
	"""
	parse GOA file and add annotated gene products to previsouly loaded graph go

	Extract of a file to be parsed:
	gaf-version: 2.1
	!GO-version: http://purl.obolibrary.org/obo/go/releases/2020-11-28/extensions/go-plus.owl
	UniProtKB       O05154  tagX            GO:0008360      GO_REF:0000043  IEA     UniProtKB-KW:KW-0133    P       Putative glycosyltransferase TagX       tagX|SAOUHSC_00644      protein 93061   20201128        UniProt 
			
	UniProtKB       O05154  tagX            GO:0016740      GO_REF:0000043  IEA     UniProtKB-KW:KW-0808    F       Putative glycosyltransferase TagX       tagX|SAOUHSC_00644      protein 93061   20201128        UniProt 
			
	UniProtKB       O05204  ahpF            GO:0000302      GO_REF:0000002  IEA     InterPro:IPR012081      P       Alkyl hydroperoxide reductase subunit F ahpF|SAOUHSC_00364      protein 93061   20201128        InterPro
			
		0        1       2   3       4             5          6        7      8             9                              10
				id    name        go_id               evidence-codes                     desc                           aliases

	GAF spec: http://geneontology.org/docs/go-annotation-file-gaf-format-2.1/
	Column 	Content 						Required? 	Cardinality 	Example
	1 		DB 								required 	1 				UniProtKB
	2 		DB Object ID 					required 	1 				P12345
	3 		DB Object Symbol 				required 	1 				PHO3
	4 		Qualifier 						optional 	0 or greater 	NOT
	5 		GO ID 							required 	1 				GO:0003993
	6 		DB:Reference (|DB:Reference) 	required 	1 or greater 	PMID:2676709
	7 		Evidence Code 					required 	1 				IMP
	8 		With (or) From 					optional 	0 or greater 	GO:0000346
	9 		Aspect 							required 	1 				F
	10 		DB Object Name 					optional 	0 or 1 			Toll-like receptor 4
	11 		DB Object Synonym (|Synonym) 	optional 	0 or greater 	hToll 	Tollbooth
	12 		DB Object Type 					required 	1 				protein
	13 		Taxon(|taxon) 					required 	1 or 2 			taxon:9606
	14 		Date 							required 	1 				20090118
	15 		Assigned By 					required 	1 				SGD
	16 		Annotation Extension 			optional 	0 or greater 	part_of(CL:0000576)
	17 		Gene Product Form ID 			optional 	0 or 1 			UniProtKB:P12345-2
	"""
	with open(filename) as f: 
		line = f.readline()
		while line:
			if not line.startswith('!'): # skip comments
				cols = line.rstrip().split('\t')
				gp_id = cols[1]
				gt_id = cols[4]
				if gt_id not in go['nodes']: # GOTerm not found search alternate ids
					while gt_id not in go['nodes'] and gt_id in go['alt_id']:
						gt_id = go['alt_id'][gt_id] # replace term by alternate
				if gt_id not in go['nodes']: # failure: warn user
					if warnings:
						print(f'Warning: could not attach a gene product ({gp_id}) to a non existing GO Term ({gt_id})')
				else: # success: GOTerm to attach to was found
					# create node for gene product if not already present
					if gp_id not in go['nodes']:
						gp_attr = add_node(go, gp_id, { 'id': gp_id, 'type': 'GeneProduct'})
					# create or update gene product attributes
					gp_attr = go['nodes'][gp_id]
					gp_attr['name'] = cols[2]
					gp_attr['desc'] = cols[9]
					gp_attr['aliases'] = cols[10].split('|')
					# attach gene product to GOTerm
					gt_attr = go['nodes'][gt_id]
					e_attr = add_edge(go, gp_id, gt_id)
					e_attr['relationship'] = 'annotation'
					if 'evidence-codes' not in e_attr:
						e_attr['evidence-codes'] = []
					e_attr['evidence-codes'].append( cols[6] )
			line = f.readline()

	

def is_goterm(go, node_id):
	"""Return True if the given node is a GOTerm."""
	return go['nodes'][node_id].get('type') == 'GOTerm'

def is_geneproduct(go, node_id):
	"""Return True if the given node is a GeneProduct."""
	return go['nodes'][node_id].get('type') == 'GeneProduct'

def go_parents(go, go_id): # just for clarity
	"""Return the direct GO parents (more general terms) of a GO term.
	child --> parent
	Parents are the neighbors reached by outgoing GO→GO edges
	(e.g. 'is a', 'part of').
	"""
	return neighbors(go,go_id)

def go_to_gps_index(go):
	"""Index: go_id -> set of GeneProducts directly annotated to it.
	Parameters
	----------
	go : dict
		Gene Ontology graph.
		
	Returns
	-------
	dict
		Mapping go_id -> set of GeneProduct identifiers."""
	index = {}
	for gp_id, adj in go["edges"].items():
		if not is_geneproduct(go, gp_id): # if it's not a GeneProduct, skip
			continue
		# adj: go_id -> edge attributes
		for go_id, attr in adj.items():
			if attr.get("relationship") == "annotation" and is_goterm(go, go_id):
				if go_id not in index:
					index[go_id] = set() # initialize set
				index[go_id].add(gp_id) # add GeneProduct to the set
	return index # we get a dict: go_id -> set of GeneProducts

def build_children_index(go):
	"""
	Index: parent_go_id -> set(child_go_id) for GO hierarchy edges (is a / part of).
	This function will be useful for recursive search of GeneProducts.
	"""
	children = {}
	for child_id, adj in go["edges"].items():
		# only GO terms
		if not is_goterm(go, child_id):
			continue

		for parent_id, attr in adj.items():
			if not is_goterm(go, parent_id):
				continue
			rel = attr.get("relationship")
			if rel in ("is a", "part of"):
				if parent_id not in children: # initialize set
					children[parent_id] = set()
				children[parent_id].add(child_id) # add child to parent

	return children

def GeneProducts(go, go_id, recursive=False, go_to_gps=None, children_index=None):
	"""
	Return GeneProducts associated with a given GO term.
	parameters
	----------
	go : dict
		Gene Ontology graph.
	go_id : str
		GO term identifier.
	recursive : bool, optional
		If False, return only GeneProducts directly annotated to go_id.
		If True, include GeneProducts annotated to all descendant GO terms
		reached via 'is a' / 'part of' edges.
	go_to_gps : dict, optional
		Precomputed index: go_id -> set of GeneProducts directly annotated to it.

	children_index : dict, optional
		Precomputed index: parent_go_id -> set(child_go_id) for GO hierarchy edges
		(is a / part of).
	Returns
	-------
	list of str
		Sorted list of GeneProduct identifiers.
		"""
	if not node_exists(go, go_id) or not is_goterm(go, go_id):
		return []

	if go_to_gps is None: # can be precomputed for efficiency
		go_to_gps = go_to_gps_index(go)
	if children_index is None:  # can be precomputed for efficiency
		children_index = build_children_index(go)

	if not recursive:
		return sorted(go_to_gps.get(go_id, set()))	
	
	# BFS/queue on GO descendants using children_index (fast)
	terms = {go_id}
	queue = [go_id]

	while queue:
		t = queue.pop(0)
		for child in children_index.get(t, ()):
			if child not in terms:
				terms.add(child)
				queue.append(child)
	
	# Union of GP sets for all collected terms
	res = set()
	for t in terms:
		res.update(go_to_gps.get(t, set()))
	return sorted(res)

def GOTerms(go, gp_id, recursive=False):
	"""Return GO terms associated with a given GeneProduct.

	Parameters
	----------
	go : dict
		Gene Ontology graph.
	gp_id : str
		GeneProduct identifier.
	recursive : bool, optional
		If False, return only directly annotated terms (successors of the GP).
		If True, include all ancestors (more general GO terms) reached
		via 'is a' / 'part of' edges.

	Returns
	-------
	list of str
		Sorted list of GO term identifiers.
	"""

	if gp_id not in go['nodes']:
		return []
	
	# GOTerms directly associated 
	directs = [v for v in neighbors(go , gp_id) if is_goterm(go,v)]
	
	if not recursive:
		return sorted(set(directs))
	
	
	# 2) recursive=True → we also add ancestors
	result = set(directs)
	attente = list(directs)

	while len(attente)!=0: 
		t = attente.pop(0)
		for parent in go_parents(go,t): # looking for parents
			if parent not in result and is_goterm(go,parent):
				result.add(parent)
				attente.append(parent)
	return sorted(result)

# Depth computation

def goterm_ids(go, namespace=None):
	"""Return list of GO term node ids (optionally filtered by namespace)."""
	res = []
	for nid, attr in go["nodes"].items():
		if attr.get("type") != "GOTerm":
			continue # skip non-GO terms
		if namespace is None or attr.get("namespace") == namespace:
			# print(nid,attr)
			res.append(nid) # keep
	return res

def induced_goterm_subgraph(go, namespace=None):
	"""Subgraph induced by GO terms only (optionally filtered by namespace)."""

	keep = set(goterm_ids(go, namespace=namespace))
	sg = create_graph(directed=True, weighted=False)
	sg["alt_id"] = {}  # keeps compatibility with GO graph structure in case

	for t in keep:
		add_node(sg, t, go["nodes"][t].copy())

	for child in keep: # for each GO term
		for parent in neighbors(go, child): 
			if parent in keep: # only GO terms
				rel = go["edges"][child][parent].get("relationship")
				if rel in ("is a", "part of"): # only hierarchy edges
					add_edge(sg, child, parent, go["edges"][child][parent].copy())
	return sg

def max_depth_go(go, namespace=None, return_path=False, reverse_path=False):
	"""Compute the maximal depth and deepest path in a GO DAG.

	Depth is defined as the length of the longest path in the direction
	child → parent (more specific → more general).
	Only is_a and part_of relationships are considered for GO→GO edges.
	Parameters
	----------
	go : dict
		Gene Ontology graph (assumed acyclic for GO→GO part).
	namespace : str, optional
		If given, only GO terms in the specified namespace
		('biological_process', 'molecular_function', 'cellular_component')
	return_path : bool, optional
		If True, also return the list of GO terms in the deepest path.
	reverse_path : bool, optional
		If True, the returned path is ordered from most specific term
		to most general term.

	Returns
	-------
	int
		Maximal depth.
	str
		Identifier of the deepest GO term.
	list of str, optional
		If return_path is True, the list of GO term identifiers in the deepest path.
	"""
	dag = induced_goterm_subgraph(go, namespace=namespace)

	if topological_sort(dag) is None:
		return False  # not a DAG!
	order = topological_sort(dag)  # works because GO DAG is acyclic
	depth = {}
	pred  = {}
	for t in nodes(dag):
		depth[t] = 0
		pred[t] = None 


	# Dynamic programming: we found a longer path to parent via child
	for child in reversed(order): # for each node in reverse topological order
		for parent in neighbors(dag, child):
			if depth[child] < depth[parent] + 1:
				depth[child] = depth[parent] + 1
				pred[child] = parent

	max_term = max(depth, key=depth.get)
	max_d = depth[max_term]

	if not return_path:
		return max_d

	path = [max_term]
	cur = max_term
	while pred[cur] is not None:
		cur = pred[cur]
		path.append(cur)

	if not reverse_path:
		path.reverse()
	return max_d, max_term, path
	
	
# Statistics functions

def count_goterm(go, by_namespace=False):
	"""Count GOTerms in the graph.

	Parameters
	----------
	go : dict
		Gene Ontology graph.
	by_namespace : bool, optional
		If False (default), return the total number of GOTerms.
		If True, return a dictionary with counts per namespace
		('biological_process', 'molecular_function', 'cellular_component').

	Returns
	-------
	int or dict
		Total number of GOTerms, or a dict of counts per namespace.
	"""

	if by_namespace: # return counts per namespace
		dico = {
			'biological_process' : 0,
			'molecular_function' : 0,
			'cellular_component' : 0
		}
		for el in dico.keys():
			dico[el] = dico.get(el,0) + len(select_nodes(go , 'namespace', el))
		return dico

	res = 0
	for u in go['nodes']:
		attr = go['nodes'][u]
		if attr.get("type") == "GOTerm":
			res += 1
	return res

def count_geneproducts(go):
	"""Count the number of GeneProduct nodes in the Gene Ontology graph.

	Parameters
	----------
	go : dict
		Gene Ontology graph.

	Returns
	-------
	int
		Number of GeneProduct nodes.
	"""
	return sum([1 for identifiant, attr in go['nodes'].items() if attr.get("type") == "GeneProduct"]) # sum of GeneProducts

def count_annotations(go):
	"""Count the total number of GeneProduct → GOTerm annotation edges.

	Parameters
	----------
	go : dict
		Gene Ontology graph.

	Returns
	-------
	int
		Number of edges with relationship == 'annotation'.
	"""

	n = 0
	for i , v in go['edges'].items():
		for j , att in v.items():
			if att.get('relationship') == 'annotation': # annotation edge
				n+=1 
	return n 

def goterms_per_geneproduct(go,recursive=False):
	"""
	Compute the number of GO terms associated with each GeneProduct.

	Parameters
	----------
	go : dict
	Gene Ontology graph.
	recursive : bool, optional
	If False, count only directly annotated GO terms.
	If True, include ancestor GO terms.

	Returns
	-------
	dict
	GeneProduct ID -> number of associated GO terms.
	"""

	counts = {}

	for node_id , attr in go['nodes'].items(): #
		if attr.get('type')=='GeneProduct': # it's a GeneProduct
			terms = GOTerms(go,node_id,recursive=recursive) # get associated GO terms (recursive or not)
			counts[node_id] = len(terms)
	return counts

def nb_geneproducts_per_goterm_fast(go):
	idx = go_to_gps_index(go) # build index go_id -> set of GeneProducts
	# For each GO term, return the size of the associated GeneProduct set
	return {
	go_id: len(idx.get(go_id, set())) for go_id, attr in go["nodes"].items() if attr.get("type") == "GOTerm"
	}

def evidence_code_distribution(go):
	"""
	Compute the distribution of evidence codes in GO annotations.

	Parameters
	----------
	go : dict
		Gene Ontology graph.

	Returns
	-------
	collections.Counter
		Mapping evidence code -> number of annotations.
	"""
	counter = Counter() # for counting evidence codes
	for u, adj in go['edges'].items():
		for v, attr in adj.items():
			if attr.get('relationship') == 'annotation': # only annotation edges
				for code in attr.get('evidence-codes',[]): # evidence codes
					counter[code] += 1
	return counter

def summary(go):
	"""Compute and print basic statistics about the Gene Ontology graph.

	Statistics include:
	- total number of GOTerms,
	- number of GOTerms per namespace,
	- number of GeneProducts,
	- total number of annotations,
	- mean number of annotations per GeneProduct.

	Parameters
	----------
	go : dict
		Gene Ontology graph.

	Returns
	-------
	dict
		Dictionary summarizing the main counts.
	"""
	nb_gp = count_geneproducts(go)
	nb_ann = count_annotations(go)
	
	sum_go = {
		"total_goterms": count_goterm(go),
		"goterms_by_namespace": count_goterm(go, by_namespace=True),
		"total_gp": count_geneproducts(go),
		"annotations": nb_ann, 
		"avg_ann_per_GP": nb_ann / nb_gp if nb_gp else 0.0,
	}
	res = relationships(go)
	for el in relationships(go).keys():
		sum_go[f"Nombre de relations {el}"] = res[el]
	return sum_go

# First attempt writing GeneProducts function
"""
def go_children(go, go_id):
    Return the direct GO children (more specific terms) of a GO term.

	Children are found as nodes having an edge child → go_id
	with relationship 'is a' or 'part of'.
	
    children = []
    for child in go['edges']:
        if go_id in go['edges'][child]:
            rel = go['edges'][child][go_id].get('relationship')
            if rel in ('is a', 'part of'):
                children.append(child)
    return children


def direct_gp_for(go,go_id):
	
	gps_directs = []
	for u in go['edges']:
		if go_id in go['edges'][u] and is_geneproduct(go,u) and go['edges'][u][go_id].get('relationship') == 'annotation':
			gps_directs.append(u) 
	return gps_directs



def GeneProducts(go, go_id, recursive=False):
	
	if not node_exists(go,go_id) or not is_goterm(go, go_id):
		return []
	
	gps_directs = direct_gp_for(go, go_id)

	if not recursive:
		return sorted(set(gps_directs))

	
	termes = set([go_id])
	attente = [go_id]

	while len(attente)!=0:
		t = attente.pop(0)
		for child in go_children(go , t):
			if child not in termes:
				termes.add(child)
				attente.append(child)
	
	res = set()
	for term in termes:
		for gp in direct_gp_for(go,term):
			res.add(gp)
	return sorted(res)"""


if __name__ == "__main__":
	print("# Gene Ontology module tests")
	