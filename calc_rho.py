#!/usr/bin/env python
import dendropy as dpy

# LIMITATION OF THIS CODE: calculates rho between a group of taxa and their most recent common ancestor,
# NOT a user-specified node deeper in the tree.

def calc_rho(tree, taxon_list, occupancy=None):
    """ 
    Returns rho statistic and its standard error (sigma) as described in Saillard et al 2000
    (Greenland Eskimo paper, doi=10.1086/303038, pmid=10924403)

    Parameters:
    tree is a dendropy tree object. You can make a dendropy tree object with:
 
          import dendropy as dpy
          tree = dpy.Tree().get_from_string( <newick_string> , 'newick' )
          #TREE MUST BE ROOTED FOR A RELIABLE RESULT!

    taxon_list is a list of the taxon labels for which rho is being calculated.
    Must be a subset of the taxa in the tree.

    occupancy is an optional dictionary describing how many sequences are
    represented by each taxon label - if not specified this defaults to 1
    for each taxon. If supplied this should be of the format:
          occupancy = { 'taxon_label_1': 2, 'taxon_label_2': 6, 'taxon_label_3': 1 ... } 
    """
    from math import sqrt
    if not occupancy: occupancy = dict(zip(taxon_list,[1 for x in taxon_list])) # Default occupancy of each node is 1
    if len(taxon_list) == 0: return 0,0 # Doesn't make sense to calculate rho for no taxa!
    
    def scan_edges(tree, occupancy, taxon_list): 
        """ Generates ni (n_pendant_taxa) and li (node.edge_length) for each edge by iterating over tree nodes and gathering
            edge information
        """      
        if len(taxon_list)==1: 
            mrca = tree.find_node_with_taxon_label(taxon_list[0]).parent_node
        else: 
            mrca = tree.mrca(taxon_labels=taxon_list)
        for node in mrca.preorder_iter():
            if node == mrca: continue # skip initial node
            # get num of pendant taxa (via tree-traversal - store in var n_pendant_taxa)
            n_pendant_taxa = 0
            for p in node.preorder_iter():
                if p.taxon and str(p.taxon) in taxon_list: # this node is terminal, and the taxon is in our list
                    n_pendant_taxa += occupancy[str(p.taxon)]
            # get edge length directly from node attribute
            yield (n_pendant_taxa, node.edge_length) # return iterator
    
    sum_nl = 0
    sum_nsql = 0
    n_taxa = sum([occupancy[taxon] for taxon in taxon_list])
    for (n,l) in scan_edges(tree, occupancy, taxon_list):
        sum_nl += n*l
        sum_nsql += n*n*l
    rho = sum_nl / n_taxa
    sigma = sqrt((sum_nsql)/n_taxa**2)
    return rho, sigma
