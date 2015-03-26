#!/usr/bin/env python
'''
Created on Mar 20, 2015

@author: GMK
'''



from multibridging_lib import Edge, Kmer, Node, Reads


Kmer.K=5
with open('Kmers','r') as f:
    for line in f:
        kmer_read=line.strip()
        Kmer.add_kmer(kmer_read)
    
Node.condense_graph()
Node.deb_out()