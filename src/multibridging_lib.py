'''
Created on Mar 2, 2015

@author: GMK
'''

class Kmer(object):
    '''
    Kmers added.
    '''
    kmer_dict={}
#     K=5


    def __init__(self, kmer_bases):
        '''
        Constructor
        '''
        kmer_prefix=kmer_bases[:-1]
        kmer_suffix=kmer_bases[1:]
        Kmer.kmer_minus_1_mer_set.add(kmer_prefix)
        Kmer.kmer_minus_1_mer_set.add(kmer_suffix)
        
    
    @staticmethod
    def add_kmer(kmer_bases, reverse_complement=False):
        assert len(kmer_bases)==Kmer.K
        assert kmer_bases not in Kmer.kmer_dict
        
        kmer_prefix=kmer_bases[:-1]
        kmer_suffix=kmer_bases[1:]
        if kmer_prefix not in Kmer.kmer_dict:
            Kmer.kmer_dict[kmer_prefix] = Node(kmer_prefix)

        
        if kmer_suffix not in Kmer.kmer_dict:
            Kmer.kmer_dict[kmer_suffix] = Node(kmer_suffix)
        
        Kmer.kmer_dict[kmer_prefix].link_to(Kmer.kmer_dict[kmer_suffix], Kmer.K-2)

        
       

class Node(object):
    """A node represents a string in the De Bruijn graph."""
    #List of nodes
    nodes = []
    #Dictionary of K-mers at the beginning
    kmer_dict = {}
    #Minimum number of reads that must see an edge between two K-mers for the edge not to be destroyed
    EDGE_READ_THRESHOLD = 2
    
    output_dir='/media/data/Code_clean/repeat_estimation/src/'

    
    def __init__(self, bases):
        self.bases = bases
        self.reads = []
        self.in_edges = []
        self.out_edges = []
        Node.nodes.append(self)
        
        
    def in_degree(self):
        """Return SELF's in-degree.
        """
        return len(self.in_edges)

    def out_degree(self):
        """Return SELF's out-degree.
        """
        return len(self.out_edges)
    
    def link_to(self, next_node, weight):
        """Make a link from SELF to NEXT_NODE. Return the resulting edge.
        """
        edge = Edge(weight, self, next_node)
        return edge

    def link_from(self, previous_node, weight):
        """Make a link from PREVIOUS_NODE to SELF. Return the resulting edge.
        """
        return previous_node.link_to(self, weight)
    
    @staticmethod
    def deb_out():
        g=open(Node.output_dir+'Graph_structure','w')
        h=open(Node.output_dir+'Input_structure','w')
        f=open(Node.output_dir+'Output_structure','w')
        ed=open(Node.output_dir+'Edge_structure','w')
        
        node_hash=0
        for node in Node.nodes:
            node_hash+=1
            node.node_hash=node_hash
        
        edge_hash=0
        for edge in Edge.edges:
            edge_hash+=1
            edge.edge_hash=edge_hash
        
        for node in Node.nodes:
            g.write(str(node.node_hash)+'\t'+str(node.bases)+'\t'+str(node.in_degree())+'\t'+str(node.out_degree())+'\n')
            h.write (str(node.node_hash)+'\t')
            for out_edge in node.out_edges:
                h.write(str(out_edge.out_node.node_hash)+'\t')
            h.write('\n')
            f.write (str(node.node_hash)+'\t')
            for in_edge in node.in_edges:
                f.write(str(in_edge.in_node.node_hash)+'\t')
            f.write('\n')
            
        for edge in Edge.edges:
            ed.write('E'+str(edge.edge_hash)+'\t'+ str(edge.in_node.node_hash)+'\t'+ str(edge.out_node.node_hash)+'\n')
            
        g.close()
        f.close()
        h.close()    
        ed.close()
        
    @staticmethod
    def condense_graph():
        '''Condense the graph'''
        no_condensed=0
        for in_node in Node.nodes:
#             print 'in outer loop'+in_node.bases
#             print str(in_node.out_degree())
            if  not hasattr(in_node, 'destroyed'):
                if in_node.out_degree()==1 and in_node.out_edges[0].out_node.in_degree()==1:
                    no_condensed+=1
                    in_node.destroyed=True
                    in_node.out_edges[0].out_node.destroyed=True
#                     print in_node.bases
                    in_node.out_edges[0].condense()
        
        #Nodes not removed from the list while looping so that the loop invariant is not disturbed.
        Node.remove_destroyed()
        return no_condensed
        
    @staticmethod
    def remove_destroyed():
        Node.nodes=[x for x in Node.nodes if not hasattr(x, 'destroyed')]
    

        
class Edge(object):
    """An edge between nodes.
    """
    edges=set()

    def __init__(self, weight, in_node, out_node):
        self.in_node = in_node
        self.out_node = out_node
        self.weight = weight
        self.in_node.out_edges.append(self)
        self.out_node.in_edges.append(self)
        self.traversed = False
        Edge.edges.add(self)
        
    def predecessors(self):
        """Return a generator for SELF's predecessors.
        """
        return (e.in_node for e in self.in_edges)

    def successors(self):
        """Return a generator for SELF's successors.
        """
        return (e.out_node for e in self.out_edges)
    
    def condense(self):
    
        '''Condenses and edge with where both the in_node and the out_node have degree 1'''
#         pass
        source = self.in_node
        destination = self.out_node
        new_bases = source.bases + destination.bases[self.weight:]
        condensed = Node(new_bases)
         
        if source is destination:
            self.destroy()
            for edge in list(source.out_edges):
                new_e = condensed.link_to(edge.out_node, edge.weight)
                edge.destroy()
            for edge in list(source.in_edges):
                new_e = condensed.link_from(edge.in_node, edge.weight)
                edge.destroy()
     
            return condensed
         
        for edge in list(source.in_edges):
            previous_node = edge.in_node
            previous_node.link_to(condensed, edge.weight)
            edge.destroy()
        for edge in list(destination.out_edges):
            next_node=edge.out_node
            next_node.link_from(condensed, edge.weight)
            edge.destroy()
        
        self.destroy()
        return condensed
#         pass
        
        
    def destroy(self): 
        ''' Destroy an edge'''

        self.in_node.out_edges.remove(self)
        self.out_node.in_edges.remove(self)
        self.in_node = None
        self.out_node = None
        Edge.edges.remove(self)
        

class Reads(object):
    reads=set()
    read_len=10
    
    @staticmethod
    def load_reads(file_name):
        with open (file_name,'r') as read_file:
            for line in read_file:
                read=line.strip()
                assert len(read)==Reads.read_len
                Reads.reads.add(read)