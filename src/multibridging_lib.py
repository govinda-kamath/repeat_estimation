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
    
    def predecessors(self):
        """Return a generator for SELF's predecessors.
        """
        return (e.in_node for e in self.in_edges)

    def successors(self):
        """Return a generator for SELF's successors.
        """
        return (e.out_node for e in self.out_edges)

    def precedes(self, node):
        """Return True if SELF is a predecessor of NODE, false otherwise.
        """
        for edge in self.out_edges:
            if edge.out_node is node:
                return True
        return False

    def succeeds(self, node):
        """Return True if NODE is a predecessor of SELF, false otherwise.
        """
        return node.precedes(self)
    
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
    
    def is_xnode(self):
        return (self.in_degree() >=2) and (self.out_degree() >=2)
    
    @staticmethod
    def all_xnodes():
        return (n for n in Node.nodes if n.is_xnode())
    
    def link_read(self, read, index):
        """Given a read READ that bridges this X-node (contains it from INDEX
        to INDEX+self.bases-1 inclusive), link READ to this X-node.
        """
        self.reads.append((read, index))
        
    def unlink_read(self, read_and_index):
        """Given a tuple READ_AND_INDEX = (read, index) that bridges this
        X-node at index, unlink it.
        """
        self.reads.remove(read_and_index)
        
    def refresh_bridging_reads(self):
        """Check whether bridging reads are still bridging reads, and
        remove if not. Reads can be invalid if they do not fully cover
        the node + 1 base at each edge, or if the bases covered at the
        edge do not correspond to nodes in the graph (e.g. if the
        original nodes have been removed).
        """
        reads = [(read, index) for (read, index) in self.reads
            if index > 0 and
            len(read.bases) > index+len(self.bases) and
            read.bases[index:index + len(self.bases)] == self.bases]
        reads = list(set(reads))
        real_reads = []

        for read, index in reads:
            bridge_in, bridge_out = False, False
            for edge in self.in_edges:
                if read.bases[index-1] == \
                    edge.in_node.bases[len(edge.in_node.bases) - edge.weight - 1]:
                    bridge_in = True
            for edge in self.out_edges:
                if read.bases[index + len(self.bases)] == \
                    edge.out_node.bases[edge.weight]:
                    bridge_out = True
            if bridge_in and bridge_out:
                real_reads.append((read, index))
        self.reads = real_reads
        
    @staticmethod
    def bridge_all():
        """Continue bridging until no bridged X-nodes remain.
        """
        while Node.bridge_some(): 
            pass

    @staticmethod
    def bridge_some():
        """Perform the bridging step on all bridged x-nodes which
        exist at the start of this method.

        Return True if bridging step performed at least once.
        Return False if no bridged X-node could be found.
        """
        #Look for bridged X-nodes
        bridged = [n for n in Node.bridged_xnodes()]
        #print "%s nodes after this run: %s" % (str(len(Node.nodes)), time.asctime())

        if len(bridged) == 0:
            return False
        for node in bridged:
            node.bridging_step()
        Node.remove_destroyed()
        return True

    def bridging_step(self):
        """Perform the bridging step for a single node.
        """
    
        self.refresh_bridging_reads()
        assert len(self.reads) > 0
        assert len(self.in_edges) >= 2 and len(self.out_edges) >= 2
        
        print len(Node.nodes)
        
        u_list, w_list = [], []
        v_back, v_forward, self_loop_weight = None, None, None
        in_edges, out_edges = list(self.in_edges), list(self.out_edges)
        
        print len(Node.nodes)
        
        for in_edge in in_edges:
            u = in_edge.extend_back()
            print 'Extend back' + str(len(Node.nodes))
            if in_edge.in_node is self:
                v_back = u
                self_loop_weight = in_edge.weight
            u_list.append(u)
        
        print len(Node.nodes)   
         
        for out_edge in out_edges:
            w = out_edge.extend_forward()
            print 'Extend for ' + str(len(Node.nodes))
            if out_edge.out_node is self:
                v_forward = w
            w_list.append(w)
        
        print len(Node.nodes)
        
        

        for edge in list(self.in_edges):
            edge.destroy()
        for edge in list(self.out_edges):
            edge.destroy()
          
        Node.deb_out()
          
        print len(Node.nodes)
        
        if v_back:
            assert v_forward is not None
            v_forward.link_to(v_back, self_loop_weight+2)

        for read, index in list(self.reads):
            bases = read.bases[index-1:index+len(self.bases)]
            matched_u = [u for u in u_list if u.bases == bases]
            bases = read.bases[index:index+len(self.bases)+1]
            matched_w = [w for w in w_list if w.bases == bases]
            self.unlink_read((read, index))

            if len(matched_u) != 1 or len(matched_w) != 1:
                continue

            u, w = matched_u[0], matched_w[0]
            u.link_read(read, index - 1)
            w.link_read(read, index)
            if not u.precedes(w):
                u.link_to(w, len(self.bases))
                u.bridged = True
                w.bridged = True
        
        Node.deb_out()
        
        
        
        unbridged_u = [u for u in u_list if not u.bridged]
        unbridged_w = [w for w in w_list if not w.bridged]
        
        
        
        if len(unbridged_u) == 1 and len(unbridged_w) == 1:
            u, w = unbridged_u[0], unbridged_w[0]
            u.link_to(w, len(self.bases))
        else:
            assert len(unbridged_u) + len(unbridged_w) == 0
        
        print len(Node.nodes)
        print len(u_list)
        print len(w_list)
        #Condense edges if necessary
        for n in u_list + w_list:
            for edge in n.in_edges + n.out_edges:
                edge.local_condense()
         
        print len(Node.nodes)
        #Destroy node
        self.destroyed=True
        self.destroy(True)
        
        

    
    @staticmethod
    def bridged_xnodes():
        for n in Node.all_xnodes():
            n.refresh_bridging_reads()
            if len(n.reads) > 0: #Can probably remove this check
                #Find if all in/out-edges are bridged
                in_bases, out_bases = set(), set()
                for read, index in n.reads:
                    in_bases.add(read.bases[index-1])
                    out_bases.add(read.bases[index+len(n.bases)])
                bridged_in = len(n.in_edges) - len(in_bases)
                bridged_out = len(n.out_edges) - len(out_bases)
                if bridged_in == 0 and bridged_out == 0:
                    yield n
                elif bridged_in == 1 and bridged_out == 1:
                    yield n
    
    def local_condense(self):
        """Try to condense every edge incident to SELF, if valid.
        """
        if len(self.out_edges) == 1:
            self.out_edges[0].local_condense()
        if len(self.in_edges) == 1:
            self.in_edges[0].local_condense()

    def destroy(self, pop):
        """Remove SELF from list of nodes.

        Iff POP is false, only mark self for destruction.
        """
        assert len(self.in_edges) == 0, "Remaining links: %s" % self
        assert len(self.out_edges) == 0, "Remaining links: %s" % self
        assert len(self.reads) == 0, "Remaining reads: %s" % self.reads
        if pop:
            Node.nodes.remove(self)
        else:
            self.destroyed = True
    
        
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
                condensed.read=source.reads
            return condensed
         
        for edge in list(source.in_edges):
            previous_node = edge.in_node
            previous_node.link_to(condensed, edge.weight)
            edge.destroy()
        for edge in list(destination.out_edges):
            next_node=edge.out_node
            next_node.link_from(condensed, edge.weight)
            edge.destroy()
        
        #Update reads
        source_reads = list(source.reads)
        destination_reads = [(r, i-(len(source.bases)-self.weight))
            for (r, i)
            in destination.reads]
        destination_reads = [(r, i) for (r, i) in destination_reads
            if (r, i) not in source_reads]
        condensed.reads = source_reads + destination_reads
        source.reads = []
        destination.reads = []
        
        self.destroy()
        source.destroyed=True
        destination.destroyed=True
        return condensed
#         pass
    
    def local_condense(self):
        """Given that SELF may or may not have already been destroyed, try to
        condense SELF if valid. If condensing succeeds, try to condense
        all in- and out-edges of the resulting node.
        """
        if self.in_node is None:
            return
        if len(self.in_node.out_edges) > 1 or len(self.out_node.in_edges) > 1:
            return
  
        condensed = self.condense()
        for edge in condensed.in_edges + condensed.out_edges:
            edge.local_condense()
        
    def destroy(self): 
        ''' Destroy an edge'''
        
        self.in_node.out_edges.remove(self)
        self.out_node.in_edges.remove(self)
        self.in_node = None
        self.out_node = None
        Edge.edges.remove(self)

    def extend_forward(self):
        """If SELF is an edge v => q, construct a node w as v extended forward
        one base along q. Then construct w => q and return w.
        """
        v, q = self.in_node, self.out_node
        extension = q.bases[self.weight]
        w = Node(v.bases + extension)
        w.link_to(q, self.weight+1)
        for (read, index) in v.reads:
            if read.bridges(w, index+1):
                w.link_read(read, index+1)
        w.bridged = False
        return w

    def extend_back(self):
        """If SELF is an edge p => v, construct a node u as v extended back
        one base along p. Then construct p => u and return u.
        """
        p, v = self.in_node, self.out_node
        extension = p.bases[-self.weight-1]
        u = Node(extension + v.bases)
        p.link_to(u, self.weight+1)
        for (read, index) in v.reads:
            if read.bridges(u, index-1):
                u.link_read(read, index-1)
        u.bridged = False
        return u

class Read(object):
    reads={}
    read_len=6
    
    def __init__(self,bases):
        self.bases = bases
        Read.reads[bases] = self
    
    @staticmethod
    def load_reads(file_name):
        with open (file_name,'r') as read_file:
            for line in read_file:
                read_seen=line.strip()
                assert len(read_seen)==Read.read_len
                read_read=Read(read_seen)
    
    @staticmethod
    def find_bridging_reads():
        found=0
        start_sequences = {}
        for node in Node.all_xnodes():
            bases = node.bases[:Kmer.K]
            if bases not in start_sequences:
                start_sequences[bases] = []
            start_sequences[bases].append(node)

        #Look for reads which contain those bases
        for read_bases, read in Read.reads.items():
            for start in range(1, len(read.bases) - Kmer.K):
                node_bases = read_bases[start:start + Kmer.K]
                #If read contains a start sequence and bridges the
                #corresponding X-node(s), link them.
                if node_bases in start_sequences:
                    xnodes = start_sequences[node_bases]
                    #Check for matching X-nodes
                    matches = (x for x in xnodes if read.bridges(x, start))
                    for xnode in matches:
                        xnode.link_read(read, start)
                        found += 1
        print "%s reads found." % found

    
    def bridges(self, node, index):
        """Return True if your bases bridge NODE at INDEX, False otherwise.
        """
        if index <= 0:
            return False
        if len(self.bases) <= index + len(node.bases):
            return False

        return (self.bases[index:index+len(node.bases)] == node.bases)

                