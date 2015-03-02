import pdb
import random
import doctest
import re
import threading
import time

class Read(object):
    """A read is a segment of DNA of length L that has been read.

    It contains information about the base pairs that were read (seq) as well
    as the nodes that were constructed from the read (nodes).

    Test bridging a repeat (T-GGGTCT-A and A-GGGTCT-C).
    >>> clear()
    >>> Read.L, Read.K = 8, 4
    >>> genome = "CCCAGGCCTGGGTCTAAGGGTCTCAACCCA"
    >>> for read_i in range(len(genome)-Read.L+1):
    ...     _ = Read.add_read(genome[read_i:read_i+Read.L])

    Test correct condensing.
    >>> Node.condense_all()
    >>> Node.nodes.sort(key = lambda n:n.bases)
    >>> [n.bases for n in Node.nodes]
    ['GGGTCT', 'TCTAAGGG', 'TCTCAACCCAGGCCTGGG']

    Test that finding bridging reads works.
    >>> Read.find_bridging_reads()
    >>> reads = Node.nodes[0].reads
    >>> len(reads)
    2
    >>> reads.sort(key=lambda r:r[0].bases)
    >>> [r[0].bases for r in reads]
    ['AGGGTCTC', 'TGGGTCTA']

    Test that bridging step works.
    >>> _ = Node.bridge_all()
    >>> len(Node.nodes)
    1
    >>> [n.bases for n in Node.nodes][0] in genome+genome[4:]+genome[4:]
    True
    """

    #Dictionary of all reads as bases:Read object
    reads = {}
    known_paths = []

    MATE_PAIR_LENGTH = 700 #Includes reads at two ends
    MATED_READS = False

    def __init__(self, bases):
        assert bases not in Read.reads
        
        self.bases = bases
        Read.reads[bases] = self #Where is Read defined ?? #gk

        #If this read is the first in a mate pair,
        #the second read in the pair; otherwise None.
        self.mate = None

        #The base from which this read starts, on the first node traversed.
        self.start_base = None
        #The base on which this read ends, on the last node traversed.
        self.end_base = None

    @staticmethod
    def add_read(read_bases, double = False):
        """Given a read READ_BASES as a string, construct and return a read
        and associated nodes and add them to the database.

        If DOUBLE, also add the reverse complement of this read.

        #Test basic reading
        >>> Read.L, Read.K = 8, 4
        >>> clear()
        >>> _ = Read.add_read('CAGAATCC')
        >>> print(Node.nodes[0])
        [] => CAGA => [AGAA3]
        >>> print(Node.nodes[1])
        [CAGA3] => AGAA => [GAAT3]

        #Try a repeat
        >>> clear()
        >>> Read.L = 15
        >>> _ = Read.add_read('CAGAATCCAGAATTA')
        >>> print(Node.nodes[2])
        [AGAA3] => GAAT => [AATC3, AATT3]
        >>> print(Node.nodes[-1])
        [AATT3] => ATTA => []
        """
        #If already in reads, only change copy count.
        def increment(read_bases):
            r = Read.reads[read_bases]
            return r

        if read_bases in Read.reads:
            read = increment(read_bases)
        else:
            read = Read(read_bases)

        #Construct nodes
        prev_k = None
        for start_i in range(Read.L-Read.K+1):
            k_bases = read_bases[start_i:start_i+Read.K]
            prev_k = Node.add_node(k_bases, prev_k)

        if double:
            Read.add_read(Read.reverse_complement(read_bases))

        return read

    @staticmethod
    def reverse_complement(bases):
        """Return the reverse complement of BASES. Assumes BASES is
        all uppercase.

        >>> Read.reverse_complement("ATCGGGG")
        'CCCCGAT'
        """
        replacements = [('A', 't'), ('T', 'a'), ('C', 'g'), ('G', 'c')]
        for ch1, ch2 in replacements:
            bases = re.sub(ch1, ch2, bases)
        return bases[::-1].upper()

    def bridges(self, node, index):
        """Return True if your bases bridge NODE at INDEX, False otherwise.
        """
        if index <= 0:
            return False
        if len(self.bases) <= index + len(node.bases):
            return False

        return (self.bases[index:index+len(node.bases)] == node.bases)

    @staticmethod
    def find_bridging_reads():
        """Find all reads which bridge an X-node and link them to their
        nodes.
        """
        found = 0
        #Construct set of first K bases of all X-nodes
        start_sequences = {}
        for node in Node.all_xnodes():
            bases = node.bases[:Read.K]
            if bases not in start_sequences:
                start_sequences[bases] = []
            start_sequences[bases].append(node)

        #Look for reads which contain those bases
        for read_bases, read in Read.reads.items():
            for start in range(1, len(read.bases) - Read.K):
                node_bases = read_bases[start:start + Read.K]
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
                
    @staticmethod
    def lay_reads():
        """Lays the reads onto the graph. Reads may lie on multiple positions
        and the positions are not guaranteed to be correct (only match the
        first 25 bases).

        >>> Read.L, Read.K = 8, 4
        >>> r1 = Read.add_read('ABCDEFGH')
        >>> r2 = Read.add_read('EFGHIJKL')
        >>> r3 = Read.add_read('2CDEFGHI')
        >>> Node.condense_all()
        >>> Read.lay_reads()
        >>> [(n.bases, i) for n, i in r1.start_positions]
        [('ABCDE', 0)]
        >>> [(n.bases, i) for n, i in r2.start_positions]
        [('CDEFGHIJKL', 2)]
        >>> [(n.bases, i) for n, i in r3.start_positions]
        [('2CDE', 0)]
        """
        reads_by_bases = {}
        for b, read in Read.reads.items():
            first_k = b[:Read.K]
            if first_k not in reads_by_bases:
                reads_by_bases[first_k] = []
            reads_by_bases[first_k].append(read)
            read.start_positions = []

        for node in Node.nodes:
            for k_i in range(len(node.bases) - Read.K + 1):
                kmer = node.bases[k_i:k_i + Read.K]
                if kmer not in reads_by_bases:
                    continue
                for read in reads_by_bases[kmer]:
                    read.start_positions.append((node, k_i))

    @staticmethod
    def find_mate_pairs():
        """For every two mate-paired reads, if they lie on the graph
        in only one way, reconstruct the read corresponding to how
        they lie on the graph, and add it to the list of reads.

        >>> clear()
        >>> Read.MATE_PAIR_LENGTH = 32
        >>> Read.L, Read.K = 8, 4
        >>> _ = Read.add_read('BLAHBLAH')
        >>> r1 = Read.add_read('HLLOBLAH')
        >>> r2 = Read.add_read('BLAHWRLD')
        >>> r1.mate = r2
        >>> Node.condense_all()
        >>> Read.find_mate_pairs()
        >>> "HLLOBLAHBLAHBLAHBLAHBLAHBLAHWRLD" in Read.reads
        True
        >>> Read.find_bridging_reads()
        >>> Node.bridge_all()
        >>> [n.bases for n in Node.nodes]
        ['HLLOBLAHBLAHBLAHBLAHBLAHBLAHWRLD']
        """
        Read.lay_reads()

        for b, read in Read.reads.items():
            if read.mate is None:
                continue

            match_reads = []
            for node, start_i in read.start_positions:
                match_reads += node.mate_search(start_i, read)
            match_reads = list(set(match_reads))

            if len(match_reads) == 1:
                new_read = match_reads[0]
                if new_read not in Read.reads:
                    _ = Read(new_read)

    def __str__(self):
        return self.bases

class Edge(object):
    """An edge between nodes.
    """

    def __init__(self, weight, in_node, out_node):
        self.in_node = in_node
        self.out_node = out_node
        self.weight = weight
        self.in_node.out_edges.append(self)
        self.out_node.in_edges.append(self)
	self.traversed = False
#	Edge.no_reads[in_node, out_node]=1
	
    def destroy(self):
        self.in_node.out_edges.remove(self)
        self.out_node.in_edges.remove(self)
        self.in_node = None
        self.out_node = None

    def condense(self, pop):
        """Given an unambiguous edge (in_node has only
        this outgoing edge, a nd out_node has only this incoming edge), condense the
        two nodes into one node, inheriting edges and read parentage.

        This method also works on self-loops but will disrupt the graph. For example, if
        there are edges (A => A, A => B), A can be condensed to itself, resulting in
        AA => B.

        >>> clear()
        >>> Read.K = 4
        >>> x = Node('AATC')
        >>> y = Node.add_node('ATCG', x)
        >>> con = x.out_edges[0].condense(True)
        >>> con.bases
        'AATCG'
        >>> [n.bases for n in Node.nodes]
        ['AATCG']
        """
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

            condensed.reads = source.reads
            source.reads = []
            source.destroy(pop)
            return condensed

        #Update edges
        for edge in list(source.in_edges):
            previous_node = edge.in_node
            edge.in_node.link_to(condensed, edge.weight)
            edge.destroy()
        for edge in list(destination.out_edges):
            edge.out_node.link_from(condensed, edge.weight)
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

        #Clean up database
        self.destroy()
        source.destroy(pop)
        destination.destroy(pop)
        return condensed

    def local_condense(self):
        """Given that SELF may or may not have already been destroyed, try to
        condense SELF if valid. If condensing succeeds, try to condense
        all in- and out-edges of the resulting node.
        """
        if self.in_node is None:
            return
        if len(self.in_node.out_edges) > 1 or len(self.out_node.in_edges) > 1:
            return

        condensed = self.condense(False)
        for edge in condensed.in_edges + condensed.out_edges:
            edge.local_condense()

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

    def to_string(self):
        return (str(self.in_node.hash) + "\t"
            + str(self.out_node.hash) + "\t"
            + str(self.weight) + "\n")

    def __str__(self):
        return self.in_node.bases[:500] + '=>' + str(self.weight) + self.out_node.bases[:500]



class Node(object):
    """A Node is a sequence of base pairs that has been pulled from a read.

    It contains information about the sequence itself (bases), any reads that
    bridge it if it is an X-node, and the connecting edges in the graph
    (in_edges and out_edges).

    >>> clear()
    >>> genome = 'CCCAGGCCTAGGGTCTCGCGCAACCCA' #Genome should be cyclic
    >>> Read.L, Read.K = 8, 4
    >>> for read_i in range(len(genome)-Read.L+1):
    ...     _ = Read.add_read(genome[read_i:read_i+Read.L])
    >>> print(Node.nodes[8])
    [CTAG3] => TAGG => [AGGG3]
    >>> Node.condense_all()
    >>> len(Node.nodes)
    1
    >>> [n.bases for n in Node.nodes][0] in genome+genome[4:]+genome[4:]
    True
    """

    #List of nodes
    nodes = []
    #Dictionary of K-mers at the beginning
    kmer_dict = {}
    #Hamming distance required for two nodes to be collapsed
    HAMMING_THRESHOLD = 3
    #Average prevalence required for a K-mer to not be destroyed
    PREVALENCE_THRESHOLD = 2

    def __init__(self, bases):
        self.bases = bases
        self.reads = []
        self.in_edges = []
        self.out_edges = []
        self.norm = 1.0
        #The number of times a K-mer in this node appeared in a read
        self.prevalence = 0.0
        #The number of K-mers condensed together to form this node
        self.count = 1.0

        #Add to database
        Node.nodes.append(self)

    def in_degree(self):
        """Return SELF's in-degree.
        """
        return len(self.in_edges)

    def out_degree(self):
        """Return SELF's out-degree.
        """
        return len(self.out_edges)

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

    def link_to(self, next_node, weight):
        """Make a link from SELF to NEXT_NODE. Return the resulting edge.
        """
        edge = Edge(weight, self, next_node)
        return edge

    def link_from(self, previous_node, weight):
        """Make a link from PREVIOUS_NODE to SELF. Return the resulting edge.
        """
        return previous_node.link_to(self, weight)

    def max_in(self):
        """Return the maximum weight of any in-edge.
        """
        return max([e.weight for e in self.in_edges]+[0])

    def max_out(self):
        """Return the maximum weight of any out-edge.
        """
        return max([e.weight for e in self.out_edges]+[0])
        
    @staticmethod
    def add_node(bases, previous_node):
        """Add a K-mer from a read with bases BASES and linked
        to previous node PREVIOUS_NODE to the database, or add the link to an
        existing node if necessary. Return the created/existing node.
        """
        #Make new node/use existing one
        if bases in Node.kmer_dict:
            n = Node.kmer_dict[bases]
        else:
            if bases in Node.kmer_dict:
                n = Node.kmer_dict[bases]
            else:
                n = Node(bases)
                Node.kmer_dict[bases] = n

        n.prevalence += 1

        #Update previous_node links
        if previous_node and not previous_node.precedes(n):
            n.link_from(previous_node, Read.K-1)

        return n

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

    @staticmethod
    def remove_destroyed():
        Node.nodes = [n for n in Node.nodes if not hasattr(n, 'destroyed')]

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

    def full_destroy(self, pop):
        """Remove all edges and reads and remove SELF from list of nodes.

        Iff POP is false, only mark self for destruction.
        """
        for edge in list(self.in_edges):
            edge.destroy()
        for edge in list(self.out_edges):
            edge.destroy()
        for read in list(self.reads):
            self.unlink_read(read)
        self.destroy(pop)

    def is_xnode(self):
        return len(self.in_edges) >= 2 and len(self.out_edges) >= 2

    '''gk'''
    def is_xnode_gk(self):
        return len(self.in_edges) >= 2 or len(self.out_edges) >= 2
    
    def is_trip_node_gk(self):
        return len(self.in_edges) >= 3 or len(self.out_edges) >= 3
    '''gk end'''
    
    
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
    def condense_all():
        """Condense edges until no more unambiguous edges remain.

        For efficiency reasons, nodes are marked as destroyed while in
        the condensing process, then all destroyed nodes are removed
        at the end.
        """
        condensed = 0
        #Find nodes with only one outgoing edge
        for source_node in (n for n in Node.nodes if len(n.out_edges) == 1):
            assert not hasattr(source_node, 'destroyed')
            #If next node also has only one incoming edge, collapse
            out_edge = source_node.out_edges[0]
            if len(out_edge.out_node.in_edges) == 1 \
                and source_node is not out_edge.out_node:
                out_edge.condense(False)
                condensed += 1
                if condensed % 1000 == 0:
                    pass #print(condensed)
        Node.remove_destroyed()

    @staticmethod
    def destroy_suspicious():
        """After the condensing stage, destroy all suspicious nodes.
        """
        while Node.destroy_some_suspicious(): pass

    @staticmethod
    def destroy_some_suspicious():
        """Destroy all suspicious nodes. This may create new
        suspicious nodes.


        Return True iff there were some suspicious nodes.
        """
        suspicious_nodes = list(n for n in Node.nodes \
            if n.is_suspicious())
        if len(suspicious_nodes) == 0:
            return False
        suspicious_nodes.sort(key = lambda n: n.average_prevalence())

        for node in (n for n in suspicious_nodes
            if not hasattr(n, 'destroyed')):
            adjacent = [e.in_node for e in node.in_edges]
            adjacent += [e.out_node for e in node.out_edges]
            node.full_destroy(False)
            for n in adjacent:
                n.local_condense()
        Node.remove_destroyed()
        return True

    def is_suspicious(self):
        """Return True iff SELF is a suspicious node.

        Suspicious nodes meet the following four conditions:
        1. The average prevalence is less than Node.PREVALENCE_THRESHOLD.
        2. The node size is less than Node.SIZE_THRESHOLD.
        3. This node has predecessors and their average out-degree
        is at least 2.
        4. This node has successors and their average in-degree
        is at least 2. 
        5. If a node has no predecessors or successors for circular genomes #gk
        """
        if self.average_prevalence() >= Node.PREVALENCE_THRESHOLD:
            return False
        if len(self.bases) >= Node.SIZE_THRESHOLD:
            return False
	if (self.in_degree() == 0 and self.out_degree() == 1 ) or (self.out_degree() == 0 and self.in_degree()==1):
            return True #gk works assuming circular genomes
        if self.in_degree() == 0 or self.out_degree() == 0:
            return False 
        total_degree = 0
        predecessors = list(self.predecessors())
        for p in predecessors:
            total_degree += p.out_degree()
        if len(predecessors) > 0:
            average_out = float(total_degree) / len(predecessors)
            if average_out < 2:
                return False

        total_degree = 0
        successors = list(self.successors())    
        for s in successors:
            total_degree += s.in_degree()
        if len(successors) > 0:
            average_in = float(total_degree) / len(successors)
            if average_in < 2:
                return False

        return True

    def average_prevalence(self):
        """Return the average prevalence of K-mers in this node.
        """
        return self.prevalence / self.count
        
    @staticmethod
    def collapse_all():
        """Collapses all similar nodes.
        """
        while Node.collapse_some():
            pass

    @staticmethod
    def collapse_some():
        """Collapses some similar nodes. Returns True if any similar nodes
        existed.
        """
        collapsed = False
        for node in Node.nodes:
            if node.collapse_out():
                collapsed = True

        Node.remove_destroyed()
        return collapsed

    def collapse(self, other):
        """Collapse two similar nodes: one is destroyed and the other gains the
        prevalence of the smaller node.
        """
        if self.prevalence < other.prevalence:
            other.collapse(self)
            return

        self.prevalence += other.prevalence
        other.full_destroy(False)

    def collapse_out(self):
        """If two of this node's successors are similar enough, collapse them.
        Return True if any two nodes were collapsed.
        """
        successors = [n for n in self.successors() if not hasattr(n, 'destroyed')]
        for i in range(len(successors)):
            for j in range(i + 1, len(successors)):
                if successors[i].similar(successors[j]):
                    successors[i].collapse(successors[j])
                    return True
        return False

    def similar(self, other):
        """Returns True iff SELF's bases are of the same length as OTHER and
        the Hamming distance is less than Node.HAMMING_THRESHOLD, and all
        in- and out-edges are identical.
        """
        if hasattr(self, 'destroyed') or hasattr(other, 'destroyed'):
            return False

        if len(self.bases) != len(other.bases):
            return False
        mismatches = len([i for i in range(len(self.bases)) if self.bases[i] != other.bases[i]])
        if mismatches >= Node.HAMMING_THRESHOLD:
            return False
        if set(self.successors()) != set(other.successors()):
            return False
        if set(self.predecessors()) != set(other.predecessors()):
            return False

        return True

    @staticmethod
    def bridged_xnodes():
        """Return a list of all bridged X-nodes.

        Check that full bridge works.
        >>> clear()
        >>> n = Node('ABCD')
        >>> _ = Node('1A').link_to(n, 1)
        >>> _ = Node('2A').link_to(n, 1)
        >>> _ = n.link_to(Node('D1'), 1)
        >>> _ = n.link_to(Node('D2'), 1)
        >>> n.reads = [(Read('1ABCD2'), 1), (Read('2ABCD1'), 1)]
        >>> len(list(Node.bridged_xnodes()))
        1

        Check that unbridged fails.
        >>> Read.reads = {}
        >>> n.reads = [(Read('1ABCD1'), 1), (Read('2ABCD1'), 1)]
        >>> len(list(Node.bridged_xnodes()))
        0

        Check that almost-bridged (1 u/w remaining) works.
        >>> Read.reads = {}
        >>> n.reads = [(Read('2ABCD1'), 1), (Read('X2ABCD1'), 2)]
        >>> len(list(Node.bridged_xnodes()))
        1

        Check bridging nodes with degree > 2.
        >>> _ = Node('3A').link_to(n, 1)
        >>> _ = n.link_to(Node('D3'), 1)
        >>> Read.reads = {}
        >>> n.reads = [(Read('1ABCD1'), 1), (Read('2ABCD1'), 1)]
        >>> n.reads += [(Read('3ABCD2'), 1), (Read('3ABCD3'), 1)]
        >>> len(list(Node.bridged_xnodes()))
        1

        Check partial bridging with greater degree.
        >>> Read.reads = {}
        >>> n.reads = [(Read('1ABCD1'), 1), (Read('2ABCD1'), 1)]
        >>> n.reads += [(Read('3ABCD2'), 1)]
        >>> len(list(Node.bridged_xnodes()))
        0

        Check that reads that don't lie on the graph
        (erroneous reads) don't count.
        >>> Read.reads = {}
        >>> n.reads = [(Read('1ABCD1'), 1), (Read('2ABCD1'), 1)]
        >>> n.reads += [(Read('3ABCD4'), 1)]
        >>> len(list(Node.bridged_xnodes()))
        0
        """
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

    @staticmethod
    def all_xnodes():
        """Return a list of all existing X-nodes.
        """
        return (n for n in Node.nodes if n.is_xnode())

    @staticmethod
    def bridge_all():
        """Continue bridging until no bridged X-nodes remain.
        """
        while Node.bridge_some(): pass

    @staticmethod
    def bridge_some():
        """Perform the bridging step on all bridged x-nodes which
        exist at the start of this method.

        Return True if bridging step performed at least once.
        Return False if no bridged X-node could be found.
        """
        #Look for bridged X-nodes
        bridged = [n for n in Node.bridged_xnodes()]
        print "%s nodes after this run: %s" % (str(len(Node.nodes)), time.asctime())

        if len(bridged) == 0:
            return False
        for node in bridged:
            node.bridging_step()
        Node.remove_destroyed()
        return True

    def bridging_step(node):
        """Perform the bridging step for a single node.
        """
        node.refresh_bridging_reads()
        assert len(node.reads) > 0
        assert len(node.in_edges) >= 2 and len(node.out_edges) >= 2

        u_list, w_list = [], []
        v_back, v_forward, self_loop_weight = None, None, None
        in_edges, out_edges = list(node.in_edges), list(node.out_edges)
        for in_edge in in_edges:
            u = in_edge.extend_back()
            if in_edge.in_node is node:
                v_back = u
                self_loop_weight = in_edge.weight
            u_list.append(u)
        for out_edge in out_edges:
            w = out_edge.extend_forward()
            if out_edge.out_node is node:
                v_forward = w
            w_list.append(w)

        for edge in list(node.in_edges):
            edge.destroy()
        for edge in list(node.out_edges):
            edge.destroy()

        if v_back:
            assert v_forward is not None
            v_forward.link_to(v_back, self_loop_weight+2)

        for read, index in list(node.reads):
            bases = read.bases[index-1:index+len(node.bases)]
            matched_u = [u for u in u_list if u.bases == bases]
            bases = read.bases[index:index+len(node.bases)+1]
            matched_w = [w for w in w_list if w.bases == bases]
            node.unlink_read((read, index))

            if len(matched_u) != 1 or len(matched_w) != 1:
                continue

            u, w = matched_u[0], matched_w[0]
            u.link_read(read, index - 1)
            w.link_read(read, index)
            if not u.precedes(w):
                u.link_to(w, len(node.bases))
                u.bridged = True
                w.bridged = True

        unbridged_u = [u for u in u_list if not u.bridged]
        unbridged_w = [w for w in w_list if not w.bridged]
        if len(unbridged_u) == 1 and len(unbridged_w) == 1:
            u, w = unbridged_u[0], unbridged_w[0]
            u.link_to(w, len(node.bases))
        else:
            assert len(unbridged_u) + len(unbridged_w) == 0

        #Condense edges if necessary
        for n in u_list + w_list:
            for edge in n.in_edges + n.out_edges:
                edge.local_condense()

        #Destroy node
        node.destroy(True)

    def connected(self):
        return len(self.in_edges) > 0 and len(self.out_edges) > 0

    @staticmethod
    def sequence():
        """Return the DNA base sequence obtained by traversing an Eulerian
        cycle of the graph.
        >>> clear()
        >>> Read.K = 1
        >>> a = Node('A')
        >>> b = Node.add_node('B', a)
        >>> c = Node.add_node('C', b)
        >>> d = Node.add_node('D', c)
        >>> _ = d.link_to(c, 0)
        >>> _ = c.link_to(a, 0)
        >>> Node.sequence()
        'ABCDCA'
        """
        start_node, edges = Node.eulerian_cycle()
        sequence = start_node.bases
        for edge in edges:
            sequence += edge.out_node.bases[edge.weight:]
        return sequence

    @staticmethod
    def eulerian_cycle():
        """Return an Eulerian circuit of the node graph as a list of edges.
        """
        for node in Node.nodes:
            assert len(node.in_edges) == len(node.out_edges)

        start_node = Node.nodes[0] #Start of Eulerian cycle
        edges = []

        def path_nodes():
            yield start_node
            for edge in edges:
                yield edge.out_node

        while True:
            #Look for untraversed edge
            splice_point = None
            for index, node in enumerate(path_nodes()):
                untraversed = [e for e in node.out_edges \
                               if not e.traversed]
                if len(untraversed) > 0:
                    splice_point = index
                    splice_node = node
                    next_edge = untraversed[0]
                    next_edge.traversed = True
                    splice_edges = [next_edge]
                    current_node = next_edge.out_node
                    break

            #All edges traversed
            if splice_point is None:
                break

            #Randomly traverse until arriving back at start
            while current_node is not splice_node:
                untraversed = [e for e in current_node.out_edges \
                               if not e.traversed]
                assert len(untraversed) > 0
                next_edge = untraversed[0]
                next_edge.traversed = True
                splice_edges.append(next_edge)
                current_node = next_edge.out_node

            #Splice your path in
            edges = edges[:splice_point] + splice_edges + edges[splice_point:]

        return start_node, edges

    def add_component(self, node_list, edge_list):
        """Add all nodes in the connected component containing SELF to the
        list NODE_LIST, if not already present. Similarly for EDGE_LIST.
        """
        if self in node_list:
            return

        node_list.add(self)

        for e in self.out_edges:
            edge_list.add(e)

        for e in self.out_edges:
            e.out_node.add_component(node_list, edge_list)
        for e in self.in_edges:
            e.in_node.add_component(node_list, edge_list)

    def mate_search(self, start_i, read):
        """Return a list of possible paired-end reads starting from
        START_I on this node and matching READ and its mate pair.

        Basic tests.
        Construct graph.
        >>> Read.K = 3
        >>> Read.L = 12
        >>> Read.MATE_PAIR_LENGTH = 12
        >>> clear()
        >>> _ = Read.add_read("HELLOABCD123")
        >>> _ = Read.add_read("ABCD1234EFGH")
        >>> _ = Read.add_read("ABCD5678EFGH")
        >>> _ = Read.add_read("ABCD5678EFGX")
        >>> _ = Read.add_read("XBCD1234EFGH")
        >>> Node.condense_all()
        >>> sorted([n.bases for n in Node.nodes])
        ['BCD', 'CD1234EF', 'CD5678EF', 'EFG', 'FGH', 'FGX', 'HELLOABC', 'XBC']

        Add mate-paired reads.
        >>> Read.L = 4
        >>> r1 = Read("ABCD")
        >>> r2 = Read("EFGH")
        >>> r1.mate = r2
        >>> n = [n for n in Node.nodes if n.bases == 'HELLOABC'][0]
        >>> n.mate_search(0, r1)
        []
        >>> n.mate_search(5, r1)
        ['ABCD1234EFGH', 'ABCD5678EFGH']
        """
        prefix = self.bases[start_i:]
        return self.extend_mate_search(prefix, read)

    def rmate_search(self, prefix, read):
        """Return a list of possible paired-end reads starting with
        PREFIX and continuing from the start of this node.
        """
        remaining = Read.MATE_PAIR_LENGTH - len(prefix)
        if remaining <= len(self.bases):
            mate_read = (prefix + self.bases)[:Read.MATE_PAIR_LENGTH]
            if Node.mate_matches(mate_read, read):
                return [mate_read]
            else:
                return []

        prefix = prefix + self.bases
        return self.extend_mate_search(prefix, read)

    def extend_mate_search(self, prefix, read):
        """Return a list of possible paired-end reads starting with
        PREFIX and continuing from this nodes' successors.
        """
        results = []
        for edge in self.out_edges:
            next_node = edge.out_node
            next_prefix = prefix[:-edge.weight]
            results += next_node.rmate_search(next_prefix, read)
        return results

    @staticmethod
    def mate_matches(bases, read):
        """Return True iff BASES is a plausible prefix for the
        paired-end path starting with READ and ending with its
        mate pair. Return False otherwise.

        >>> Read.L = 4
        >>> Read.MATE_PAIR_LENGTH = 12
        >>> clear()
        >>> read = Read("ABCD")
        >>> r2 = Read("EFGH")
        >>> read.mate = r2
        >>> Node.mate_matches("", read)
        True
        >>> Node.mate_matches("ABC", read)
        True
        >>> Node.mate_matches("BCD", read)
        False
        >>> Node.mate_matches("ABCDXXXX", read)
        True
        >>> Node.mate_matches("ABCDXXXXEFGH", read)
        True
        >>> Node.mate_matches("ABCDXXXXXXX", read)
        False
        >>> Node.mate_matches("ABCDXXXXXEFGH", read)
        False
        >>> Node.mate_matches("ABCDXXXEFGH", read)
        False
        """
        if len(bases) <= Read.L:
            return read.bases[:len(bases)] == bases
        elif len(bases) <= Read.MATE_PAIR_LENGTH - Read.L:
            return read.bases == bases[:Read.L]
        elif len(bases) <= Read.MATE_PAIR_LENGTH:
            first_match = (read.bases == bases[:Read.L])
            second = bases[Read.MATE_PAIR_LENGTH - Read.L:]
            second_match = (read.mate.bases[:len(second)] == second)
            return (first_match and second_match)
        else:
            return False

    @staticmethod
    def sequence():
        """Return the DNA base sequence obtained by traversing an Eulerian
        cycle of the graph.
        >>> clear()
        >>> Read.K = 1
        >>> a = Node('A')
        >>> b = Node.add_node('B', a)
        >>> c = Node.add_node('C', b)
        >>> d = Node.add_node('D', c)
        >>> _ = d.link_to(c, 0)
        >>> _ = c.link_to(a, 0)
        >>> Node.sequence()
        'ABCDCA'
        """
        for node in Node.nodes:
            if len(node.in_edges) != len(node.out_edges):
                print "Edge mismatch: no Eulerian cycle"
                return ""

        start_node, edges = Node.eulerian_cycle()
        sequence = start_node.bases
        for edge in edges:
            sequence += edge.out_node.bases[edge.weight:]
        return sequence

    @staticmethod
    def eulerian_cycle():
        """Return an Eulerian circuit of the node graph as a list of edges.
        """
        start_node = Node.nodes[0] #Start of Eulerian cycle
        edges = []

        def path_nodes():
            yield start_node
            for edge in edges:
                yield edge.out_node

        while True:
            #Look for untraversed edge
            splice_point = None
            for index, node in enumerate(path_nodes()):
                untraversed = [e for e in node.out_edges \
                               if not e.traversed]
                if len(untraversed) > 0:
                    splice_point = index
                    splice_node = node
                    next_edge = untraversed[0]
                    next_edge.traversed = True
                    splice_edges = [next_edge]
                    current_node = next_edge.out_node
                    break

            #All edges traversed
            if splice_point is None:
                break

            #Randomly traverse until arriving back at start
            while current_node is not splice_node:
                untraversed = [e for e in current_node.out_edges \
                               if not e.traversed]
                assert len(untraversed) > 0
                next_edge = untraversed[0]
                next_edge.traversed = True
                splice_edges.append(next_edge)
                current_node = next_edge.out_node

            #Splice your path in
            edges = edges[:splice_point] + splice_edges + edges[splice_point:]

        return start_node, edges

    def local_condense(self):
        """Try to condense every edge incident to SELF, if valid.
        """
        if len(self.out_edges) == 1:
            self.out_edges[0].local_condense()
        if len(self.in_edges) == 1:
            self.in_edges[0].local_condense()

    def to_string(self):
        return (str(self.hash) + "\t" + self.bases + "\n")

    def __str__(self):
        in_string = "["+", ".join(sorted(e.in_node.bases[:500]+str(e.weight) for e in self.in_edges))+"]"
        out_string = "["+", ".join(sorted(e.out_node.bases[:500]+str(e.weight) for e in self.out_edges))+"]"
        return in_string + " => " + self.bases[:500] + " => " + out_string

def clear():
    Read.reads = {}
    Read.known_paths = []
    Node.kmer_dict = {}
    Node.nodes = []

if __name__ == '__main__':
    doctest.testmod()