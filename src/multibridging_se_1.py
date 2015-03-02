import random
import re
import doctest
import profile
import time
from mbgraph_1 import Read, Edge, Node, clear

def load_reads_fasta(reads_file, double_stranded):
    #There must be a newline at the end of the file!
    print("Building graph...")

    with open(reads_file) as f:
        for line in f:
            Read.add_read(line, double_stranded)

def run(files, output_dir, double_stranded = False):
    """Runs the algorithm on input files, writing the results to OUTPUT_DIR.
    """
    Read.L = 100
    Node.SIZE_THRESHOLD = Read.L
    Read.K = 15

    clear()
    print(time.asctime())
    load_reads_fasta(files[0], double_stranded)

    Node.kmer_dict = {}
    
    print(time.asctime())
    print("Deleting Suspicious Edges")
    print(str(len(Node.nodes)) + " initial nodes.")
    print(str(len(Edge.no_reads)) + " initial edges.")
    initial_edges= len(Edge.no_reads)
    
    Node.remove_suspicious_edges()
    
    final_edges= len(Edge.no_reads)
    edge_del=initial_edges-final_edges
    print(str(edge_del) + " edges deleted.")

    print(time.asctime())
    print("Condensing graph...")
    print(str(len(Node.nodes)) + " initial nodes.")

    Node.condense_all()
    print(str(len(Node.nodes)) + " nodes after condensing.")

    '''gk'''
    node_hash=0
    f=open('nodes_after_cond','w')
    for node in Node.nodes:
	node.hash=node_hash
	node_hash+=1
    
    for node in Node.nodes:
	print  str(node.hash)+'\t'+str(len(node.bases))+'\t'+str(node.in_degree())+'\t'+str(node.out_degree())+'\n'
	f.write(str(node.hash)+'\t'+node.bases+'\t'+str(node.in_degree())+'\t'+str(node.out_degree())+'\n')
    '''end gk'''
    
    Node.destroy_suspicious()
    print(str(len(Node.nodes)) + " nodes after destroying suspicious nodes.")
    
    
    '''gk'''
    node_hash=0
    f=open('nodes_after_dest','w')
    for node in Node.nodes:
	node.hash=node_hash
	node_hash+=1
    
    for node in Node.nodes:
	print  str(node.hash)+'\t'+str(len(node.bases))+'\t'+str(node.in_degree())+'\t'+str(node.out_degree())+'\n'
	f.write(str(node.hash)+'\t'+node.bases+'\t'+str(node.in_degree())+'\t'+str(node.out_degree())+'\n')

    f.close()
    '''end gk'''

    Node.collapse_all()
    print(str(len(Node.nodes)) + " nodes after collapsing similar nodes.")

    print "Finding bridging reads: %s" % time.asctime()
    Read.find_bridging_reads()
    print "Bridging graph: %s" % time.asctime()
    Node.bridge_all()
    print "%s nodes after bridging: %s" % (str(len(Node.nodes)), time.asctime())

    print("Exporting graph...")
    print(str(len(Node.nodes)) + " final nodes.")

    output_components(output_dir) #Be careful - this will destroy the graph

def output_components(output_dir):
    """Outputs a textual representation of the graph as components, one
    component per set of node/edge/known-paths files. Destroys the graph.

    Also outputs the sequence.
    """
    if output_dir is None:
        return

    component = 0

    with open(output_dir + '/sequence.txt', 'w') as seq_file:
        sequence = Node.sequence()
        i = 0
        while i < len(sequence):
            seq_file.write(sequence[i:i+50] + '\n')
            i += 50

    with open(output_dir + '/single_nodes.txt', 'w') as single_file:
        single_file.write("ID\tBases\tCopycount\tNormalization\n")
	
	all_node_f= open(output_dir+'/all_nodes.txt','w')
        x_node_f=open(output_dir+'/x_nodes.txt','w')
        trip_node_f=open(output_dir+'/trip_nodes.txt','w')
        while len(Node.nodes) > 0:
            with open(output_dir + '/nodes'+str(component)+'.txt', 'w') as nodefile, \
                open(output_dir + '/edges'+str(component)+'.txt', 'w') as edgefile:
                component_nodes = set()
                component_edges = set()
                Node.nodes[0].add_component(component_nodes, component_edges)

                if len(component_nodes) == 1:
                    node = Node.nodes[0]
                    node.hash = -1
                    single_file.write(node.to_string())
                    Node.nodes.remove(node)
                    continue

                node_hash = 0
                nodefile.write("ID\tBases\tCopycount\tNormalization\n")
                for node in component_nodes:
                    node.hash = node_hash
                    node_hash += 1
                    nodefile.write(node.to_string())
                    all_node_f.write(str(node.hash)+'\t'+str(len(node.bases))+'\t'+str(node.in_degree())+'\t'+str(node.out_degree())+'\n')
                    if node.is_xnode_gk():
			x_node_f.write(str(node.hash)+'\t'+str(len(node.bases))+'\t'+str(node.in_degree())+'\t'+str(node.out_degree())+'\n')
		    if node.is_trip_node_gk():
			trip_node_f.write(str(node.hash)+'\t'+str(len(node.bases))+'\t'+str(node.in_degree())+'\t'+str(node.out_degree())+'\n')
                    Node.nodes.remove(node)

                edgefile.write("InID\tOutID\tWeight\tCopycount\tNormalization\n")
                for edge in component_edges:
                    edgefile.write(edge.to_string())

                component += 1
        x_node_f.close()
        trip_node_f.close()
        all_node_f.close()
##gk relevant output for x-node detection              
    #with open(output_dir+'/x_nodes.txt','w') as x_node_f:
	#for n in component_nodes:
	    ##if n.is_xnode_gk():
	    #x_node_f.write(str(n.hash)+'\t'+str(len(n.bases))+'\t'+str(n.in_degree())+'\t'+str(n.out_degree())+'\n')

    #with open(output_dir+'/trip_nodes.txt','w') as trip_node_f:
	#for n in component_nodes:
	    #if n.is_trip_node_gk():
		#trip_node_f.write(str(n.hash)+'\t'+str(len(n.bases))+'\t'+str(n.in_degree())+'\t'+str(n.out_degree())+'\n')
##end gk


def main():
    import sys
    if len(sys.argv) == 1:
        run(['input/hc19_10k_spectrum_1.fasta'], 'output')
    else:
        names = sys.argv[1:]
        run(names[:-1], names[-1])

if __name__ == '__main__':
    main()
