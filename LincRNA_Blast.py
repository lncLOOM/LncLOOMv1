import subprocess
from Bio.Blast import NCBIXML
import networkx as nx
import re
import math
from pulp import *
from operator import itemgetter
from collections import defaultdict, Counter
import random
import itertools
import subprocess
from gurobipy import *
from itertools import combinations
import numpy as np          


def run_blast(outdir,project_name,file_location,suffix):
    blast = True
    try:
        subprocess.call(["makeblastdb", "-in", file_location, "-dbtype", "nucl","-out",outdir+'/'+project_name+"/Run_Files/Blast/Sequences_BLAST_DB"+suffix])
        
        subprocess.call(["blastn", "-task", "blastn", "-word_size", "11","-query", file_location, "-db", outdir+'/'+project_name+"/Run_Files/Blast/Sequences_BLAST_DB"+suffix,
                         "-strand", "plus", "-out", outdir+'/'+project_name+"/Run_Files/Blast/Sequences_blast_results"+suffix+".xml", "-outfmt", "5","-evalue", "1e-5"])
        #subprocess.call(["blastn", "-task", "blastn","-query", file_location, "-db", project_name+"/Run_Files/Blast/Sequences_BLAST_DB","-out", project_name+"/Run_Files/Blast/Sequences_blast_results.xml", "-outfmt", "5"])

    except:
        blast = False
    return blast


def filter_hsps(hsps):
    
    
    #Firstly: if there are hsps that contain overlapping regions - included only the highest scoring HSP for that region
    
    #The hsp list contains a list of each hsp for each top alignment between consecutive layers
    for alignment_hsps in hsps:
        hsps_to_remove = []
        for i in range(len(alignment_hsps)-1):
            hsp_i = alignment_hsps[i]
            if hsp_i in hsps_to_remove:
                continue            
            for j in range(i+1,len(alignment_hsps)):
                hsp_j = alignment_hsps[j]
                
                #Check if hsp_j overlaps with hsp_i
                if (hsp_j[0][0]>=hsp_i[0][0] and hsp_j[0][0]<=hsp_i[0][1]) or (hsp_i[0][0]>=hsp_j[0][0] and hsp_i[0][0]<=hsp_j[0][1]):
                    if not hsp_j in hsps_to_remove:
                        hsps_to_remove.append(hsp_j)
                elif (hsp_j[1][0]>=hsp_i[1][0] and hsp_j[1][0]<=hsp_i[1][1]) or (hsp_i[1][0]>=hsp_j[1][0] and hsp_i[1][0]<=hsp_j[1][1]):
                    if not hsp_j in hsps_to_remove:
                        hsps_to_remove.append(hsp_j)


        for hsp in hsps_to_remove:

            if hsp in alignment_hsps:
                alignment_hsps.remove(hsp)
              
    #Secondly: if one region aligns with mutiple sections to the query/target sequence - include only the highest scoring HSP for that region

    for alignment_hsps in hsps:
        hsps_to_remove = []
        for i in range(len(alignment_hsps)-1):
            hsp_i = alignment_hsps[i]
            if hsp_i in hsps_to_remove:
                continue            
            for j in range(i+1,len(alignment_hsps)):
                hsp_j = alignment_hsps[j]
                
                #Check if hsp_j aligns regions that contradict the order of hsps aligned in hsp_i
                if hsp_i[0][0]>=hsp_j[0][0] and hsp_i[1][0]<hsp_j[1][0]:
                    if not hsp_j in hsps_to_remove:
                        hsps_to_remove.append(hsp_j)
                elif hsp_i[0][0]<=hsp_j[0][0] and hsp_i[1][0]>hsp_j[1][0]:
                    if not hsp_j in hsps_to_remove:
                        hsps_to_remove.append(hsp_j)
        for hsp in hsps_to_remove:
            alignment_hsps.remove(hsp)

    return hsps

def parse_blast(outdir,project_name,file_location,headers,sequences,seq_lengths,intron_indices,suffix):

    #Create lists of ordered sequences, headers and lengths
    ordered_sequences = []
    ordered_headers = []
    ordered_lengths = []
    ordered_introns = []
    ordered_hsps = [] #list of hsps between each ordered layer - if no hsp found then a blank list is return for the pair of consecutive layers

    if 1==1:
    #try:
        #Blast records can only be accessed once - therefore store data into lists
        blast_xml = open(outdir+'/'+project_name+"/Run_Files/Blast/Sequences_blast_results"+suffix+".xml")
        records = NCBIXML.parse(blast_xml)
        query_titles = []
        query_alignments = []
        query_hsps = []
        for record in records:
            query_titles.append(record.query)
            alignment_titles = []
            record_hsps = []
            for alignment in record.alignments:
                alignment_titles.append(int(alignment.hit_id.split('|')[2].strip()))
                alignment_hsps = []
                for hsp in alignment.hsps:
                    alignment_hsps.append([(hsp.query_start,hsp.query_end),(hsp.sbjct_start,hsp.sbjct_end)])
                record_hsps.append(alignment_hsps)
            

            query_alignments.append(alignment_titles)
            query_hsps.append(record_hsps)
        blast_xml.close()
        #for i in range(len(query_titles)):
            #print query_titles[i]
            #print query_alignments[i]
            #print '---------------------'



        #use data in lists to sort sequences and obtain HSPs
        ordered_sequences_index = [0] #list of layers ordered by similarity - original top layer remains top layer
        len_alignments = [len(x) for x in query_alignments[1:]]
        align_by_len = zip(len_alignments, range(1,len(query_titles)))
        align_by_len.sort(key=lambda x: x[0],reverse=True)

        ordered_sequences_index.extend([x[1] for x in align_by_len])
        for i in range(len(query_titles)):
            current_title_index = ordered_sequences_index[i]
            current_alignments = query_alignments[current_title_index]
            ordered_set = ordered_sequences_index[0:i+1]
            No_HSP = True
            for h,best_alignment in enumerate(current_alignments): #use h as an index to obtain the correct HSP for the ordered layers
                if best_alignment in ordered_set or best_alignment==current_title_index:
                    continue
                else:
                    ordered_sequences_index.remove(best_alignment)
                    ordered_sequences_index.insert(i+1, best_alignment)
                    ordered_hsps.append(query_hsps[current_title_index][h])
                    No_HSP = False
                    break
            if No_HSP:
                ordered_hsps.append([])



        for i in ordered_sequences_index:
            ordered_sequences.append(sequences[i])
            ordered_headers.append(headers[i])
            ordered_lengths.append(seq_lengths[i])
            ordered_introns.append(intron_indices[i])
     
        ordered_hsps = filter_hsps(ordered_hsps)
    else:
    #except:
        print "ERROR! BLAST cannot read fasta file of sequences!"
    return ordered_headers,ordered_sequences,ordered_lengths,ordered_introns,ordered_hsps



##########################################################

#Possible alternative to blast

###################################
def solve_lpp(G,method):

    nodes_pruned = list(G.nodes())
    #Get all edges of the graph
    graph_egdes=list(G.edges())
    #Get first and last layers of graph:
    graph_layers = [G.node[x]['layer'] for x in G.nodes()]
    graph_layers.sort()
    #Define and solve the lp problem

    #Define variables (edges are variables, where each edge can = 0 or 1 only)
    vars = LpVariable.dicts("edge",[i for i in G.edges()],0,1,LpInteger)

    #Define a maximisation problem
    prob = LpProblem("problem", LpMaximize)

    #Define the objectiove: maximization of the sum of the edges
    prob += lpSum(vars), "Sum_of_Edges" #ADDS THE OBJECTIVE ATTRIBUTE TO PROB

    #Define the constraints

    for i in nodes_pruned:
        
        if int(G.node[i]['layer'])!=graph_layers[-1] and int(G.node[i]['layer'])!=graph_layers[0]:
            prob += lpSum([vars[(i,a)] for a in G.successors(i)])<= 20*lpSum([vars[(b,i)] for b in G.predecessors(i)])          
            prob += lpSum([vars[(a,i)] for a in G.predecessors(i)])<= 20*lpSum([vars[(i,b)] for b in G.successors(i)])
        #Third constraint: Edges cannot intersect each other across layers - this constraint ensures that the order of the motifs are maintained

        #if the node is not in the last layer
        if float(G.node[i]['layer'])!=graph_layers[-1]:

            #iterate through the successors of that node:            
            for r in G.successors(i):
              
                #iterate through the edges of other nodes in the same layer:
                for x in [e for e in graph_egdes if G.node[e[0]]['layer']==G.node[i]['layer'] if e!=(i,r)]:
                    if G.node[i]['location']<G.node[x[0]]['location']:
                        if G.node[r]['location']>G.node[x[1]]['location']:
                            prob += lpSum(vars[m] for m in [(i,r),x])<=1

                    elif G.node[i]['location']>G.node[x[0]]['location']:
                        if G.node[r]['location']<G.node[x[1]]['location']:
                            prob += lpSum(vars[m] for m in [(i,r),x])<=1

    #solve the lp-problem with the lp-solver of choice
    


    prob.solve(pulp.GUROBI())



    #Construct Solution Graph

    F=nx.MultiDiGraph()
    

    for v in prob.variables():
        #if the edge got a score of 1 (meaning, it is within the solution):
        if v.varValue==1.0:

            #Add both nodes and the edge to the solution graph
            edge=(v.name.split('edge_')[1]).strip('(').strip(')').split(',_')
            node1 = edge[0].strip("'")
            node2 = edge[1].strip("'")
            F.add_node(node1)
            F.add_node(node2)
            F.add_edge(node1,node2)

            F.node[node1]['layer']=G.node[node1]['layer']
            F.node[node1]['seq']=G.node[node1]['seq']
            F.node[node1]['location']=G.node[node1]['location']
            F.node[node1]['end']=G.node[node1]['end']
            F.node[node1]['posx']=G.node[node1]['posx']
            F.node[node1]['posy']=G.node[node1]['posy']
            F.node[node1]['Label']=G.node[node1]['Label']
  
            F.node[node2]['layer']=G.node[node2]['layer']
            F.node[node2]['seq']=G.node[node2]['seq']
            F.node[node2]['location']=G.node[node2]['location']
            F.node[node2]['end']=G.node[node2]['end']
            F.node[node2]['posx']=G.node[node2]['posx']
            F.node[node2]['posy']=G.node[node2]['posy']
            F.node[node2]['Label']=G.node[node2]['Label']

    return F




def order_and_segment(outdir,project_name,headers,sequences,seq_lengths,intron_indices,kmers_len):

    #for each sequence_run pairwise comparison.
    #print to file if wanted
    #order sequences, headers, lens, introns
    totalSequences=len(sequences)
    PairwiseGraphs = {}
    MaxNodes = {}
    kmers_len=10
    print "Initialising Segments"
    for i in range(totalSequences):
        PairwiseGraphs[i] = {}
        MaxNodes[i] = []
    for i in range(totalSequences-1):
        layer1=i+1
        sequence1 = sequences[i]
        seq_len1=len(sequence1)
        for j in range(i+1,totalSequences):
            if i!=j:
                layer2=j+1
                sequence2 = sequences[j]
                seq_len2=len(sequence2)
                #Build graph
                G = nx.DiGraph()
                seq1_nodes = []
                seq2_nodes = []
                for x in range(seq_len1):
                   if x<=(seq_len1-kmers_len):
                       location1 = x+1
                       posx1 = location1/float(seq_len1)
                       seq1 = sequence1[x:x+kmers_len]
                       end1 = location1+len(seq1)-1 
                       nodeLabel = str(layer1)+'_'+str(location1)+'_'+seq1
                       G.add_node(nodeLabel,layer=layer1,location=location1,seq=seq1,end=end1,posx=posx1,posy=layer1,Label=nodeLabel)
                       seq1_nodes.append(nodeLabel)

                for x in range(seq_len2):
                   if x<=(seq_len2-kmers_len):
                       location2 = x+1
                       posx2 = location2/float(seq_len2)
                       seq2 = sequence2[x:x+kmers_len]
                       end2 = location2+len(seq2)-1
                       nodeLabel = str(layer2)+'_'+str(location2)+'_'+seq2
                       G.add_node(nodeLabel,layer=layer2,location=location2,seq=seq2,end=end2,posx=posx2,posy=layer2,Label=nodeLabel)
                       seq2_nodes.append(nodeLabel)

                #Add edges for the segment
                for node1 in seq1_nodes:
                    seq1 = G.node[node1]['seq']
                    successors = [x for x in seq2_nodes if G.node[x]['seq']==seq1]
                    for node2 in successors:
                        G.add_edge(node1,node2)

                
                #Create a dictionary of all kmers defined thus far from all sequences
                dict_kmers = {}
                for node in list(G.nodes()):
                    kmer = node.split('_')[2]
                    dict_kmers[kmer]=[]
                for n in list(G.nodes()):
                    seq = n.split('_')[2]
                    dict_kmers[seq].append(n)

                #Determine which nodes to remove
                nodes_to_remove = []
                prune = 4   
                nodes_to_remove = []
                for k in dict_kmers: #Iterate over KMERS as there are fewer kmers than nodes
                    kmer_nodes = dict_kmers[k]
                    for n in kmer_nodes:
                        if len(list(G.successors(n)))>prune:
                            nodes_to_remove.extend(kmer_nodes)
                            break #As soon as the first node with this KMER has too many successors remove the full list of nodes and move to check next KMER  
                for nr in nodes_to_remove:
                    G.remove_node(nr)
                    
                nodes_to_remove = []
                for nr in G.nodes():
                    if not list(nx.all_neighbors(G,nr)):
                        nodes_to_remove.append(nr)
                           
                for nr in nodes_to_remove:
                    G.remove_node(nr)
                #Solve graph
                
                F = solve_lpp(G,'CBC')
                numNodes = len(F.nodes())
                MaxNodes[i].append((numNodes,j))
                MaxNodes[j].append((numNodes,i)) 
                PairwiseGraphs[i][j]=F
                PairwiseGraphs[j][i]=F
                
  
                  

    
    #Order sequences to maximise layer by layer connections
    #Create lists of ordered sequences, headers and lengths
    ordered_sequences = []
    ordered_headers = []
    ordered_lengths = []
    ordered_introns = []

    
    ordered_indexes = [0] #top layer automatically defined to remain top layer
    i=0
    while len(ordered_indexes)<totalSequences:
        kmerMatches = sorted(MaxNodes[i],reverse=True,key=itemgetter(0))
        no_continuation = True
        k=0
        while no_continuation and k<len(kmerMatches):
            topMatch = kmerMatches[k][1]
            if not topMatch in ordered_indexes:
                #To increase lenght of connection path, check if the topMatch has any matches
                subsequent_layer_matches = sorted([x for x in MaxNodes[topMatch] if x[1] not in ordered_indexes],reverse=True,key=itemgetter(0))
                if subsequent_layer_matches!=0:
                    ordered_indexes.append(topMatch)
                    i = topMatch
                    no_continuation = False
            k+=1
        if no_continuation:#if there is no next layer that allows continuation, then keep the sequence in current position - change will make no difference
            ordered_indexes.append(topMatch)
            i = topMatch    
        

    for i in ordered_indexes:
        ordered_sequences.append(sequences[i])
        ordered_headers.append(headers[i])
        ordered_lengths.append(seq_lengths[i])
        ordered_introns.append(intron_indices[i])
        
    ordered_hsps = [] 
    #Calculate enriched kmer segments between ordered layers
    for i in range(totalSequences-1):
        layer1 = ordered_indexes[i]
        len1 = ordered_lengths[i]
        layer2 = ordered_indexes[i+1]
        len2 = ordered_lengths[i+1]
        F = PairwiseGraphs[layer1][layer2]
        if len(F.nodes())!=0:
            nodes_layer1 = [(F.node[n]['location'],F.node[n]['end'],F.node[n]['seq'],F.node[n]['Label']) for n in F.nodes if F.node[n]['layer']==layer1+1]
            nodes_layer2 = [(F.node[n]['location'],F.node[n]['end'],F.node[n]['seq'],F.node[n]['Label']) for n in F.nodes if F.node[n]['layer']==layer2+1]
            ke_segments1 = get_kmer_enriched_segments(nodes_layer1,len1)
            ke_segments2 = get_kmer_enriched_segments(nodes_layer2,len2)
            hsps = align_segments(ke_segments1,ke_segments2,F)
            ordered_hsps.append(hsps)
        else:
            ordered_hsps.append([]) 
        
       
    print 'Segments completed'
    return ordered_headers,ordered_sequences,ordered_lengths,ordered_introns,ordered_hsps



def get_kmer_enriched_segments(nodes_of_layer,seq_len):
    segments = []
    total_nodes = len(nodes_of_layer)
    #nodes_of_layer is list of tuples: (start,end,seq) of each node
    #order this list according to start position
    nodes_of_layer = sorted(nodes_of_layer,key=itemgetter(0))
    #print nodes_of_layer
    #caluclate distance between each consecutive node
    #only spaces between non-overlapping nodes are considered
    spaced_distances = []
    all_distances =[]
    space = nodes_of_layer[0][0]-1 #include distance to first node in stdev calc
    spaced_distances.append(space)
    all_distances.append(space)
    for i in range(total_nodes-1):
        end_node1 = nodes_of_layer[i][1]
        start_node2 = nodes_of_layer[i+1][0]
        space = (start_node2-end_node1)-1
        all_distances.append(space)
        if space>0:#nodes are not overlapping
            spaced_distances.append(space)
    space = (seq_len-nodes_of_layer[-1][1])-1 #include distance from last node to end of seq in stdev calc
    spaced_distances.append(space)
    all_distances.append(space)
    
    st_dev = np.std(spaced_distances)
    #separate sequence into enriched segments based on standard deviation
    n=0
    while n<total_nodes:
        segment = [nodes_of_layer[n]]
        space = all_distances[n+1]
        while space<=st_dev*1.5 and n+1<total_nodes:
            segment.append(nodes_of_layer[n+1])
            n+=1
            space = all_distances[n+1]

        segments.append(segment)
        n+=1
    return segments

def count_unique(segment,F):
    prev = F.node[segment[0][3]]['seq']
    unique = 1
    for kmer in segment[1:]:
        current = F.node[kmer[3]]['seq']
        if prev!=current:
            unique+=1
            prev=current
    return unique    
    

def align_segments(segments1,segments2,F):
    #align based on edges and return hsps
    aligned_clusters = []
    all_nodes = list(F.nodes())
    for seg in segments1:
        unique1=count_unique(seg,F)#count unique kmers because repeats can result in bias length
        extended_seg2 = []
        for kmer in seg:
            label = kmer[3]
            for seg2 in segments2:
                if seg2[0][3] in list(nx.all_neighbors(F,label)):
                    extended_seg2.extend(seg2)
                    break
        if extended_seg2:
            unique2 = count_unique(extended_seg2,F)
            if unique1==unique2:
                cluster = (seg,extended_seg2)
            elif unique1>unique2:
                extended_seg2=[]
                for n in seg:
                    extended_seg2.extend([(F.node[x]['location'],F.node[x]['end'],F.node[x]['seq'],F.node[x]['Label'])for x in sorted(nx.all_neighbors(F,n[3]),key=itemgetter(0))])
                    cluster = (seg,extended_seg2)
            else:
                seg=[]
                for n in extended_seg2:
                    seg.extend([(F.node[x]['location'],F.node[x]['end'],F.node[x]['seq'],F.node[x]['Label'])for x in sorted(nx.all_neighbors(F,n[3]),key=itemgetter(0))])
                    cluster = (seg,extended_seg2)

            aligned_clusters.append(cluster)
        
             
    aligned_hsps = []
   
    for cluster in aligned_clusters:
        seq1 = (cluster[0][0][0],cluster[0][-1][1])
        seq2 = (cluster[1][0][0],cluster[1][-1][1])
        aligned_hsps.append([seq1,seq2])
    print aligned_hsps
    return aligned_hsps
    


###########################################################







