#####################
#LncLOOM coded with recursive segmentation of graph, based on simple paths of conserved kmers
#22 December 2019 created
#Updated: 4 March 2020
#####################

#import all the relevant modules
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

def run_mafft_msa(outdir,project_name):
    print "Running Mafft MSA..."
    mafft = True
    in_file = outdir+"/"+project_name+"/Run_Files/Sequences.fasta"
    out_file = open(outdir+"/"+project_name+"/Run_Files/Mafft.fasta","w",0)
    try:
        subprocess.call(["mafft","--maxiterate","1000","--inputorder","--quiet","--6merpair",in_file],stdout=out_file)        
    except:
        mafft = False
    out_file.close()
    return mafft


def extract_msa(outdir,project_name):
    f = open(outdir+"/"+project_name+"/Run_Files/Mafft.fasta","r")
    seq_data = f.readlines()
    f.close()
    aligned_seqs = []

    #Extract aligned sequences
    for i,line in enumerate(seq_data):
        if line[0]=='>':
            seq = ''
            iterSeq = i+1
            while iterSeq<len(seq_data) and seq_data[iterSeq][0]!= '>':
                seq+=seq_data[iterSeq].strip().upper()
                iterSeq = iterSeq+1
            aligned_seqs.append(seq)

    return aligned_seqs


def calculate_similarity(seq1,seq2):
    match = 0.00
    length = len(seq1)
    for i in range(length):
        base1 = seq1[i]
        base2 = seq2[i]
        if base1==base2:
            match+=1

    match = (match/length)*100
    return match


def getIdenticalSequences(sequences):
    print "Excluding identical sequences..."
    total = len(sequences)

    similar = []
    for i in [x for x in range(total-1) if x not in similar]:
        for j in [x for x in range(i+1,total) if x not in similar]:
            if sequences[i]==sequences[j]:
                similar.append(j)

    return similar


def getSimilarSequences(aligned,similarity_cutoff):
    print 'Excluding sequences with >'+str(similarity_cutoff)+' identity...'
    total = len(aligned)

    similar = []
    for i in [x for x in range(total-1) if x not in similar]:
        for j in [x for x in range(i+1,total) if x not in similar]:
            if aligned[i]==aligned[j]:
                similar.append(j)


    for i in [x for x in range(total-1) if x not in similar]:
        for j in [x for x in range(i+1,total) if x not in similar]:
            seqI = aligned[i]
            seqJ = aligned[j]
            similarity = calculate_similarity(seqI,seqJ)
            if similarity>similarity_cutoff:
                similar.append(j)

    return similar         


def exlcudeSimilarSequences(outdir,project_name,sequences,headers,seq_lengths,intron_indices,mafft,similarity):

    if mafft and similarity<100:
        aligned_seqs = extract_msa(outdir,project_name)
        similar_indexes = getSimilarSequences(aligned_seqs,similarity)
        similar_indexes.sort(reverse=True)
        for si in similar_indexes:
            sequences.pop(si)
            headers.pop(si)
            seq_lengths.pop(si)
            intron_indices.pop(si)
            aligned_seqs.pop(si)
        #rewrite mafft file
        if similar_indexes:
            w = open(outdir+"/"+project_name+"/Run_Files/Mafft.fasta","w",0)
            for i in range(len(aligned_seqs)):
                w.write(headers[i]+'\n')
                seq = aligned_seqs[i]
                countBases = 0
                #write sequence to fasta file with 200 bases per line
                for base in seq:
                    w.write(base)
                    countBases+=1
                    if countBases%200==0:
                        w.write('\n')
                w.write('\n')
            w.close

    elif similarity==100.00 and mafft:
        aligned_seqs = extract_msa(outdir,project_name)
        similar_indexes = getIdenticalSequences(sequences)
        similar_indexes.sort(reverse=True)
        for si in similar_indexes:
            sequences.pop(si)
            headers.pop(si)
            seq_lengths.pop(si)
            intron_indices.pop(si)
            aligned_seqs.pop(si)
        #rewrite mafft file
        if similar_indexes:
            w = open(outdir+"/"+project_name+"/Run_Files/Mafft.fasta","w",0)
            for i in range(len(aligned_seqs)):
                w.write(headers[i]+'\n')
                seq = aligned_seqs[i]
                countBases = 0
                #write sequence to fasta file with 200 bases per line
                for base in seq:
                    w.write(base)
                    countBases+=1
                    if countBases%200==0:
                        w.write('\n')
                w.write('\n')
            w.close



    elif similarity==100.00:
        similar_indexes = getIdenticalSequences(sequences)
        similar_indexes.sort(reverse=True)
        for si in similar_indexes:
            sequences.pop(si)
            headers.pop(si)
            seq_lengths.pop(si)
            intron_indices.pop(si)
    else:
        print "Mafft Sequence Alignment Failed: No similar sequences were excluded"



    


    #Write selected sequences (those marked by # in original file excluded) with introns removed to a fasta file for blast purposes
    #Store file in Run_Files folder generated for the current project
    w= open(outdir+"/"+project_name+"/Run_Files/Sequences.fasta",'w',0) #set buffereing to false so that file is written immediately 
    for i,head in enumerate(headers):
        w.write(head+'\n')
        seq = sequences[i]
        countBases = 0
        #write sequence to fasta file with 200 bases per line
        for base in seq:
            w.write(base)
            countBases+=1
            if countBases%200==0:
                w.write('\n')
        w.write('\n')
    w.close

    max_len = 0
    min_len = 0
    if seq_lengths:
        max_len = max(seq_lengths) 
        min_len = min(seq_lengths)   
    print "\n***************************************\nFASTA FILE PARSED:\nTotal sequences extracted = "+str(len(sequences))+\
          "\nLongest sequence = "+str(max_len)+" bases\nShortest sequence = "+str(min_len)+\
          " bases\n***************************************"

    return headers,sequences,seq_lengths,intron_indices




#function 1 - parses the fasta file input to extract sequences, intron indices, headers and lenghts of sequences
def parse_fasta_file(outdir,project_name,file_location,uracil):

    fasta_file = open(file_location,'r')
    seq_data = fasta_file.readlines()
    fasta_file.close()

    headers = list()
    sequences = list()

    

    #Populate list of headers and list of sequences
    for i,line in enumerate(seq_data):
        if line[0]=='>' and line[1]!='#':
            headers.append(line.strip().upper())
            seq = ''
            iterSeq = i+1
            
            while iterSeq<len(seq_data) and seq_data[iterSeq][0]!= '>':
                seq_line = "".join(seq_data[iterSeq].split())
                seq += seq_line.strip().upper()
                iterSeq = iterSeq+1
            
            if uracil:
                seq = seq.replace('T','U')
            else:
                seq = seq.replace('U','T')
            '''else:#check if user did not upload RNA sequence in which case uracil must be set to true
                if 'U' in seq:
                    print "Uracil detected in sequence - dataset now regarded as RNA"
                    uracil = True'''
            sequences.append(seq)

    #Get indexes of introns and remove introns from sequences for analysis
    
    intron_indices = []
    seq_lengths = []
    for i,seq in enumerate(sequences):
        seq=re.sub('X+','|',seq)
        introns=[]
        seq_exons=''
        for n,base in enumerate(seq):
            if base !='|':
                seq_exons+=base
            else:
                index = n-len(introns)
                introns.append(index)
        sequences[i]=seq_exons
        intron_indices.append(introns)
        seq_lengths.append(len(seq_exons))

    


    #Write selected sequences (those marked by # in original file excluded) with introns removed to a fasta file for mafft/blast purposes
    #Store file in Run_Files folder generated for the current project
    w= open(outdir+"/"+project_name+"/Run_Files/Sequences.fasta",'w',0) #set buffereing to false so that file is written immediately 
    for i,head in enumerate(headers):
        w.write(head+'\n')
        seq = sequences[i]
        countBases = 0
        #write sequence to fasta file with 200 bases per line
        for base in seq:
            w.write(base)
            countBases+=1
            if countBases%200==0:
                w.write('\n')
        w.write('\n')
    w.close()   
    return headers,sequences,seq_lengths,intron_indices,uracil



#function 2 - splits two sequences into segments based on hsps from blast - these segments are used as contraints for edge construction
def splice_hsps(query,subject,hsps):
    hsps.sort(key=lambda x: x[0])

    query_segments=[]
    subject_segments=[]
    

    qstart = 0
    sstart = 0
    for hsp in hsps:
        query_start = hsp[0][0]
        query_end = hsp[0][1]

        subject_start = hsp[1][0]
        subject_end = hsp[1][1]

        non_hsp_seg = query[qstart:query_start-1]
        hsp_seg = query[query_start-1:query_end]
        query_segments.append(non_hsp_seg)
        query_segments.append(hsp_seg)
        qstart=query_end

        non_hsp_seg = subject[sstart:subject_start-1]
        hsp_seg = subject[subject_start-1:subject_end]
        subject_segments.append(non_hsp_seg)
        subject_segments.append(hsp_seg)
        sstart=subject_end
        
    non_hsp_seg = query[qstart:]
    query_segments.append(non_hsp_seg)
    non_hsp_seg = subject[sstart:]
    subject_segments.append(non_hsp_seg)

    return query_segments,subject_segments


#function 3 - constructs the current graph given a set of sequences
def build_graph(seq_subset,lens_subset,hsps,kmers_len,prune):
    print 'Building graph...kmer length = '+str(kmers_len)


    #Coded only for blast (code without blast option later)
    #STAGE ONE: BUILDING

    num_seqs_in_set = len(seq_subset)
    #for s in seq_subset:
        #print len(s)

    G=nx.DiGraph()
    high_repetitive_kmers=[]
    breakpoints = []
    #Construct graph, layer by layer - based on HSPs between consecutive layers
    
    for i in range(num_seqs_in_set-1):
        sequence1 = seq_subset[i]
        seq_len1 = lens_subset[i]        

        sequence2 = seq_subset[i+1]
        seq_len2 = lens_subset[i+1]
        #print "Layer "+str(i+1)+' and Layer '+str(i+2)
 
        hsps_seq1_seq2 = hsps[i]

        

        
        if hsps_seq1_seq2: #if hsps exist between the two layers, splice into segments
            

            seq1_segments,seq2_segments = splice_hsps(sequence1,sequence2,hsps_seq1_seq2)
            
            #Logically the number of segments in each sequence should be the same
            num_hsps_segments = len(seq1_segments)

            range1 = seq_len1-kmers_len+1
            seq1_nodes = [sequence1[x:x+kmers_len] for x in range(range1)]
            all_kmers1 = dict.fromkeys(seq1_nodes,0)

            range2 = seq_len2-kmers_len+1
            seq2_nodes = [sequence2[x:x+kmers_len] for x in range(range2)]
            all_kmers2 = dict.fromkeys(seq2_nodes,0)
            
            #Connect all identical kmers that lie in the aligned segments of the consecutive layers
            current_length_seq1 = 0 
            current_length_seq2 = 0
            for si in range(num_hsps_segments):
                seg1 = seq1_segments[si]
                seg2 = seq2_segments[si]

                len_seg1 = len(seg1)
                len_seg2 = len(seg2)

                range1 = len_seg1-kmers_len+1
                range2 = len_seg2-kmers_len+1
                seq1_nodes = [seg1[x:x+kmers_len] for x in range(range1)]
                seq2_nodes = [seg2[x:x+kmers_len] for x in range(range2)]


                kmers1 = dict.fromkeys(seq1_nodes)
                for kmer in kmers1:
                    kmers1[kmer]=[]

                kmers2 = dict.fromkeys(seq2_nodes)
                for kmer in kmers2:
                    kmers2[kmer]=[]

                for x in range(range1):
                    kmer = seq1_nodes[x]
                    loc = current_length_seq1+x+1
                    G.add_node(str(i+1)+'_'+str(loc)+'_'+kmer,layer=i+1,location=loc,seq=kmer,end=loc+kmers_len-1,posx=(loc+kmers_len-1)/float(seq_len1)) 
                    kmers1[kmer].append(str(i+1)+'_'+str(loc)+'_'+kmer)
                    all_kmers1[kmer]+=1
                
                
                for x in range(range2):
                    kmer = seq2_nodes[x]
                    loc = current_length_seq2+x+1
                    G.add_node(str(i+2)+'_'+str(loc)+'_'+kmer,layer=i+2,location=loc,seq=kmer,end=loc+kmers_len-1,posx=(loc+kmers_len-1)/float(seq_len2)) 
                    kmers2[kmer].append(str(i+2)+'_'+str(loc)+'_'+kmer)
                    all_kmers2[kmer]+=1



                common_kmers = [x for x in kmers1 if x in kmers2]

                for kmer in common_kmers:
                    nodes1 = kmers1[kmer]
                    nodes2 = kmers2[kmer]
                    for n1 in nodes1:
                        for n2 in nodes2:
                            G.add_edge(n1,n2)
                
                        
                current_length_seq1+=len_seg1
                current_length_seq2+=len_seg2
            high_repetitive_kmers.extend([x for x in all_kmers1 if all_kmers1[x]>prune])
            if i+2==num_seqs_in_set:
                high_repetitive_kmers.extend([x for x in all_kmers2 if all_kmers2[x]>prune])
            #print high_repetitive_kmers
            

        else: #else connect any match
            range1 = seq_len1-kmers_len+1
            range2 = seq_len2-kmers_len+1

            seq1_nodes = [sequence1[x:x+kmers_len] for x in range(range1)]
            seq2_nodes = [sequence2[x:x+kmers_len] for x in range(range2)]
            kmers1 = dict.fromkeys(seq1_nodes)
            for kmer in kmers1:
                kmers1[kmer]=[]

            kmers2 = dict.fromkeys(seq2_nodes)
            for kmer in kmers2:
                kmers2[kmer]=[]

            
            for x in range(range1):
                kmer = seq1_nodes[x]
                G.add_node(str(i+1)+'_'+str(x+1)+'_'+kmer,layer=i+1,location=x+1,seq=kmer,end=x+1+kmers_len-1,posx=(x+1+kmers_len-1)/float(seq_len1)) 
                kmers1[kmer].append(str(i+1)+'_'+str(x+1)+'_'+kmer)

            for x in range(range2):
                kmer = seq2_nodes[x]
                G.add_node(str(i+2)+'_'+str(x+1)+'_'+kmer,layer=i+2,location=x+1,seq=kmer,end=x+1+kmers_len-1,posx=(x+1+kmers_len-1)/float(seq_len2)) 
                kmers2[kmer].append(str(i+2)+'_'+str(x+1)+'_'+kmer)

            common_kmers = [x for x in kmers1 if x in kmers2]          
            for kmer in common_kmers:
                nodes1 = kmers1[kmer]
                nodes2 = kmers2[kmer]
                for n1 in nodes1:
                    for n2 in nodes2:
                        G.add_edge(n1,n2)

            high_repetitive_kmers.extend([x for x in kmers1 if len(kmers1[x])>prune])
            if i+2==num_seqs_in_set:
                high_repetitive_kmers.extend([x for x in kmers2 if len(kmers2[x])>prune])     
            #print high_repetitive_kmers

            #remove breakpoints
            remove_breaks = [x for x in breakpoints if x in kmers1]
            for kmer in remove_breaks:
                G.remove_nodes_from(kmers1[kmer])
            remove_breaks = [x for x in breakpoints if x in kmers2]
            for kmer in remove_breaks:
                G.remove_nodes_from(kmers2[kmer])
            breakpoints.extend([x for x in kmers1 if x not in kmers2])
            breakpoints.extend([x for x in kmers2 if x not in kmers1])
    #remove breakpoints between secondlast to last layer

    remove_breaks = [x for x in breakpoints if x in kmers2]
    for kmer in remove_breaks:
        G.remove_nodes_from(kmers2[kmer])

                    
    #prune high_repetitive_kmers
    remove = [x for x in G.nodes() if G.node[x]['seq'] in high_repetitive_kmers]
    G.remove_nodes_from(remove) 
    return G

def prune_graph_level(G,level):
    
    #print 'Pruning Level...'
    #ADDITIONAL PRUNING STEP
    #Remove all nodes that do not have more than 100% depth throughout sequence layers
    F=nx.DiGraph()
    F.add_nodes_from(G.nodes(data=True))
    F.add_edges_from(G.edges())

    #First Remove all nodes that extend past the specified level-faster
    remove_nodes =[n for n in F.nodes() if F.node[n]['layer']>level]
    F.remove_nodes_from(remove_nodes)
  
    FR = F.reverse(copy=False)

    for n in F.nodes():   
        maxSuccessors = max(nx.single_source_shortest_path_length(F,n).values())
        maxPredecessors = max(nx.single_source_shortest_path_length(FR,n).values())            
        if maxSuccessors+maxPredecessors<level-1:
            remove_nodes.append(n)

        
    F.remove_nodes_from(remove_nodes)
    #print "Final Prune Complete"
    return F

def prune_graph_level_plus(G,level):
    
    #print 'Pruning Level...'
    #ADDITIONAL PRUNING STEP
    #Remove all nodes that do not have more than 100% depth throughout sequence layers: keeps nodes greater than current level
    F=nx.DiGraph()
    F.add_nodes_from(G.nodes(data=True))
    F.add_edges_from(G.edges())


    remove_nodes =[]

  
    FR = F.reverse(copy=False)

    for n in F.nodes():   
        maxSuccessors = max(nx.single_source_shortest_path_length(F,n).values())
        maxPredecessors = max(nx.single_source_shortest_path_length(FR,n).values())            
        if maxSuccessors+maxPredecessors<level-1:
            remove_nodes.append(n)

        
    F.remove_nodes_from(remove_nodes)


    #Next remove unconnected paths that have the correct number of S and P fro example 3-4 being pruned at level 2
    remove_nodes =[]  
    #FR = F.reverse(copy=False)

    for n in [x for x in F.nodes() if F.node[x]['layer']>1]:
        
        if len(list(F.predecessors(n)))==0:  
            n_aux_neighbors = nx.single_source_shortest_path(F,n).keys()            
            remove_nodes.extend(n_aux_neighbors)        
    F.remove_nodes_from(remove_nodes)

    #print "Final Prune Complete"
    return F




def hasOverlap(path1,path2,G):

    overlap = False
    
    if len(path2)>len(path1):
        r = len(path1)
    else:
        r = len(path2)
    for i in range(r):
        end1 = G.node[path1[i][0]]['end']
        start2 = G.node[path2[i][0]]['location']
        if start2<=end1:
            overlap=True
    return overlap






#function 4 - lpp solver

def solve_lpp(G,solver):

    F=nx.DiGraph()

    graph_kmers = list(dict.fromkeys([G.node[x]['seq'] for x in G.nodes()]))
        
    if len(graph_kmers)==1:

        nodes_pruned = list(G.nodes())
        #Get all edges of the graph
        graph_egdes=list(G.edges())
        #print len(graph_egdes)
        #Get first and last layers of graph:
        graph_layers = list(dict.fromkeys([G.node[x]['layer'] for x in G.nodes()]))
        graph_layers.sort()
        #Define and solve the lp problem

        #Define variables (edges are variables, where each edge can = 0 or 1 only)
        vars = LpVariable.dicts("edge",[i for i in G.edges()],0,1,LpInteger)

        #Define a maximisation problem
        prob = LpProblem("problem", LpMaximize)

        #Define the objectiove: maximization of the sum of the edges
        prob += lpSum(vars), "Sum_of_Edges" #ADDS THE OBJECTIVE ATTRIBUTE TO PROB


        for i in nodes_pruned:
            #First contraint - each node must have at least one predessors if it has a successor and vice versa         
            if int(G.node[i]['layer'])!=graph_layers[-1] and int(G.node[i]['layer'])!=graph_layers[0]:
                prob += lpSum([vars[(i,a)] for a in G.successors(i)])<= 100*lpSum([vars[(b,i)] for b in G.predecessors(i)])          
                prob += lpSum([vars[(a,i)] for a in G.predecessors(i)])<= 100*lpSum([vars[(i,b)] for b in G.successors(i)])

            #Second contraint - edges must not intersect
            #if the node is not in the last layer
            if float(G.node[i]['layer'])!=graph_layers[-1]:

                #iterate through the successors of that node:            
                for r in G.successors(i):

                    for x in [e for e in graph_egdes if G.node[i]['layer']==G.node[e[0]]['layer'] if G.node[i]['end']>=G.node[e[0]]['location'] and G.node[i]['location']<G.node[e[0]]['location']]:
                        if G.node[i]['end']-G.node[x[0]]['location']!=G.node[r]['end']-G.node[x[1]]['location']:
                            prob += lpSum(vars[m] for m in [(i,r),x])<=1


                    for x in [e for e in graph_egdes if G.node[i]['layer']==G.node[e[0]]['layer'] if G.node[r]['end']>=G.node[e[1]]['location'] and G.node[r]['location']<G.node[e[1]]['location']]:
                        if G.node[i]['end']-G.node[x[0]]['location']!=G.node[r]['end']-G.node[x[1]]['location']:
                            prob += lpSum(vars[m] for m in [(i,r),x])<=1


        #solve the lp-problem with the lp-solver of choice
        if solver=='GUROBI':
            prob.solve(pulp.GUROBI(msg=0))
        elif solver=='CBC':
            prob.solve()
        #Construct Solution Graph
        for v in prob.variables():
            #if the edge got a score of 1 (meaning, it is within the solution):
            if round(v.varValue)==1.0:
 
                #Add both nodes and the edge to the solution graph
                edge=(v.name.split('edge_')[1]).strip('(').strip(')').split(',_')
                node1 = edge[0].strip("'")
                node2 = edge[1].strip("'")
                F.add_node(node1)
                F.add_node(node2)
                F.add_edge(node1,node2)
                for attribute in G.node[node1].keys():
                    F.node[node1][attribute]=G.node[node1][attribute]
                    F.node[node2][attribute]=G.node[node2][attribute]
    else:
        nodes_pruned = list(G.nodes())
        #Get all edges of the graph
        graph_egdes=list(G.edges())
        #print len(graph_egdes)
        #Get first and last layers of graph:
        graph_layers = list(dict.fromkeys([G.node[x]['layer'] for x in G.nodes()]))
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
            #First contraint - each node must have at least one predessors if it has a successor and vice versa         
            if int(G.node[i]['layer'])!=graph_layers[-1] and int(G.node[i]['layer'])!=graph_layers[0]:
                prob += lpSum([vars[(i,a)] for a in G.successors(i)])<= 100*lpSum([vars[(b,i)] for b in G.predecessors(i)])          
                prob += lpSum([vars[(a,i)] for a in G.predecessors(i)])<= 100*lpSum([vars[(i,b)] for b in G.successors(i)])


            #Second contraint - edges must not intersect
            #if the node is not in the last layer
            if float(G.node[i]['layer'])!=graph_layers[-1]:

                #iterate through the successors of that node:            
                for r in G.successors(i):
             
                    all_edges_layer = [e for e in graph_egdes if G.node[e[0]]['layer']==G.node[i]['layer'] if e!=(i,r)]
                    constraint_set = [e for e in all_edges_layer if G.node[i]['location']<G.node[e[0]]['location'] and G.node[r]['location']>G.node[e[1]]['location']]
                    constraint_kmers = [G.node[e[0]]['seq'] for e in constraint_set]
                    constraint_kmers.append(G.node[i]['seq'])
                    constraint_kmers = list(dict.fromkeys(constraint_kmers))
                    if len(constraint_kmers)>1:
                        for x in constraint_set:               
                            prob += lpSum(vars[m] for m in [(i,r),x])<=1


                    constraint_set = [e for e in all_edges_layer if G.node[i]['location']>G.node[e[0]]['location'] and G.node[r]['location']<G.node[e[1]]['location']]    
                    constraint_kmers = [G.node[e[0]]['seq'] for e in constraint_set]
                    constraint_kmers.append(G.node[i]['seq'])
                    constraint_kmers = list(dict.fromkeys(constraint_kmers))
                    if len(constraint_kmers)>1:
                        for x in constraint_set:               
                            prob += lpSum(vars[m] for m in [(i,r),x])<=1


                    #Third contraint - overlapping nodes maintained

                    for x in [e for e in graph_egdes if G.node[i]['layer']==G.node[e[0]]['layer'] if G.node[i]['end']>=G.node[e[0]]['location'] and G.node[i]['location']<G.node[e[0]]['location']]:
                        if G.node[i]['end']-G.node[x[0]]['location']!=G.node[r]['end']-G.node[x[1]]['location']:
                            prob += lpSum(vars[m] for m in [(i,r),x])<=1


                    for x in [e for e in graph_egdes if G.node[i]['layer']==G.node[e[0]]['layer'] if G.node[r]['end']>=G.node[e[1]]['location'] and G.node[r]['location']<G.node[e[1]]['location']]:
                        if G.node[i]['end']-G.node[x[0]]['location']!=G.node[r]['end']-G.node[x[1]]['location']:
                            prob += lpSum(vars[m] for m in [(i,r),x])<=1  

        

        #solve the lp-problem with the lp-solver of choice
        if solver=='GUROBI':
            prob.solve(pulp.GUROBI(msg=0))


        elif solver=='CBC':
            prob.solve()
  

        #Construct Solution Graph
        F=nx.DiGraph()
        for v in prob.variables():
            #if the edge got a score of 1 (meaning, it is within the solution):
            if round(v.varValue)==1.0:

                #Add both nodes and the edge to the solution graph
                edge=(v.name.split('edge_')[1]).strip('(').strip(')').split(',_')
                node1 = edge[0].strip("'")
                node2 = edge[1].strip("'")
                F.add_node(node1)
                F.add_node(node2)
                F.add_edge(node1,node2)
            

                for attribute in G.node[node1].keys():
                    F.node[node1][attribute]=G.node[node1][attribute]
                    F.node[node2][attribute]=G.node[node2][attribute]

    return F


def getComplexPaths(G):

    F=nx.DiGraph()
    F.add_nodes_from(G.nodes(data=True))
    F.add_edges_from(G.edges())

    #Get first and last layers of graph:
    top_nodes = [x for x in F.nodes() if F.node[x]['layer']==1]
    UF = F.to_undirected()
       
    
    Paths = {}    
    while len(top_nodes)>0:

        top_node = top_nodes[0]
        repeat = False
        kmer = F.node[top_node]['seq']
        kmer_edges = []
        

        
        tn_aux_neighbors = nx.single_source_shortest_path(UF,top_node).keys()#gets a list of all auxially neighbors with a path to the top_node
        tn_aux_neighbors.sort(key=lambda x: UF.node[x]['layer'])

        i=0
        while i<len(tn_aux_neighbors):
            aux = tn_aux_neighbors[i]
            kmer_edges.extend(list(F.edges(aux)))
            successors = len(list(F.successors(aux)))
            predecessors = len(list(F.predecessors(aux)))
            if successors>1 or predecessors>1:
                repeat = True
            i+=1

        if repeat:
            #Allow for more than one path of the same kmer that has been separated by a path of non-identical kmer
            if kmer in Paths:
                Paths[kmer].append(kmer_edges)
            else:
                Paths[kmer] = [kmer_edges]
        UF.remove_nodes_from(tn_aux_neighbors)
        for x in [n for n in top_nodes if n in tn_aux_neighbors]:
            top_nodes.remove(x)
      
    return Paths

#function 5 - calculate all simple paths (paths without repeats) of a graph

def getSimplePaths(G):
    
    #create copy of mutable graph G
    F=nx.DiGraph()
    F.add_nodes_from(G.nodes(data=True))
    F.add_edges_from(G.edges())


    top_nodes = [x for x in F.nodes() if F.node[x]['layer']==1]
    UF = F.to_undirected()
       

    Paths = {}    
    while len(top_nodes)>0:

        top_node = top_nodes[0]
        repeat = False
        kmer = F.node[top_node]['seq']
        kmer_edges = []
        
        tn_aux_neighbors = nx.single_source_shortest_path(UF,top_node).keys()#gets a list of all auxially neighbors with a path to the top_node
        tn_aux_neighbors.sort(key=lambda x: UF.node[x]['layer'])
        i=0
        while i<len(tn_aux_neighbors) and repeat==False:
            aux = tn_aux_neighbors[i]
            kmer_edges.extend(list(F.edges(aux)))
            successors = len(list(F.successors(aux)))
            predecessors = len(list(F.predecessors(aux)))
            if successors>1 or predecessors>1:
                repeat = True
            i+=1
        if not repeat:
            #Allow for more than one path of the same kmer that has been separated by a path of non-identical kmer
            if kmer in Paths:
                Paths[kmer].append(kmer_edges)
            else:
                Paths[kmer] = [kmer_edges]
        UF.remove_nodes_from(tn_aux_neighbors)
        for x in [n for n in top_nodes if n in tn_aux_neighbors]:
            top_nodes.remove(x)
            
    return Paths




#function 6 - splits a graph into multiple graphs based on simple paths - edges were initally established based on hsps - this contrainst is implicitly kept here
def split_graph(graph,simple_paths,CN,LongCN):
    #print '*****'
    #print simple_paths
    #print LongCN.nodes()
    #print '*****'
    #print "SPLITTING"
    #Convert complex_paths to a graph for easy split
    #graph = nx.DiGraph()
    #graph.add_nodes_from(main_graph.nodes(data=True))
    #graph.add_edges_from(main_graph.edges())
    
  
    #graph.add_nodes_from(LongCN.nodes(data=True))
    #graph.add_edges_from(LongCN.edges())
    
    #graph.add_nodes_from(CN.nodes(data=True))
    #graph.add_edges_from(CN.edges())
 
    graph_layers = list(dict.fromkeys([graph.node[x]['layer'] for x in graph.nodes()]))
    #extract all paths from the simple_paths (dict keyed by k-mer)
    all_simple_paths = []
    for kmer in simple_paths:
        for path in simple_paths[kmer]:
            all_simple_paths.append(path)

       
    all_simple_paths.sort(key=lambda x: graph.node[x[0][0]]['location'])
    #all_simple_paths.sort(key=lambda x: int(x[0][0].split('_')[1]))
    total_paths = len(all_simple_paths)
    overlapped = []
    i=0
    while i<total_paths:
        path1 = all_simple_paths[i]
        cluster=[path1]
        while i+1<total_paths and hasOverlap(path1,all_simple_paths[i+1],graph):
            cluster.append(all_simple_paths[i+1])
            i+=1    
            path1 = all_simple_paths[i]
        overlapped.append(cluster)            
        i+=1
    
       
    #overlapped list contains lists of paths that overlap - these lists are ordered - so can easily extract indices where graph must be split

    #graph_layers = list(dict.fromkeys([graph.node[x]['layer'] for x in graph.nodes()]))
    graph_layers.sort()
    split_indices = {l:[[1],[]] for l in graph_layers}
    split_indices_CN = {l:[[1],[]] for l in graph_layers}
    number_of_subgraphs = len(overlapped)+1




    for cluster in overlapped:
        first_path = cluster[0]
        last_path = cluster[-1]
        for edge in first_path:
            if edge[0] in graph.nodes():
                end = graph.node[edge[0]]['location']-1
                layer = graph.node[edge[0]]['layer']
                split_indices[layer][1].append(end)
                split_indices_CN[layer][1].append(end)
                split_indices_CN[layer][0].append(end+1)
        for edge in last_path:  
            if edge[0] in graph.nodes():
                start = graph.node[edge[0]]['end']+1
                #startCN = graph.node[edge[0]]['location']+1
                layer = graph.node[edge[0]]['layer']
                split_indices[layer][0].append(start)
                #split_indices_CN[layer][0].append(startCN)

    #Construct Sub-graphs based on boundaries
    #list of all segmented graphs
    subgraph_array =[]


    #print split_indices_CN 
    #print split_indices
    #print overlapped
    #print CN.nodes

    for i in range(number_of_subgraphs-1):
        G=nx.DiGraph()
        cp = nx.DiGraph()
        longer_cp = nx.DiGraph() #graph of complex paths from previous graph contructed of longer kmers
        for layer in split_indices:
            s = split_indices[layer][0][i]
            e = split_indices[layer][1][i]
            sc = split_indices_CN[layer][0][i]
            ec = split_indices_CN[layer][1][i]
            G.add_nodes_from([x for x in graph.nodes() if graph.node[x]['layer']==layer and graph.node[x]['location']>=s and graph.node[x]['end']<=e])
            #G.add_nodes_from([x for x in CN.nodes() if CN.node[x]['layer']==layer and CN.node[x]['location']>=sc and CN.node[x]['location']<=ec])
            #G.add_nodes_from([x for x in LongCN.nodes() if LongCN.node[x]['layer']==layer and LongCN.node[x]['location']>=sc and LongCN.node[x]['location']<=ec])
            cp.add_nodes_from([x for x in CN.nodes() if CN.node[x]['layer']==layer and CN.node[x]['location']>=sc and CN.node[x]['location']<=ec])
            longer_cp.add_nodes_from([x for x in LongCN.nodes() if LongCN.node[x]['layer']==layer and LongCN.node[x]['location']>=sc and LongCN.node[x]['location']<=ec])
        G.add_edges_from([x for x in graph.edges() if G.has_node(x[0]) and G.has_node(x[1])])
        cp.add_edges_from([x for x in CN.edges() if cp.has_node(x[0]) and cp.has_node(x[1])])
        longer_cp.add_edges_from([x for x in LongCN.edges() if longer_cp.has_node(x[0]) and longer_cp.has_node(x[1])])

        for x in cp.nodes():
            for attribute in CN.node[x].keys():
                cp.node[x][attribute]=CN.node[x][attribute]
 
        #add edges that continue past the current level
        last_cp = [x for x in cp.nodes() if cp.node[x]['layer']==layer]
        while last_cp:
            for x in last_cp:
                successors = list(CN.successors(x))
                for s in successors:
                    cp.add_node(s,layer=CN.node[s]['layer'],location=CN.node[s]['location'],seq=CN.node[s]['seq'],end=CN.node[s]['end'],posx=CN.node[s]['posx'])
                    cp.add_edge(x,s)
            layer+=1          
            last_cp = [x for x in cp.nodes() if cp.node[x]['layer']==layer]


        for x in longer_cp.nodes():
            for attribute in LongCN.node[x].keys():
                longer_cp.node[x][attribute]=LongCN.node[x][attribute]
 
        #add edges that continue past the current level
        last_cp = [x for x in longer_cp.nodes() if longer_cp.node[x]['layer']==layer]
        while last_cp:
            for x in last_cp:
                successors = list(LongCN.successors(x))
                for s in successors:
                    longer_cp.add_node(s,layer=LongCN.node[s]['layer'],location=LongCN.node[s]['location'],seq=LongCN.node[s]['seq'],end=LongCN.node[s]['end'],posx=LongCN.node[s]['posx'])
                    longer_cp.add_edge(x,s)
            layer+=1          
            last_cp = [x for x in longer_cp.nodes() if longer_cp.node[x]['layer']==layer]

        
        #cp=getNonIntersectingPaths(cp)
        if len(G.nodes())>0 or len(cp.nodes())>0 or len(longer_cp.nodes())>0:
            subgraph_array.append([G,cp,longer_cp])


            
            
    i+=1
    G=nx.DiGraph()
    cp = nx.DiGraph()
    longer_cp = nx.DiGraph()
    for layer in split_indices:
        s = split_indices[layer][0][i]
        sc = split_indices_CN[layer][0][i]
        G.add_nodes_from([x for x in graph.nodes() if graph.node[x]['layer']==layer and graph.node[x]['location']>=s])
        #G.add_nodes_from([x for x in CN.nodes() if CN.node[x]['layer']==layer and CN.node[x]['location']>=sc])
        #G.add_nodes_from([x for x in LongCN.nodes() if LongCN.node[x]['layer']==layer and LongCN.node[x]['location']>=sc])
        cp.add_nodes_from([x for x in CN.nodes() if CN.node[x]['layer']==layer and CN.node[x]['location']>=sc])
        longer_cp.add_nodes_from([x for x in LongCN.nodes() if LongCN.node[x]['layer']==layer and LongCN.node[x]['location']>=sc])
    G.add_edges_from([x for x in graph.edges() if G.has_node(x[0]) and G.has_node(x[1])])
    cp.add_edges_from([x for x in CN.edges() if cp.has_node(x[0]) and cp.has_node(x[1])])
    longer_cp.add_edges_from([x for x in LongCN.edges() if longer_cp.has_node(x[0]) and longer_cp.has_node(x[1])])
    
        
    for x in cp.nodes():
        for attribute in CN.node[x].keys():
            cp.node[x][attribute]=CN.node[x][attribute] 
    #add edges that continue past the current level
    last_cp = [x for x in cp.nodes() if cp.node[x]['layer']==layer]
    while last_cp:
        for x in last_cp:
            successors = list(CN.successors(x))
            for s in successors:
                cp.add_node(s,layer=CN.node[s]['layer'],location=CN.node[s]['location'],seq=CN.node[s]['seq'],end=CN.node[s]['end'],posx=CN.node[s]['posx'])
                cp.add_edge(x,s)
        layer+=1          
        last_cp = [x for x in cp.nodes() if cp.node[x]['layer']==layer]


    for x in longer_cp.nodes():
        for attribute in LongCN.node[x].keys():
            longer_cp.node[x][attribute]=LongCN.node[x][attribute] 
    #add edges that continue past the current level
    last_cp = [x for x in longer_cp.nodes() if longer_cp.node[x]['layer']==layer]
    while last_cp:
        for x in last_cp:
            successors = list(LongCN.successors(x))
            for s in successors:
                longer_cp.add_node(s,layer=LongCN.node[s]['layer'],location=LongCN.node[s]['location'],seq=LongCN.node[s]['seq'],end=LongCN.node[s]['end'],posx=LongCN.node[s]['posx'])
                longer_cp.add_edge(x,s)
        layer+=1          
        last_cp = [x for x in longer_cp.nodes() if longer_cp.node[x]['layer']==layer]


   
    #cp=getNonIntersectingPaths(cp)
    if len(G.nodes())>0 or len(cp.nodes())>0 or len(longer_cp.nodes())>0:
         subgraph_array.append([G,cp,longer_cp])

    #add attributes
    for sub in subgraph_array:
        G=sub[0]
        for x in G.nodes():
            for attribute in graph.node[x].keys():
                G.node[x][attribute]=graph.node[x][attribute]
        '''for x in [n for n in G.nodes() if n in graph.nodes()]:
            for attribute in graph.node[x].keys():
                G.node[x][attribute]=graph.node[x][attribute]
        for x in [n for n in G.nodes() if n in LongCN.nodes()]:
            for attribute in LongCN.node[x].keys():
                G.node[x][attribute]=LongCN.node[x][attribute]'''


    
    return subgraph_array

def get_neighbourhoods(G):

    graph_layers = list(dict.fromkeys([G.node[x]['layer'] for x in G.nodes()]))
    #create copy of mutable graph G
    F=nx.DiGraph()
    F.add_nodes_from(G.nodes(data=True))
    F.add_edges_from(G.edges())


    top_nodes = [x for x in F.nodes() if F.node[x]['layer']==1]
    top_nodes.sort(key=lambda x: F.node[x]['location'])

    UF = F.to_undirected()
    paths = []
      
    while len(top_nodes)>0:
        
        top_node = top_nodes[0]
        nodes = [top_node]
        
        tn_aux_neighbors = nx.single_source_shortest_path(UF,top_node).keys()#gets a list of all auxially neighbors with a path to the top_node
        UF.remove_nodes_from(tn_aux_neighbors)

        for x in [n for n in top_nodes if n in tn_aux_neighbors]:
            top_nodes.remove(x)
        while top_node in tn_aux_neighbors:
            tn_aux_neighbors.remove(top_node)

        tn_aux_neighbors.sort(key=lambda x: (F.node[x]['layer'],F.node[x]['location']))
        nodes.extend(tn_aux_neighbors)
        paths.append(nodes)


    #combined overlapping paths into neighbourhoods

    neighbourhoods = {}
    ni = 1
    i = 0
    total = len(paths)
    while i<total:
        neighbours = [paths[i]]
        tandems = [x for x in paths[i] if G.node[x]['layer']==1]
        end = G.node[tandems[-1]]['end']
        depth = G.node[paths[i][-1]]['layer']
        j = i+1
        break_overlap = True
        while j<total and break_overlap:
            if end>=G.node[paths[j][0]]['location']:
                neighbours.append(paths[j])
                tandems = [x for x in paths[j] if G.node[x]['layer']==1]
                if G.node[tandems[-1]]['end']>end:
                    end = G.node[tandems[-1]]['end']
                if G.node[paths[j][-1]]['layer']>depth:
                    depth = G.node[paths[j][-1]]['layer']
                i+=1
            else:
                break_overlap = False 
            j+=1
        i+=1
        neighbourhoods[ni]=(neighbours,depth)
        ni+=1


    processed_neighbourhoods = {ni:{} for ni in neighbourhoods}
    #neighbourhood_paths = {ni:{} for ni in neighbourhoods}
    for nh in neighbourhoods:
        depth = neighbourhoods[nh][1]
        #neighbourhood_paths[nh]={d:[] for d in [x for x in graph_layers if x>1 and x<=depth]}

    #process neighbourhoods
    for nh in neighbourhoods:
        paths = neighbourhoods[nh][0]
        depth = neighbourhoods[nh][1]

        nodes_per_layer = {l:[] for l in [x for x in graph_layers if x<=depth]}

        processed_neighbourhoods[nh]['Depth']=depth
        for path in paths:
            path_depth = G.node[path[-1]]['layer']
            #np = []
            for n in path:
                layer = G.node[n]['layer']
                nodes_per_layer[layer].append(n)
                #np.append((G.node[n]['layer'],G.node[n]['location'],G.node[n]['end'],G.node[n]['seq']))
            #for pd in [x for x in  graph_layers if x>1 and x<=path_depth]: 
                #neighbourhood_paths[nh][pd].append(np)

        for layer in nodes_per_layer:

            nodes = nodes_per_layer[layer]
            nodes.sort(key=lambda x: G.node[x]['end'])
            limit2 = G.node[nodes[-1]]['end']
            nodes.sort(key=lambda x: (G.node[x]['location'],len(G.node[x]['seq'])))
            limit1 = G.node[nodes[0]]['location']
            combined_nodes = []

            i = 0
            total = len(nodes)
            while i<total:
                cn = G.node[nodes[i]]['seq']
                start = G.node[nodes[i]]['location']
                end = G.node[nodes[i]]['end'] 
                j = i+1
                break_overlap = True
                while j<total and break_overlap:
                    if end>=G.node[nodes[j]]['location'] and G.node[nodes[j]]['end']>end:
                        cn = cn[0:G.node[nodes[j]]['location']-start]+G.node[nodes[j]]['seq']
                        end = G.node[nodes[j]]['end']
                        i+=1
                    elif end>=G.node[nodes[j]]['location'] and G.node[nodes[j]]['end']<=end:
                        i+=1
                    else:
                        break_overlap = False 

                    j+=1
                i+=1
                combined_nodes.append((start,end,cn))

            processed_neighbourhoods[nh][layer] = (combined_nodes,limit1,limit2)

    return processed_neighbourhoods





# function getNonIntersectingPaths - returns a dictionary of all paths in a graph
def getNonIntersectingPaths(F):
    #Create a copy of G - G is mutable

    G=nx.DiGraph()
    G.add_nodes_from(F.nodes(data=True))
    G.add_edges_from(F.edges())
    UF = G.to_undirected()

    Paths = {}

    #print G.nodes()
    #determine the top layer of the graph
    if len(G.nodes())>0:
        layerID = sorted(list(dict.fromkeys([G.node[x]['layer'] for x in G.nodes()])))[0] 
    while len(G.nodes())>0:
        layerID = sorted(list(dict.fromkeys([G.node[x]['layer'] for x in G.nodes()])))[0] 
        top_node = [x for x in G.nodes() if G.node[x]['layer']==layerID][0]
        kmer = G.node[top_node]['seq']
        kmer_edges = []
        tn_aux_neighbors = nx.single_source_shortest_path(UF,top_node).keys()

        for aux in tn_aux_neighbors:
            kmer_edges.extend(list(G.edges(aux)))

        #Allow for more than one path of the same kmer that has been separated by a path of non-identical kmer
        if kmer in Paths:
            Paths[kmer].append(kmer_edges)
        else:
            Paths[kmer] = [kmer_edges]
        G.remove_nodes_from(tn_aux_neighbors)        
    return Paths





def if_intersect(path1,path2):
    intersect = False
    for edge1 in path1:
        layer1=int(edge1[0].split('_')[0])
        pos11 = int(edge1[0].split('_')[1])
        end11 = pos11+len(edge1[0].split('_')[2])-1

        layer2=int(edge1[1].split('_')[0])
        pos12 = int(edge1[1].split('_')[1])
        end12 = pos12+len(edge1[1].split('_')[2])-1

        path2_edges = [e for e in path2 if int(e[0].split('_')[0])==layer1 and int(e[1].split('_')[0])==layer2]
        for edge2 in path2_edges:
            pos21 = int(edge2[0].split('_')[1])
            pos22 = int(edge2[1].split('_')[1])
            end21 = pos21+len(edge2[0].split('_')[2])-1
            end22 = pos22+len(edge2[1].split('_')[2])-1

            if pos11<pos21 and pos12>pos22:
                intersect = True
            elif pos11<pos21 and end11>=pos21:
                intersect = True
            elif pos12<pos22 and end12>=end22:
                intersect = True          
            elif pos11>pos21 and pos12<pos22: 
                intersect = True
            elif pos21<pos11 and end21>=pos11:
                intersect = True
            elif pos22<pos12 and end22>=pos12:
                intersect = True

            #also cheack for overlapping nodes between new paths and complex paths            
            elif pos11<=pos21 and pos21<=end11:
                intersect = True
            elif pos21<=pos11 and pos11<=end21:
                intersect = True
            elif pos12<=pos22 and pos22<=end12:
                intersect = True
            elif pos22<=pos12 and pos12<=end22:
                intersect = True

    return intersect 


#function to calculate and return an cluster of paths that intersect an indiviual deeper path
def getLocalOptimsedPath(complex_path,new_paths):

    Cluster = {'Deep':complex_path,'Shallow':{}}
    deep_path = complex_path.values()[0][0]  
    for kmer in new_paths:
        for path in new_paths[kmer]:
            intersect = if_intersect(path,deep_path)
            if intersect:
                if kmer in Cluster['Shallow']:
                    Cluster['Shallow'][kmer].append(path)
                else:
                    Cluster['Shallow'][kmer]=[path]

    return Cluster

#function to calculate and return local clusters of previous complex_paths with new shallower paths
def getLocalIntersectingPaths(deep_paths,new_paths):

    #print deep_paths

    #print deep_paths
    #print new_paths

    #Hierarchical clustering
    #First cluster deep nodes based on which shallow nodes they intersect.
    Unique_deep_paths = []#list of deep paths that have repeats but do not intersect a new shallow path
    
    ClustersShallow = []
    ClustersDeep = []
    for kmer in new_paths:
        for path in new_paths[kmer]:
            deep_intersects = []
            for deep_kmer in deep_paths:
                for deep_path in deep_paths[deep_kmer]:
                    intersect = if_intersect(path,deep_path)
                    if intersect:
                        deep_intersects.append(deep_path)            
            ClustersShallow.append([path])
            ClustersDeep.append(deep_intersects)

    #Update Unique_deep_paths
    for deep_kmer in deep_paths:
        for deep_path in deep_paths[deep_kmer]:
            intersect = False
            for cluster in ClustersDeep:
                if deep_path in cluster:
                    intersect = True
                    break
            if not intersect:
                Unique_deep_paths.append(deep_path)


    #Second group shallow nodes that have common deep path intersections
    Change = True
    while Change:
        i = 0
        Change = False
        while i<len(ClustersShallow)-1:
            paths = ClustersDeep[i]
            j = i+1
            while j<len(ClustersShallow):
                next_paths = ClustersDeep[j]
                Match = False
                for p in paths:
                    for np in next_paths:
                        if p==np:
                            Match = True
                if Match:
                    ClustersShallow[i].extend(ClustersShallow[j])
                    for np in next_paths:
                        if not np in ClustersDeep[i]:
                            ClustersDeep[i].append(np)
                    Change = True
                    del ClustersShallow[j]
                    del ClustersDeep[j]
                else:
                    j=j+1    
            i=i+1              

    
    #Construct a single dictionary of all clusters to return
    LocalClusters = {}
   
    for i in range(len(ClustersShallow)):
        shallow_paths = ClustersShallow[i]
        deep_paths = ClustersDeep[i]
        Cluster = {'Deep':{},'Shallow':{}} 
        for p in shallow_paths:            
            kmer = p[0][0].split('_')[2]
            if kmer in Cluster['Shallow']:
                Cluster['Shallow'][kmer].append(p)
            else:
                Cluster['Shallow'][kmer]=[p]
        for p in deep_paths:
            kmer = p[0][0].split('_')[2]
            if kmer in Cluster['Deep']:
                Cluster['Deep'][kmer].append(p)
            else:
                Cluster['Deep'][kmer]=[p]

        LocalClusters[i+1]=Cluster

    #print Unique_deep_paths
    for p in Unique_deep_paths:
        i=i+1
        Cluster = {'Deep':{},'Shallow':{}}
        #print p
        kmer = p[0][0].split('_')[2]
        if kmer in Cluster['Deep']:
            Cluster['Deep'][kmer].append(p)
        else:
            Cluster['Deep'][kmer]=[p]    
        LocalClusters[i+1]=Cluster

    return LocalClusters





     

def split_at_current_level(LOOM_Level,G,C,LC):
 
    simple_paths = getSimplePaths(LOOM_Level)
    

    #extract all paths from the simple_paths (dict keyed by k-mer)
    all_simple_paths = []
    
    for kmer in simple_paths:
        for path in simple_paths[kmer]:
            all_simple_paths.append(path)
    #all_simple_paths.sort(key=lambda x: LOOM_Level.node[x[0][0]]['location'])
    all_simple_paths.sort(key=lambda x: int(x[0][0].split('_')[1]))
    total_paths = len(all_simple_paths)
    overlapped = []
    i=0
    
    while i<total_paths:
        path1 = all_simple_paths[i]
        cluster=[path1]
        while i+1<total_paths and hasOverlap(path1,all_simple_paths[i+1],LOOM_Level):
            cluster.append(all_simple_paths[i+1])
            i+=1    
            path1 = all_simple_paths[i]
        overlapped.append(cluster)            
        i+=1

    
    graph_layers = list(dict.fromkeys([LOOM_Level.node[x]['layer'] for x in LOOM_Level.nodes()]))
    graph_layers.sort()
    split_indices = {l:[[1],[]] for l in graph_layers}
    split_indices_CN = {l:[[1],[]] for l in graph_layers}
    number_of_subgraphs = len(overlapped)+1

    
    for cluster in overlapped:
        first_path = cluster[0]
        last_path = cluster[-1]
        #add first layer split indices
        edge = first_path[0]
        end = LOOM_Level.node[edge[0]]['location']-1
        layer = LOOM_Level.node[edge[0]]['layer']
        split_indices[layer][1].append(end)
        split_indices_CN[layer][0].append(end+1)
        
        for edge in first_path:
            end = LOOM_Level.node[edge[1]]['location']-1
            layer = LOOM_Level.node[edge[1]]['layer']
            split_indices[layer][1].append(end)
            split_indices_CN[layer][0].append(end+1)

        edge = last_path[0]
        #start = LOOM_Level.node[edge[0]]['location']+1
        start = LOOM_Level.node[edge[0]]['location']
        layer = LOOM_Level.node[edge[0]]['layer']
        split_indices[layer][0].append(start)
            
        for edge in last_path:  
            #start = LOOM_Level.node[edge[1]]['location']+1
            start = LOOM_Level.node[edge[1]]['location']
            layer = LOOM_Level.node[edge[1]]['layer']
            split_indices[layer][0].append(start)

    #print split_indices_CN 
    #print split_indices
    #print overlapped


    #Split graph
    
    #Construct Sub-graphs based on boundaries
    #list of all segmented graphs
    subgraph_array =[]
    #print split_indices
    for i in range(number_of_subgraphs-1):
        S=nx.DiGraph()
        cp = nx.DiGraph()
        lcp = nx.DiGraph()
        for layer in split_indices:
            s = split_indices[layer][0][i]
            e = split_indices[layer][1][i]
            sc = split_indices_CN[layer][0][i]
            S.add_nodes_from([x for x in G.nodes() if G.node[x]['layer']==layer and G.node[x]['location']>=s and G.node[x]['location']<=e])
            #S.add_nodes_from([x for x in C.nodes() if C.node[x]['layer']==layer and C.node[x]['location']>=sc and C.node[x]['location']<=e])
            #S.add_nodes_from([x for x in LC.nodes() if LC.node[x]['layer']==layer and LC.node[x]['location']>=sc and LC.node[x]['location']<=e])
            cp.add_nodes_from([x for x in C.nodes() if C.node[x]['layer']==layer and C.node[x]['location']>=sc and C.node[x]['location']<=e])
            lcp.add_nodes_from([x for x in LC.nodes() if LC.node[x]['layer']==layer and LC.node[x]['location']>=sc and LC.node[x]['location']<=e])
        S.add_edges_from([x for x in G.edges() if S.has_node(x[0]) and S.has_node(x[1])])
        cp.add_edges_from([x for x in C.edges() if cp.has_node(x[0]) and cp.has_node(x[1])])
        lcp.add_edges_from([x for x in LC.edges() if lcp.has_node(x[0]) and lcp.has_node(x[1])])

        for x in cp.nodes():
            for attribute in C.node[x].keys():
                cp.node[x][attribute]=C.node[x][attribute]
 
        #add edges that continue past the current level
        last_cp = [x for x in cp.nodes() if cp.node[x]['layer']==layer]
        while last_cp:
            for x in last_cp:
                successors = list(C.successors(x))
                for s in successors:
                    cp.add_node(s,layer=C.node[s]['layer'],location=C.node[s]['location'],seq=C.node[s]['seq'],end=C.node[s]['end'],posx=C.node[s]['posx'])
                    cp.add_edge(x,s)
            layer+=1          
            last_cp = [x for x in cp.nodes() if cp.node[x]['layer']==layer]


        for x in lcp.nodes():
            for attribute in LC.node[x].keys():
                lcp.node[x][attribute]=LC.node[x][attribute]
 
        #add edges that continue past the current level
        last_cp = [x for x in lcp.nodes() if lcp.node[x]['layer']==layer]
        while last_cp:
            for x in last_cp:
                successors = list(LC.successors(x))
                for s in successors:
                    lcp.add_node(s,layer=LC.node[s]['layer'],location=LC.node[s]['location'],seq=LC.node[s]['seq'],end=LC.node[s]['end'],posx=LC.node[s]['posx'])
                    lcp.add_edge(x,s)
            layer+=1          
            last_cp = [x for x in lcp.nodes() if lcp.node[x]['layer']==layer]

        #cp=getNonIntersectingPaths(cp)
        if len(S.nodes())>0 or len(cp.nodes())>0 or len(lcp.nodes())>0:
            subgraph_array.append([S,cp,lcp])

    
    #print split_indices            
    i+=1
    S=nx.DiGraph()
    cp = nx.DiGraph()
    lcp = nx.DiGraph()
    for layer in split_indices:
        s = split_indices[layer][0][i]
        sc = split_indices_CN[layer][0][i]
        S.add_nodes_from([x for x in G.nodes() if G.node[x]['layer']==layer and G.node[x]['location']>=s])
        #S.add_nodes_from([x for x in C.nodes() if C.node[x]['layer']==layer and C.node[x]['location']>=sc])
        #S.add_nodes_from([x for x in LC.nodes() if LC.node[x]['layer']==layer and LC.node[x]['location']>=sc])
        cp.add_nodes_from([x for x in C.nodes() if C.node[x]['layer']==layer and C.node[x]['location']>=sc])
        lcp.add_nodes_from([x for x in LC.nodes() if LC.node[x]['layer']==layer and LC.node[x]['location']>=sc])
    S.add_edges_from([x for x in G.edges() if S.has_node(x[0]) and S.has_node(x[1])])
    cp.add_edges_from([x for x in C.edges() if cp.has_node(x[0]) and cp.has_node(x[1])])
    lcp.add_edges_from([x for x in LC.edges() if lcp.has_node(x[0]) and lcp.has_node(x[1])])
    
    
    for x in cp.nodes():
        for attribute in C.node[x].keys():
            cp.node[x][attribute]=C.node[x][attribute] 
    #add edges that continue past the current level
    last_cp = [x for x in cp.nodes() if cp.node[x]['layer']==layer]
    while last_cp:
        for x in last_cp:
            successors = list(C.successors(x))
            for s in successors:
                cp.add_node(s,layer=C.node[s]['layer'],location=C.node[s]['location'],seq=C.node[s]['seq'],end=C.node[s]['end'],posx=C.node[s]['posx'])
                cp.add_edge(x,s)
        layer+=1          
        last_cp = [x for x in cp.nodes() if cp.node[x]['layer']==layer]

    for x in lcp.nodes():
        for attribute in LC.node[x].keys():
            lcp.node[x][attribute]=LC.node[x][attribute] 
    #add edges that continue past the current level
    last_cp = [x for x in lcp.nodes() if lcp.node[x]['layer']==layer]
    while last_cp:
        for x in last_cp:
            successors = list(LC.successors(x))
            for s in successors:
                lcp.add_node(s,layer=LC.node[s]['layer'],location=LC.node[s]['location'],seq=LC.node[s]['seq'],end=LC.node[s]['end'],posx=LC.node[s]['posx'])
                lcp.add_edge(x,s)
        layer+=1          
        last_cp = [x for x in lcp.nodes() if lcp.node[x]['layer']==layer]
    
    
    #cp=getNonIntersectingPaths(cp)
    if len(S.nodes())>0 or len(cp.nodes())>0 or len(lcp.nodes())>0:
         subgraph_array.append([S,cp,lcp])

    #add attributes
    for sub in subgraph_array:
        S=sub[0]

        for x in [n for n in S.nodes() if n in G.nodes()]:
            for attribute in G.node[x].keys():
                S.node[x][attribute]=G.node[x][attribute]
        '''for x in [n for n in S.nodes() if n in LC.nodes()]:
            for attribute in LC.node[x].keys():
                S.node[x][attribute]=LC.node[x][attribute]
        for x in [n for n in S.nodes() if n in C.nodes()]:
            for attribute in C.node[x].keys():
                S.node[x][attribute]=C.node[x][attribute]'''    

    return subgraph_array




    
################################
def shrinkGraph(F,OptimisedComplex,solver,level,maxedges):
    print "Shrinking"
    optimised_complex_paths = getNonIntersectingPaths(OptimisedComplex)
    all_paths_F = getNonIntersectingPaths(F)
    all_paths_list = []
    for kmer in all_paths_F:
        for path in all_paths_F[kmer]:
            all_paths_list.append((kmer,path))
            all_paths_list.sort(key=lambda x: len(x[1]))
    a = 0
    while len(F.edges())>maxedges and a<len(all_paths_list):
        special_case = nx.DiGraph()
        Clusters = getLocalIntersectingPaths(optimised_complex_paths,{kmer:[all_paths_list[a][1]]})

        for cl in Clusters:
            deep_paths = Clusters[cl]['Deep']
            shallow_paths = Clusters[cl]['Shallow']
            if not deep_paths:
                special_case = add_paths_to_graph(special_case,F,shallow_paths)
            if not shallow_paths:
                special_case = add_paths_to_graph(special_case,OptimisedComplex,deep_paths)            

            elif deep_paths and shallow_paths:
                CG = nx.DiGraph()
                CG = add_paths_to_graph(CG,F,shallow_paths)
                
                CGComplex = nx.DiGraph()
                CGComplex = add_paths_to_graph(CGComplex,OptimisedComplex,deep_paths)
                CG = solve_lpp_optimise(CG,CGComplex,solver,level,maxedges)

                #Restore previous edges between allowed repeated kmers of F that have now been selected
                '''for e in F.edges():
                    if e[0] in CG.nodes() and e[1] in CG.nodes():
                        CG.add_edge(e[0],e[1])'''
                        
                special_case.add_nodes_from(CG.nodes(data=True))
                special_case.add_edges_from(CG.edges())

        OptimisedComplex = special_case
        optimised_complex_paths = getNonIntersectingPaths(OptimisedComplex)

        to_remove = []
        for edge in all_paths_list[a][1]:
            to_remove.extend([edge[0],edge[1]])
            F.remove_nodes_from(to_remove)
        a+=1
    return F,OptimisedComplex,optimised_complex_paths
        

def imposeComplexConstraints(G,C,level):

    #impose complex path constraints at current depth
    G=remove_intersecting_edges(G,C)
    G = prune_graph_level(G,level)

    return G







def solve_lpp_optimise(G,C,solver,level,maxedges):
    
    F=nx.DiGraph()
    #if graph size exceeds limit after stricter constraints, return complex paths 
    if len(G.edges())+len(C.edges())>maxedges:
        print "limit 1 imposed"
        G = imposeComplexConstraints(G,C,level)


    if len(G.edges())+len(C.edges())>maxedges:
        print "limit 2 imposed"
        F.add_nodes_from(C.nodes(data=True))
        F.add_edges_from(C.edges())

    else:

        #if graph only has one kmer do not apply any constraints
        G.add_nodes_from(C.nodes(data=True))
        G.add_edges_from(C.edges())
        graph_kmers = list(dict.fromkeys([G.node[x]['seq'] for x in G.nodes()]))

        
        if len(graph_kmers)==1:

            graph_nodes = list(G.nodes())
            graph_egdes=list(G.edges())

            #Get first and last layers of graph:
            graph_layers = list(dict.fromkeys([G.node[x]['layer'] for x in G.nodes()]))
            graph_layers.sort()

            #Define and solve the lp problem
            #Define variables (edges are variables, where each edge can = 0 or 1 only)
            vars = LpVariable.dicts("edge",[i for i in G.edges()],0,1,LpInteger)

            #Define a maximisation problem
            prob = LpProblem("problem", LpMaximize)

            #Define the objectiove: maximization of the sum of the edges
            prob += lpSum(vars), "Sum_of_Edges" #ADDS THE OBJECTIVE ATTRIBUTE TO PROB
    
    
            for i in graph_nodes:

                if int(G.node[i]['layer'])<graph_layers[level-1] and int(G.node[i]['layer'])!=graph_layers[0]:
                    prob += lpSum([vars[(i,a)] for a in G.successors(i)])<= 100*lpSum([vars[(b,i)] for b in G.predecessors(i)])          
                    prob += lpSum([vars[(a,i)] for a in G.predecessors(i)])<= 100*lpSum([vars[(i,b)] for b in G.successors(i)])


                #if the node is not in the last layer
                if float(G.node[i]['layer'])<graph_layers[level-1]:

                    #iterate through the successors of that node:            
                    for r in G.successors(i):

                        for x in [e for e in graph_egdes if G.node[i]['layer']==G.node[e[0]]['layer'] if G.node[i]['end']>=G.node[e[0]]['location'] and G.node[i]['location']<G.node[e[0]]['location']]:
                            if G.node[i]['end']-G.node[x[0]]['location']!=G.node[r]['end']-G.node[x[1]]['location']:
                                prob += lpSum(vars[m] for m in [(i,r),x])<=1


                        for x in [e for e in graph_egdes if G.node[i]['layer']==G.node[e[0]]['layer'] if G.node[r]['end']>=G.node[e[1]]['location'] and G.node[r]['location']<G.node[e[1]]['location']]:
                            if G.node[i]['end']-G.node[x[0]]['location']!=G.node[r]['end']-G.node[x[1]]['location']:
                                prob += lpSum(vars[m] for m in [(i,r),x])<=1


            #solve the lp-problem with the lp-solver of choice
            if solver=='GUROBI':
                prob.solve(pulp.GUROBI(msg=0))
           

            elif solver=='CBC':
                prob.solve()
  

            #Construct Solution Graph
            for v in prob.variables():
                #if the edge got a score of 1 (meaning, it is within the solution):
                if round(v.varValue)==1.0:
 
                    #Add both nodes and the edge to the solution graph
                    edge=(v.name.split('edge_')[1]).strip('(').strip(')').split(',_')
                    node1 = edge[0].strip("'")
                    node2 = edge[1].strip("'")
                    F.add_node(node1)
                    F.add_node(node2)
                    F.add_edge(node1,node2)

                    for attribute in G.node[node1].keys():
                        F.node[node1][attribute]=G.node[node1][attribute]
                        F.node[node2][attribute]=G.node[node2][attribute]


        else:

            graph_nodes = list(G.nodes())
            graph_egdes=list(G.edges())
    

            #Get first and last layers of graph:
            graph_layers = list(dict.fromkeys([G.node[x]['layer'] for x in G.nodes()]))
            graph_layers.sort()

            #Define and solve the lp problem
            #Define variables (edges are variables, where each edge can = 0 or 1 only)
            vars = LpVariable.dicts("edge",[i for i in G.edges()],0,1,LpInteger)

            #Define a maximisation problem
            prob = LpProblem("problem", LpMaximize)

            #Define the objectiove: maximization of the sum of the edges
            prob += lpSum(vars), "Sum_of_Edges" #ADDS THE OBJECTIVE ATTRIBUTE TO PROB

            #Define the constraints
            special_kmers = list(dict.fromkeys([C.node[x]['seq'] for x in C.nodes() if C.node[x]['layer']==1]))
    
            for kmer in special_kmers:
                special_edges = [x for x in G.edges() if G.node[x[0]]['seq']== kmer and G.node[x[0]]['layer']== 1]
                if special_edges:
                    prob += lpSum([vars[e] for e in special_edges])>=1

            for i in graph_nodes:
                if int(G.node[i]['layer'])<graph_layers[level-1] and int(G.node[i]['layer'])!=graph_layers[0]:
                    prob += lpSum([vars[(i,a)] for a in G.successors(i)])<= 100*lpSum([vars[(b,i)] for b in G.predecessors(i)])          
                    prob += lpSum([vars[(a,i)] for a in G.predecessors(i)])<= 100*lpSum([vars[(i,b)] for b in G.successors(i)])



                #if the node is not in the last layer
                if float(G.node[i]['layer'])<graph_layers[level-1]:

                    #iterate through the successors of that node:            
                    for r in G.successors(i):
              
                        #iterate through the edges of other nodes in the same layer:
                        all_edges_layer = [e for e in graph_egdes if G.node[e[0]]['layer']==G.node[i]['layer'] if e!=(i,r)]
                        constraint_set = [e for e in all_edges_layer if G.node[i]['location']<G.node[e[0]]['location'] and G.node[r]['location']>G.node[e[1]]['location']]
                        constraint_kmers = [G.node[e[0]]['seq'] for e in constraint_set]
                        constraint_kmers.append(G.node[i]['seq'])
                        constraint_kmers = list(dict.fromkeys(constraint_kmers))
                
                        if len(constraint_kmers)>1:
                            for x in constraint_set:               
                                prob += lpSum(vars[m] for m in [(i,r),x])<=1
                
                        constraint_set = [e for e in all_edges_layer if G.node[i]['location']>G.node[e[0]]['location'] and G.node[r]['location']<G.node[e[1]]['location']]    
                        constraint_kmers = [G.node[e[0]]['seq'] for e in constraint_set]
                        constraint_kmers.append(G.node[i]['seq'])
                        constraint_kmers = list(dict.fromkeys(constraint_kmers))
                
                        if len(constraint_kmers)>1:
                            for x in constraint_set:               
                                prob += lpSum(vars[m] for m in [(i,r),x])<=1


                        for x in [e for e in graph_egdes if G.node[i]['layer']==G.node[e[0]]['layer'] if G.node[i]['end']>=G.node[e[0]]['location'] and G.node[i]['location']<G.node[e[0]]['location']]:
                            if G.node[i]['end']-G.node[x[0]]['location']!=G.node[r]['end']-G.node[x[1]]['location']:
                                prob += lpSum(vars[m] for m in [(i,r),x])<=1


                        for x in [e for e in graph_egdes if G.node[i]['layer']==G.node[e[0]]['layer'] if G.node[r]['end']>=G.node[e[1]]['location'] and G.node[r]['location']<G.node[e[1]]['location']]:
                            if G.node[i]['end']-G.node[x[0]]['location']!=G.node[r]['end']-G.node[x[1]]['location']:
                                prob += lpSum(vars[m] for m in [(i,r),x])<=1


            #solve the lp-problem with the lp-solver of choice
            if solver=='GUROBI':
                prob.solve(pulp.GUROBI(msg=0))
           

            elif solver=='CBC':
                prob.solve()
  

            #Construct Solution Graph
            for v in prob.variables():
                #if the edge got a score of 1 (meaning, it is within the solution):
                if round(v.varValue)==1.0:
 
                    #Add both nodes and the edge to the solution graph
                    edge=(v.name.split('edge_')[1]).strip('(').strip(')').split(',_')
                    node1 = edge[0].strip("'")
                    node2 = edge[1].strip("'")
                    F.add_node(node1)
                    F.add_node(node2)
                    F.add_edge(node1,node2)

                    for attribute in G.node[node1].keys():
                        F.node[node1][attribute]=G.node[node1][attribute]
                        F.node[node2][attribute]=G.node[node2][attribute]

    return F




def get_shallow_overlaps(graphS,graphD):
    #Returns a list of nodes from shallower graph, with longer kmers, where kmers start at same base as shorter kmers with more depth
    #In this case the shallower path must be removed and longer paths are favours (for now this is done for complex paths only
    #Shallow simple paths retained for efficiency

    F=nx.DiGraph()
    F.add_nodes_from(graphS.nodes(data=True))
    F.add_edges_from(graphS.edges())


    overlapping_nodes = []
    for n in graphD.nodes():
        overlap = [x for x in graphS.nodes() if graphS.node[x]['layer']==graphD.node[n]['layer'] and graphS.node[x]['location']<=graphD.node[n]['location'] and graphD.node[n]['location']<=graphS.node[x]['end']]
        overlapping_nodes.extend(overlap)
        overlap = [x for x in graphS.nodes() if graphS.node[x]['layer']==graphD.node[n]['layer'] and graphS.node[x]['location']>=graphD.node[n]['location'] and graphS.node[x]['location']<=graphD.node[n]['end']]
        overlapping_nodes.extend(overlap)


    #get all neigbours of overlapping kmers
    neighbors = [] 
    UF = F.to_undirected()
    for s in overlapping_nodes:
        neighbors.extend(nx.single_source_shortest_path(UF,s).keys())

    overlapping_nodes.extend(neighbors)
    overlapping_nodes = list(dict.fromkeys(overlapping_nodes))
    return overlapping_nodes



def recursive_segmentation(sub_segments,kmers_len,solver,LevelGraphs,prune,level,min_depth,maxedges,noconstraints):
    FinalGraph = nx.DiGraph()
    if level>=2:
        for segment in sub_segments:

            segments_of_segment = [segment]

            #Step 1: perform initial splitting of graph based on longer kmers (the main graph and complex_paths are split into independent segments)
            #Determine if additional splitting is required
            #graph = segment[0]
            #Construct graph for pruning and solving

            #Z= nx.DiGraph()
            #Z.add_nodes_from(graph.nodes(data=True))
            #Z.add_edges_from(graph.edges())
            #Z = prune_graph_level(Z,level)
          
                
            #if level in LevelGraphs and noconstraints==False and len(segment[0].nodes())>0:
                #segments_of_segment = split_at_current_level(LevelGraphs[level],segment[0],segment[1],segment[2])


            
            #Step 2: run framework on each subgraph of the current subgraph
            for seg_of_seg in segments_of_segment:
                #Graph of nodes of length = current kmer length            
                graph = seg_of_seg[0]
                #Construct graph for pruning and solving
                F=nx.DiGraph()
                F.add_nodes_from(graph.nodes(data=True))
                F.add_edges_from(graph.edges())

                #Graph of complex paths from solved at current kmer length - may contain unconnected nodes due to splitting
                C=seg_of_seg[1]

                #Graph of complex paths solved previously at longer kmer length (note this graph contains paths shallower than current depth)
                LC = seg_of_seg[2]
                #Correct and prune complex paths that do not have ful connection at current level
                C = prune_graph_level_plus(C,level)

                F.remove_nodes_from(C.nodes())
                nodes_to_remove = cut_and_prune(kmers_len,C)
                F.remove_nodes_from(nodes_to_remove)

                PrunedLC = prune_graph_level_plus(LC,level) #gets longer complex paths at current depth
                F.remove_nodes_from(PrunedLC.nodes())
                nodes_to_remove = cut_and_prune(kmers_len,PrunedLC)
                F.remove_nodes_from(nodes_to_remove)

                #if level==9 and kmers_len==8:
                    #print level,C.nodes(),LC.nodes(),PrunedLC.nodes()

                #Check for direct overlap (start==start of kmers)
                to_remove=get_shallow_overlaps(PrunedLC,C)
                PrunedLC.remove_nodes_from(to_remove)
                LC.remove_nodes_from(to_remove)


                #Remove nodes that belong to longer simple kmers found previously at the current level of layers
                if level in LevelGraphs and noconstraints==False:
                    nodes_to_remove = cut_and_prune(kmers_len,LevelGraphs[level])
                    F.remove_nodes_from(nodes_to_remove)
                

                #Calculate an optimised graph of complex_paths: Favour shorter kmers with greater depth
                OptimisedComplex = nx.DiGraph()
                if len(C.nodes())>0 and len(PrunedLC.nodes())>0:
                    #if total number of edges is large: split into clusters
                    if len(C.edges())+len(PrunedLC.edges())>maxedges:
                        longer_complex_paths = getNonIntersectingPaths(PrunedLC)
                        complex_paths = getNonIntersectingPaths(C)
                        Clusters = getLocalIntersectingPaths(complex_paths,longer_complex_paths)
                        for cl in Clusters:
                            deep_paths = Clusters[cl]['Deep']
                            shallow_paths = Clusters[cl]['Shallow']

                            if not deep_paths:
                                OptimisedComplex = add_paths_to_graph(OptimisedComplex,PrunedLC,shallow_paths)

                            elif not shallow_paths:
                                OptimisedComplex = add_paths_to_graph(OptimisedComplex,C,deep_paths)
                            else:
                                CG = nx.DiGraph()
                                CG = add_paths_to_graph(CG,PrunedLC,shallow_paths)

                                CGComplex = nx.DiGraph()
                                CGComplex = add_paths_to_graph(CGComplex,C,deep_paths)

                                                              
                                CG = solve_lpp_optimise(CG,CGComplex,solver,level,maxedges)
             
                                OptimisedComplex.add_nodes_from(CG.nodes(data=True))
                                OptimisedComplex.add_edges_from(CG.edges())
            
                    else: #if total graph is smaller enough, simply compute in single step
                        OptimisedComplex = solve_lpp_optimise(PrunedLC,C,solver,level,maxedges)
                        #OptimisedComplex = prune_graph_level_plus(OptimisedComplex,level)
                        #ed = [e for e in OptimisedComplex.edges() if OptimisedComplex.node[e[0]]['seq']=='TGTACATT']

                elif len(C.nodes())>0 and len(PrunedLC.nodes())==0:
                    OptimisedComplex = C
                elif len(C.nodes())==0 and len(PrunedLC.nodes())>0:
                    OptimisedComplex = PrunedLC



                F = prune_graph_level(F,level)


                if len(F.nodes())>0: #Proceed with solving at current level
                    optimised_complex_paths = getNonIntersectingPaths(OptimisedComplex)

                    if len(F.edges())>maxedges:
                        F,OptimisedComplex,optimised_complex_paths = shrinkGraph(F,OptimisedComplex,solver,level,maxedges)
                    OptimisedGraph = nx.DiGraph() 
                    if len(F.nodes())>0:
                        if len(F.edges())+len(OptimisedComplex.edges())>maxedges:
                            F = solve_lpp(F,solver)
                            if len(F.nodes())>0:
                                if len(F.edges())+len(OptimisedComplex.edges())>maxedges:
                                    
                                    all_paths_F = getNonIntersectingPaths(F)
                                    Clusters = getLocalIntersectingPaths(optimised_complex_paths,all_paths_F)
                                    for cl in Clusters:
                                        deep_paths = Clusters[cl]['Deep']
                                        shallow_paths = Clusters[cl]['Shallow']
                                        if not deep_paths:
                                            OptimisedGraph = add_paths_to_graph(OptimisedGraph,F,shallow_paths)

                                        elif not shallow_paths:
                                            OptimisedGraph = add_paths_to_graph(OptimisedGraph,OptimisedComplex,deep_paths)

                                        else:
                                            CG = nx.DiGraph()
                                            CG = add_paths_to_graph(CG,F,shallow_paths)
                                            
                                            CGComplex = nx.DiGraph()
                                            CGComplex = add_paths_to_graph(CGComplex,OptimisedComplex,deep_paths)
                                            
                                            CG = solve_lpp_optimise(CG,CGComplex,solver,level,maxedges)
                                        
                                            OptimisedGraph.add_nodes_from(CG.nodes(data=True))
                                            OptimisedGraph.add_edges_from(CG.edges())
          
                                else:
                                    OptimisedGraph = solve_lpp_optimise(F,OptimisedComplex,solver,level,maxedges)

                            
                        else:
                            OptimisedGraph = solve_lpp_optimise(F,OptimisedComplex,solver,level,maxedges)
                    else:
                        OptimisedGraph.add_nodes_from(OptimisedComplex.nodes(data=True))
                        OptimisedGraph.add_edges_from(OptimisedComplex.edges())

                    

                    #get new simple and complex paths of optimised graph
                    if len(OptimisedGraph.nodes())>0:
                           
                        simple_paths = getSimplePaths(OptimisedGraph)
                        complex_paths = getComplexPaths(OptimisedGraph)
                        #Remove PrunedLC.nodes from LC, complex_paths now contains the optimised Longer Complex paths    
                        to_remove = [x for x in PrunedLC.nodes()]
                        LC.remove_nodes_from(to_remove)

                        #Add optimised simple paths to graph (longer complex paths may have become simple)
                        OptimisedSimple = nx.DiGraph()
                        OptimisedSimple = add_paths_to_graph(OptimisedSimple,OptimisedGraph,simple_paths)
                        #if level==9 and kmers_len==8:
                            #print complex_paths
                            #print LC.nodes()
                        LC=remove_intersecting_edges(LC,OptimisedSimple)

                        #OptimisedComplex = nx.DiGraph()
                        #OptimisedComplex = add_paths_to_graph(OptimisedComplex,OptimisedGraph,complex_paths)                       
                        #LC=remove_intersecting_edges(LC,OptimisedComplex)

                        CN = nx.DiGraph()
                        CN = add_paths_to_graph(CN,OptimisedGraph,complex_paths) 


                        #add simple paths to final graph
                        FinalGraph = add_paths_to_graph(FinalGraph,OptimisedGraph,simple_paths)

                        last_layer = [x for x in OptimisedSimple.nodes() if OptimisedSimple.node[x]['layer']>level]
                        OptimisedSimple.remove_nodes_from(last_layer)
                        simple_paths = getSimplePaths(OptimisedSimple)

                        #add simple paths to graph incase a longer kmers have now become a simple path
                        graph = add_paths_to_graph(graph,OptimisedSimple,simple_paths)
                        #eliminate last layernd layers greater than current level a and split graph based on simple paths
                        last_layer = [x for x in graph.nodes() if graph.node[x]['layer']>=level]
                        graph.remove_nodes_from(last_layer)
 
    
                        if simple_paths:
                            graph_segments = split_graph(graph,simple_paths,CN,LC)
                            #print "Simple path found, segments = "+str(len(graph_segments))
                        else:
                            graph_segments = [[graph,CN,LC]]
                            #print "No Simple path found, segments = "+str(len(graph_segments))                            


                    else:
                                                   
                        last_layer = [x for x in graph.nodes() if graph.node[x]['layer']==level]
                        graph.remove_nodes_from(last_layer)
                        #Remove PrunedLC.nodes from LC, complex_paths now contains the optimised Longer Complex paths    
                        to_remove = [x for x in PrunedLC.nodes()]
                        LC.remove_nodes_from(to_remove)
                        graph_segments = [[graph,OptimisedComplex,LC]]
                        
                else:#pruned graph is empty
                    last_layer = [x for x in graph.nodes() if graph.node[x]['layer']==level]
                    #Remove PrunedLC.nodes from LC, complex_paths now contains the optimised Longer Complex paths
                    to_remove = [x for x in PrunedLC.nodes()]
                    LC.remove_nodes_from(to_remove)
                    graph_segments = [[graph,OptimisedComplex,LC]]
                    

            
                del graph           
                RF = recursive_segmentation(graph_segments,kmers_len,solver,LevelGraphs,prune,level-1,min_depth,maxedges,noconstraints)
                FinalGraph.add_nodes_from(RF.nodes(data=True)) 
                FinalGraph.add_edges_from(RF.edges())

    else: #level is less than min_depth

        # add complex paths for each sub_segment of graph: note that all simple paths have already beed added at this point

        for segment in sub_segments:
            C = segment[1]         
            C = prune_graph_level_plus(C,2)

            FinalGraph.add_nodes_from(C.nodes(data=True))
            FinalGraph.add_edges_from(C.edges())

            #FinalGraph.add_nodes_from(LC.nodes(data=True))
            #FinalGraph.add_edges_from(LC.edges())

            
    return FinalGraph
####################################


def remove_intersecting_edges(F,LOOM):

    graph_layers = list(dict.fromkeys([F.node[x]['layer'] for x in F.nodes()]))
    graph_layers.extend(list(dict.fromkeys([LOOM.node[x]['layer'] for x in LOOM.nodes()])))
    graph_layers = list(dict.fromkeys(graph_layers))
    dict_of_edges = {l:{} for l in graph_layers}
    for node in F.nodes():
        node_successors = F.successors(node)
        for ns in node_successors:
            dict_of_edges[F.node[node]['layer']][(node,ns)]={}
            dict_of_edges[F.node[node]['layer']][(node,ns)]['location'] = F.node[node]['location']
            dict_of_edges[F.node[node]['layer']][(node,ns)]['slocation'] = F.node[ns]['location']
            dict_of_edges[F.node[node]['layer']][(node,ns)]['end'] = F.node[node]['end']
            dict_of_edges[F.node[node]['layer']][(node,ns)]['send'] = F.node[ns]['end']

            #dict_of_edges[F.node[node]['layer']][F.node[node]['location']]['edge']=(node,ns)

    remove_edges = []           
    constraint_edges = LOOM.edges()
    for ce in constraint_edges:
        remove_edges = []
        layer = LOOM.node[ce[0]]['layer']
        location1 = LOOM.node[ce[0]]['location']
        location2 = LOOM.node[ce[1]]['location']
        end1 = LOOM.node[ce[0]]['end']
        end2 = LOOM.node[ce[1]]['end']
        remove_edges.extend([e for e in dict_of_edges[layer] if dict_of_edges[layer][e]['location']>location1 and dict_of_edges[layer][e]['slocation']<location2])
        remove_edges.extend([e for e in dict_of_edges[layer] if dict_of_edges[layer][e]['location']<location1 and dict_of_edges[layer][e]['slocation']>location2])
        remove_edges = list(dict.fromkeys(remove_edges))
        for e in remove_edges:
            del dict_of_edges[layer][e]
            F.remove_edge(e[0],e[1])

        overlaps = []
        
        overlaps.extend([e for e in dict_of_edges[layer] if dict_of_edges[layer][e]['location']<=end1 and dict_of_edges[layer][e]['location']>=location1])
        overlaps.extend([e for e in dict_of_edges[layer] if dict_of_edges[layer][e]['slocation']<=end2 and dict_of_edges[layer][e]['slocation']>=location2])
        overlaps.extend([e for e in dict_of_edges[layer] if dict_of_edges[layer][e]['send']>=location2 and dict_of_edges[layer][e]['slocation']<=location2])
        overlaps.extend([e for e in dict_of_edges[layer] if dict_of_edges[layer][e]['end']>=location1 and dict_of_edges[layer][e]['location']<=location1])
        
        overlaps = list(dict.fromkeys(overlaps)) 
        for e in overlaps:
            if end1-dict_of_edges[layer][e]['location']!=end2-dict_of_edges[layer][e]['slocation']:
                remove_edges.append(e)
                del dict_of_edges[layer][e]
                F.remove_edge(e[0],e[1])
                         
    return F



def cut_and_prune(kmers_len,LG):
    nodes_to_remove=[]
    all_nodes = list(LG.nodes())
    for x in all_nodes:
        sequence = LG.node[x]['seq']
        start = LG.node[x]['location']
        layer = LG.node[x]['layer']
        for i in range(len(sequence)):
            if i<=(len(sequence)-kmers_len):
                kmer = sequence[i:i+kmers_len]
                loc = start+i
                nodes_to_remove.append(str(layer)+'_'+str(loc)+'_'+kmer)

    return nodes_to_remove
                

def get_graphs_per_level(LOOM,min_depth):
    LevelGraphs = {}
    total_levels = sorted(list(dict.fromkeys([LOOM.node[x]['layer'] for x in LOOM.nodes()])))
    if total_levels:
        total_levels = total_levels[-1]
        for i in range(total_levels,min_depth-1,-1):
            F = prune_graph_level(LOOM,i)
            LevelGraphs[i]=F
    return LevelGraphs



def combine_kmers_in_layer(nodes):
    combined_nodes = []
    i = 0
    while i<len(nodes):
        start_location = nodes[i][0]
        combined_kmer = nodes[i][2]
        end_combined_kmer = nodes[i][1]
        while i+1<len(nodes) and nodes[i+1][0]<=end_combined_kmer:
            kmer2 = nodes[i+1][2]
            end_kmer2 = nodes[i+1][1]
            overhangIndex =end_kmer2-end_combined_kmer
            if overhangIndex>0: #account for case where a shorter k-mer is contained in longer k-mer  
                combined_kmer+=kmer2[-overhangIndex:]
                end_combined_kmer = end_kmer2
            i+=1
        combined_nodes.append((start_location,end_combined_kmer,combined_kmer))            
        i+=1    
    return combined_nodes




def combine_kmers_with_segments(nodes_dict):
    
    combined_nodes_dict={}
    for layer in nodes_dict:
        combined_nodes_of_layer = []
        nodes_of_layer = nodes_dict[layer]
        i = 0
        while i<(len(nodes_of_layer)):
            start_location = nodes_of_layer[i][0]
            combined_kmer = nodes_of_layer[i][2]
            combined_kmer_end = nodes_of_layer[i][1]
            segments = [nodes_of_layer[i]]
            while i+1<len(nodes_of_layer) and nodes_of_layer[i+1][0]<=combined_kmer_end:
                i+=1
                kmer2 = nodes_of_layer[i][2]
                end_kmer2 = nodes_of_layer[i][1]
                overhangIndex =end_kmer2-combined_kmer_end
                segments.append(nodes_of_layer[i])
                if overhangIndex>0: #account for case where a shorter k-mer is contained in longer k-mer           
                    combined_kmer+=kmer2[-overhangIndex:]
                    combined_kmer_end=end_kmer2
            #sort segments by position and length (for output of overlapping k-mers)
            #segments = sorted(segments,key=lambda x: (x[0],-len(x[2])))
            segments = sorted(segments,key=lambda x: x[3])
            combined_nodes_of_layer.append(((start_location,combined_kmer_end,combined_kmer),segments))
            i+=1
        combined_nodes_dict[layer] = combined_nodes_of_layer    
    
    return combined_nodes_dict

def get_kmers_per_layer(All_Graphs,number_of_layers):


    #Return a dictionary of combined k-mers per layer, where the internal segments of each combined k-mer are combined at each specific depth
    #LOOM_layers = list(dict.fromkeys([LOOM.node[x]['layer'] for x in LOOM.nodes()]))
    All_Combined_Kmers =  {i+1:[] for i in range(number_of_layers)}
    Checked_Combined_Kmers =  {i+1:[] for i in range(number_of_layers)}
    levels = All_Graphs.keys()
    levels.sort(reverse=True)

    for level in levels:
        G = All_Graphs[level]
        G_layers = list(dict.fromkeys([G.node[x]['layer'] for x in G.nodes()]))
        for layer in G_layers:
            nodes_in_layer = [(G.node[x]['location'],G.node[x]['end'],G.node[x]['seq']) for x in G.nodes() if G.node[x]['layer']==layer]
            nodes_in_layer = sorted(nodes_in_layer,key=lambda x: (x[0],-len(x[2])))   
            #Combine the k-mers of this specific layer
            nodes_in_layer = combine_kmers_in_layer(nodes_in_layer)
            for node in nodes_in_layer:
                if not node in Checked_Combined_Kmers[layer]:
                    Checked_Combined_Kmers[layer].append(node)
                    node = (node[0],node[1],node[2],level)
                    All_Combined_Kmers[layer].append(node)

        #print '***************'    
    #sort final list of combined nodes in each layer, first by position and then by length
    
    for layer in All_Combined_Kmers:
        All_Combined_Kmers[layer] = sorted(All_Combined_Kmers[layer],key=lambda x: (x[0],-len(x[2])))
    All_Combined_Kmers = combine_kmers_with_segments(All_Combined_Kmers)
    return All_Combined_Kmers


def calculate_modules(Graph_at_levels,number_of_layers):

    levels = Graph_at_levels.keys()
    levels.sort(reverse=True)

    #Modules = {level:{} for level in levels}
    Modules = {}
    reference = []
    for level in levels:
        mod_id = 1
        all_deep_levels = [x for x in levels if x>= level]
        current_graphs = {x:Graph_at_levels[x] for x in all_deep_levels}
        kmers_dict = get_kmers_per_layer(current_graphs,number_of_layers)
        #kmers_dict = get_kmers_per_layer({level:Graph_at_levels[level]},number_of_layers)
        temp = [x[0] for x in kmers_dict[1]]
        if temp!=reference:
            reference = [x[0] for x in kmers_dict[1]]
            #temp.sort(key = lambda x:x[0])
            #proceed with module calulation
            #if len(reference)<4: #do not split into groups of modules
            Modules_at_level={}
            pattern_list = [site[2] for site in reference]
            pattern = {site[2]:0 for site in reference}
            layers = [x for x in kmers_dict.keys() if len(kmers_dict[x])>0 and x<=level]
            layers.sort()
            #if layers:
                #Modules[level][mod_id] = {}
            mod_in_layers = {}
            for layer in layers:
                
                kmers_layer = [x[0] for x in kmers_dict[layer]]
                inner_kmers_layer = [x[1] for x in kmers_dict[layer]]
                count = {site[2]:0 for site in kmers_layer}
                for site in kmers_layer:
                    count[site[2]]+=1
                for site in count:
                    if count[site]>pattern[site]:
                        pattern[site] = count[site]

                if len(count)>1:
                    mod_in_layers[layer]=(kmers_layer,inner_kmers_layer)
            if mod_in_layers:
                Modules_at_level[mod_id] = mod_in_layers
                regular_exp = [pattern_list[0]]
                if pattern[pattern_list[0]]>1:
                    regular_exp.append('*(n)')
                for site in pattern_list[1:]:
                    regular_exp.append('-'+site)
                    if pattern[site]>1:
                        regular_exp.append('*(n)')
                Modules_at_level[mod_id]['RE'] = regular_exp
            if Modules_at_level:
                Modules[level]=Modules_at_level
    return Modules

def get_modules(Main_LOOM_Level, Prime5_Levels, Prime3_Levels,number_of_layers):

    
    mainModules = calculate_modules(Main_LOOM_Level,number_of_layers)
    prime5Modules = calculate_modules(Prime5_Levels,number_of_layers)
    prime3Modules = calculate_modules(Prime3_Levels,number_of_layers)

    all_levels = []
    all_levels.extend(mainModules.keys())
    all_levels.extend(prime5Modules.keys())
    all_levels.extend(prime3Modules.keys())

    all_levels = dict.fromkeys(all_levels)
    AllModules = {l:{'Main':{},'5prime':{},'3prime':{}} for l in all_levels}
    for level in mainModules:
        AllModules[level]['Main'] = mainModules[level]
    for level in prime5Modules:
        AllModules[level]['5prime'] = prime5Modules[level]
    for level in prime3Modules:
        AllModules[level]['3prime'] = prime3Modules[level]
   
    return AllModules



def kmers_by_conservation(topNodes,number_of_layers):

    All_Combined_Kmers =  {i+1:[] for i in range(number_of_layers)}
    for node in topNodes:
        layer = node[3]-1
        for i in range(layer,0,-1):
            All_Combined_Kmers[i].append(node)

    for layer in All_Combined_Kmers:
        All_Combined_Kmers[layer] = sorted(All_Combined_Kmers[layer],key=lambda x: (x[0],-len(x[2])))
    All_Combined_Kmers = combine_kmers_with_segments(All_Combined_Kmers)
   
    return All_Combined_Kmers



      
def add_paths_to_graph(graphA,graphB,paths):

    for kmer in paths:
        for p in paths[kmer]:
            for edge in p:
                for x in edge:
                    graphA.add_node(x)
                    for attribute in graphB.node[x].keys():
                        graphA.node[x][attribute] = graphB.node[x][attribute]
                graphA.add_edge(edge[0],edge[1])

    return graphA


def calculateMissing5Prime(F,total_seq,tolerance):
    layersToRemove = []
    firstNodeLocations = []
    firstNodesDict ={}
    for i in range(total_seq):
        layer = i+1
        positions_of_nodes_in_layer =[F.node[x]['location'] for x in F.nodes() if F.node[x]['layer']==layer]
        positions_of_nodes_in_layer.sort()
        if positions_of_nodes_in_layer: #If the layer contains any nodes     
            firstNode = positions_of_nodes_in_layer[0]
            firstNodeLocations.append(firstNode)
            firstNodesDict[layer] = firstNode
        else:
            layersToRemove.append(layer) #if the layer has no nodes, remove it for subsequent search 

    median_pos_node = 0
    firstNodeLocations.sort()
  
    num_nodes = len(firstNodeLocations)
    if num_nodes%2!=0:
        median = float(firstNodeLocations[num_nodes/2])
    else:
        median = (firstNodeLocations[num_nodes/2]+firstNodeLocations[(num_nodes/2)-1])/2.0

    for layer in firstNodesDict:
        pos = firstNodesDict[layer]
        if pos<=median*tolerance[0]:
            layersToRemove.append(layer)     
    return layersToRemove



def calculateMissing3Prime(F,total_seq,tolerance):
    layersToRemove = []
    lastNodeLocations = []
    lastNodesDict ={}
    for i in range(total_seq):
        layer = i+1
        positions_of_nodes_in_layer =[F.node[x]['posx'] for x in F.nodes() if F.node[x]['layer']==layer]       
        positions_of_nodes_in_layer.sort()

        if positions_of_nodes_in_layer: #If the layer contains any nodes     
            lastNode = positions_of_nodes_in_layer[-1]
            lastNodeLocations.append(lastNode)
            lastNodesDict[layer] = lastNode
        else:
            layersToRemove.append(layer) #if the layer has no nodes, remove it for subsequent search 
    
    median_pos_node = 0
    lastNodeLocations.sort()
    num_nodes = len(lastNodeLocations)
    if num_nodes%2!=0:
        median = float(lastNodeLocations[num_nodes/2])
    else:
        median = (lastNodeLocations[num_nodes/2]+lastNodeLocations[(num_nodes/2)-1])/2.0

    for layer in lastNodesDict:
        pos = lastNodesDict[layer]
        if pos>=median+(median*tolerance[1]):
            layersToRemove.append(layer)   
    return layersToRemove


def cutFor5PrimeIteration(F,sequences):
    cut_sequences = []
    for i,seq in enumerate(sequences):
        layer = i+1
        positions_of_nodes_in_layer =[F.node[x]['location'] for x in F.nodes() if F.node[x]['layer']==layer]
        positions_of_nodes_in_layer.sort()
        if positions_of_nodes_in_layer: #If the layer contains any nodes   
            firstNode = positions_of_nodes_in_layer[0]
            sub_seq = seq[0:firstNode-1] 
            cut_sequences.append(sub_seq)
        else:
            cut_sequences.append(seq)
    return cut_sequences

def cutFor3PrimeIteration(F,sequences):
    marked_sequences = []
    for i,seq in enumerate(sequences):
        layer = i+1
        positions_of_nodes_in_layer =[F.node[x]['end'] for x in F.nodes() if F.node[x]['layer']==layer]
        positions_of_nodes_in_layer.sort()
        if positions_of_nodes_in_layer: #If the layer contains any nodes     
            lastNode = positions_of_nodes_in_layer[-1]
            marked_sequences.append((seq[lastNode:],lastNode))
        else:
            marked_sequences.append((seq,0))

    return marked_sequences




def build_LOOM(sequences,seq_lengths,kmers_len,stopw,hsps,solver,prune,min_depth,shorttandem,maxedges,noconstraints):
    LOOM = nx.DiGraph()
    LOOM_at_Levels={}

    Complex_Graph = nx.DiGraph()

     
    level = len(sequences)#currently as in, cannot skip layers without messing layer IDs

    
    #Build sparse graph with simple paths of longer kmers - (shorttandem is boolean option to also exclude complex paths)
    for kl in range(kmers_len,stopw,-1):

        F = build_graph(sequences,seq_lengths,hsps,kl,prune)
        if len(LOOM.edges)>0 and noconstraints==False:
            F = remove_intersecting_edges(F,LOOM)
      


        C = nx.DiGraph() # an empty graph which is populated during recursive segmentation
        FinalGraph = recursive_segmentation([[F,C,Complex_Graph]],kl,solver,LOOM_at_Levels,prune,level,min_depth,maxedges,noconstraints)

        if len(FinalGraph.nodes())>0:
            Complex_Graph = nx.DiGraph()
            if noconstraints:
                paths_to_add = getSimplePaths(FinalGraph)
                Complex_Graph = add_paths_to_graph(Complex_Graph,FinalGraph,paths_to_add)
            else:
                paths_to_add = getSimplePaths(FinalGraph)
                LOOM = add_paths_to_graph(LOOM,FinalGraph,paths_to_add)
                LOOM_at_Levels=get_graphs_per_level(LOOM,2)

            #If user has set option to exclude complex paths of longer k-mers then do not add them to graph
            if not shorttandem:
                paths_to_add = getComplexPaths(FinalGraph)
                Complex_Graph = add_paths_to_graph(Complex_Graph,FinalGraph,paths_to_add)

        
        
    #Run framework for kmers=stopw in length    
    kl=stopw
    F = build_graph(sequences,seq_lengths,hsps,kl,prune)
    
    
    
    if len(LOOM.edges)>0 and noconstraints==False:
        F = remove_intersecting_edges(F,LOOM)

    C = nx.DiGraph() # an empty graph which is populated during recursive segmentation
    FinalGraph = recursive_segmentation([[F,C,Complex_Graph]],kl,solver,LOOM_at_Levels,prune,level,min_depth,maxedges,noconstraints)
    
    LOOM.add_nodes_from(FinalGraph.nodes(data=True))
    LOOM.add_edges_from(FinalGraph.edges())
    #print [x for x in LOOM.edges() if LOOM.node[x[0]]['layer']==1]

    return LOOM


def get_extenstion_barriers_start(SG):

    graph_layers = list(dict.fromkeys([SG.node[x]['layer'] for x in SG.nodes()]))
    graph_layers.sort()
    barriers = {}
    for i in range(len(graph_layers)-1):
        start_edges1 = [SG.node[x[0]]['location'] for x in SG.edges() if SG.node[x[0]]['layer']==graph_layers[i]]
        start_edges2 = [SG.node[x[1]]['location'] for x in SG.edges() if SG.node[x[0]]['layer']==graph_layers[i]]
        barriers[(graph_layers[i],graph_layers[i+1])] = [start_edges1,start_edges2]

    return barriers


def get_extenstion_barriers_end(SG):

    graph_layers = list(dict.fromkeys([SG.node[x]['layer'] for x in SG.nodes()]))
    graph_layers.sort()
    barriers = {}
    for i in range(len(graph_layers)-1):
        end_edges1 = [SG.node[x[0]]['end'] for x in SG.edges() if SG.node[x[0]]['layer']==graph_layers[i]]
        end_edges2 = [SG.node[x[1]]['end'] for x in SG.edges() if SG.node[x[0]]['layer']==graph_layers[i]]
        barriers[(graph_layers[i],graph_layers[i+1])] = [end_edges1,end_edges2]

    return barriers








def refine_LOOM(SG,sequences,min_depth):


    print "Refining Graph..."

    extenstion_barriers_start = get_extenstion_barriers_start(SG)
    extenstion_barriers_end = get_extenstion_barriers_end(SG)


    graph_layers = list(dict.fromkeys([SG.node[x]['layer'] for x in SG.nodes()]))
    graph_layers.sort()
    last_layer = graph_layers[-1]
    connected_bases_end = {l:[] for l in graph_layers}
    connected_bases_start = {l:[] for l in graph_layers}
    for layer in graph_layers:
        #connected_bases_end[layer] = [SG.node[x]['end'] for x in SG.nodes() if SG.node[x]['layer']==layer]
        #connected_bases_start[layer] = [SG.node[x]['location'] for x in SG.nodes() if SG.node[x]['layer']==layer]
        #for pos in node_positions:
            #connected_bases[layer].extend(range(pos[0],pos[1]+1))
        nodes_in_layer = [(SG.node[x]['location'],SG.node[x]['end'],SG.node[x]['seq']) for x in SG.nodes() if SG.node[x]['layer']==layer]
        nodes_in_layer = sorted(nodes_in_layer,key=lambda x: (x[0],-len(x[2])))
        #print layer,nodes_in_layer
        combined_nodes = combine_kmers_in_layer(nodes_in_layer)
        
        connected_bases_end[layer] = [x[1] for x in combined_nodes]
        connected_bases_start[layer] = [x[0] for x in combined_nodes]

        connected_bases_end[layer] = list(dict.fromkeys(connected_bases_end[layer]))
        connected_bases_start[layer] = list(dict.fromkeys(connected_bases_start[layer]))

    #get individual simple and complex paths of SG, sort by location.     
    simple_paths = getSimplePaths(SG)
    all_simple_paths = []
    for kmer in simple_paths:
        for path in simple_paths[kmer]:
            all_simple_paths.append(path)
    all_simple_paths.sort(key=lambda x: SG.node[x[0][0]]['location'])

    complex_paths = getComplexPaths(SG)
    all_complex_paths = []
    for kmer in complex_paths:
        for path in complex_paths[kmer]:
            all_complex_paths.append(path)
    all_complex_paths.sort(key=lambda x: SG.node[x[0][0]]['location'])


    #process simple paths
    simple_extend = nx.DiGraph()
    for path in all_simple_paths:
        path.sort(key=lambda x: SG.node[x[0]]['layer'])
        connect = True
        for edge in path:
            start1 = SG.node[edge[0]]['location']
            start2 = SG.node[edge[1]]['location']
            layer1 = SG.node[edge[0]]['layer']
            layer2 = SG.node[edge[1]]['layer']

            Allowed = True
            extenstion_barriers = extenstion_barriers_end[(layer1,layer2)]

            if start1-1 in extenstion_barriers[0] and start2-1 not in extenstion_barriers[1]:
                Allowed = False
            elif start1-1 not in extenstion_barriers[0] and start2-1 in extenstion_barriers[1]:
                Allowed = False

 
            if (start1==1 or start2==1) and Allowed:
                connect=False
            else:
                end1 = SG.node[edge[0]]['end']
                end2 = SG.node[edge[1]]['end']
                        

                base1 = sequences[layer1-1][start1-2]
                base2 = sequences[layer2-1][start2-2]            
                if base1==base2:
                    if start1-1 not in connected_bases_end[layer1] and start2-1 not in connected_bases_end[layer2]:
                        seq1 = base1+SG.node[edge[0]]['seq']
                        seq2 = base2+SG.node[edge[1]]['seq']
                        posx1 = SG.node[edge[0]]['posx']
                        posx2 = SG.node[edge[1]]['posx']                        
                        simple_extend.add_node(str(layer1)+'_'+str(start1-1)+'_'+seq1,layer=layer1,location = start1-1,seq=seq1,end=end1,posx=posx1)
                        simple_extend.add_node(str(layer2)+'_'+str(start2-1)+'_'+seq2,layer=layer2,location = start2-1,seq=seq2,end=end2,posx=posx2)
                        simple_extend.add_edge(str(layer1)+'_'+str(start1-1)+'_'+seq1,str(layer2)+'_'+str(start2-1)+'_'+seq2)
                    else:
                        connect = False
   
                else:
                    connect = False

            if not connect:
                break

    #simple_extend = remove_intersecting_edges(simple_extend,SG)
    #simple_extend = prune_graph_level_plus(simple_extend,min_depth)
    #SG.add_nodes_from(simple_extend.nodes()(data=True))
    #SG.add_edges_from(simple_extend.edges())

    for node in simple_extend.nodes():
        layer = simple_extend.node[node]['layer']
        start = simple_extend.node[node]['location']
        if start+1 in connected_bases_start[layer]:
            connected_bases_start[layer].append(start)

    #simple_extend = nx.DiGraph()
    for path in all_simple_paths:
        path.sort(key=lambda x: SG.node[x[0]]['layer'])
        connect = True
        for edge in path:
            end1 = SG.node[edge[0]]['end']
            end2 = SG.node[edge[1]]['end']
            layer1 = SG.node[edge[0]]['layer']
            layer2 = SG.node[edge[1]]['layer']

            Allowed = True
            extenstion_barriers = extenstion_barriers_start[(layer1,layer2)]

            if end1+1 in extenstion_barriers[0] and end2+1 not in extenstion_barriers[1]:
                Allowed = False
            elif end1+1 not in extenstion_barriers[0] and end2+1 in extenstion_barriers[1]:
                Allowed = False


            if (end1==len(sequences[layer1-1]) or end2==len(sequences[layer2-1])) and Allowed:
                connect = False
        
            else:
                base1 = sequences[layer1-1][end1]
                base2 = sequences[layer2-1][end2]      
                if base1==base2:
                    if end1+1 not in connected_bases_start[layer1] and end2+1 not in connected_bases_start[layer2]:
                        seq1 = SG.node[edge[0]]['seq']+base1
                        seq2 = SG.node[edge[1]]['seq']+base2
                        start1 = SG.node[edge[0]]['location']
                        start2 = SG.node[edge[1]]['location']
                        simple_extend.add_node(str(layer1)+'_'+str(start1)+'_'+seq1,layer=layer1,location = start1,seq=seq1,end=end1+1,posx=float(end1+1)/len(sequences[layer1-1]))
                        simple_extend.add_node(str(layer2)+'_'+str(start2)+'_'+seq2,layer=layer2,location = start2,seq=seq2,end=end2+1,posx=float(end2+1)/len(sequences[layer2-1]))
                        simple_extend.add_edge(str(layer1)+'_'+str(start1)+'_'+seq1,str(layer2)+'_'+str(start2)+'_'+seq2)
                    else:
                        connect = False
   
                else:
                    connect = False

            if not connect:
                break

    
    simple_extend = remove_intersecting_edges(simple_extend,SG)
    simple_extend = prune_graph_level_plus(simple_extend,min_depth)
    SG.add_nodes_from(simple_extend.nodes()(data=True))
    SG.add_edges_from(simple_extend.edges())

    for node in simple_extend.nodes():
        layer = simple_extend.node[node]['layer']
        end = simple_extend.node[node]['end']
        if end-1 in connected_bases_end[layer]:
            connected_bases_end[layer].append(end)


    local_connect = False
    EP= nx.DiGraph()

    #process complex paths
    for path in all_complex_paths:
        E = nx.DiGraph()
        C = nx.DiGraph()

        nodes_per_layer = {}
        for edge in [e for e in path if SG.node[e[0]]['layer']==1]:
            node = edge[0]
            node2 = edge[1]
            layer1 = SG.node[edge[0]]['layer']
            layer2 = SG.node[edge[1]]['layer']
            start1 = SG.node[node]['location']
            if start1!=1:
                end1 = SG.node[node]['end']
                base1 = sequences[layer1-1][start1-2]          
                if start1-1 not in connected_bases_end[layer1]:
                    seq1 = base1+SG.node[node]['seq']
                    posx1 = SG.node[node]['posx']

                    start2 = SG.node[node2]['location']

                    Allowed = True
                    extenstion_barriers = extenstion_barriers_end[(layer1,layer2)]

                    if start1-1 in extenstion_barriers[0] and start2-1 not in extenstion_barriers[1]:
                        Allowed = False
                    elif start1-1 not in extenstion_barriers[0] and start2-1 in extenstion_barriers[1]:
                        Allowed = False




                    if start2!=1 and Allowed:
                        end2 = SG.node[node2]['end']
                        base2 = sequences[layer2-1][start2-2]          
                        if start2-1 not in connected_bases_end[layer2] and base1==base2:
                            seq2 = base2+SG.node[node2]['seq']
                            posx2 = SG.node[node2]['posx']
                            E.add_node(str(layer1)+'_'+str(start1-1)+'_'+seq1,layer=layer1,location = start1-1,seq=seq1,end=end1,posx=posx1)
                            E.add_node(str(layer2)+'_'+str(start2-1)+'_'+seq2,layer=layer2,location = start2-1,seq=seq2,end=end2,posx=posx2)
                            E.add_edge(str(layer1)+'_'+str(start1-1)+'_'+seq1,str(layer2)+'_'+str(start2-1)+'_'+seq2)

                            C.add_node(node,layer=layer1,location=start1,seq=SG.node[node]['seq'],end=end1,posx=posx1)
                            C.add_node(node2,layer=layer2,location=start2,seq=SG.node[node2]['seq'],end=end2,posx=posx2)
                            C.add_edge(node,node2)      
                            local_connect = True

        i = 2
        while i<last_layer and local_connect:
            local_connect = False
            nodes1 = [x for x in C.nodes() if C.node[x]['layer']==i]
            edges2 = [x for x in path if x[0] in nodes1]
            layer1 = i
            layer2 = i+1
            for edge in edges2:
                node1 = edge[0]
                node2 = edge[1]
                start1 = SG.node[node1]['location']
                end1 = SG.node[node1]['end']
                base1 = sequences[layer1-1][start1-2]          
                seq1 = base1+SG.node[node]['seq']
                posx1 = SG.node[node]['posx']       
                start2 = SG.node[node2]['location']

                Allowed = True
                extenstion_barriers = extenstion_barriers_end[(layer1,layer2)]

                if start1-1 in extenstion_barriers[0] and start2-1 not in extenstion_barriers[1]:
                    Allowed = False
                elif start1-1 not in extenstion_barriers[0] and start2-1 in extenstion_barriers[1]:
                    Allowed = False

                if start2!=1 and Allowed:
                    end2 = SG.node[node2]['end']
                    base2 = sequences[layer2-1][start2-2]
                    seq2 = base2+SG.node[node2]['seq']          
                    if start2-1 not in connected_bases_end[layer2] and seq1==seq2:
                        posx2 = SG.node[node2]['posx']
                        E.add_node(str(layer1)+'_'+str(start1-1)+'_'+seq1,layer=layer1,location = start1-1,seq=seq1,end=end1,posx=posx1)
                        E.add_node(str(layer2)+'_'+str(start2-1)+'_'+seq2,layer=layer2,location = start2-1,seq=seq2,end=end2,posx=posx2)
                        E.add_edge(str(layer1)+'_'+str(start1-1)+'_'+seq1,str(layer2)+'_'+str(start2-1)+'_'+seq2)
                        #for cn in [x for x in C.nodes() if C.node[x]['layer']==layer1]:
                        C.add_node(node2,layer=layer2,location = start2,seq=SG.node[node2]['seq'],end=end2,posx=posx2)
                        C.add_edge(node1,node2)      
                        local_connect = True
            if not local_connect:
                #print 'done'
                #add remaining layers to C
                for j in range(i,last_layer):
                    nodes1 = [x for x in C.nodes() if C.node[x]['layer']==j]
                    edges2 = [x for x in path if x[0] in nodes1]
                    for edge in edges2:
                        node = edge[1]
                        pn = edge[0]
                        C.add_node(node,layer=j+1,location=SG.node[node]['location'],seq=SG.node[node]['seq'],end=SG.node[node]['end'],posx=SG.node[node]['posx'])
                        C.add_edge(pn,node)
            i+=1
        E = prune_graph_level_plus(E,min_depth)


        if len(E.nodes())>0:
            to_remove = []
            for edge in path:
                to_remove.append(edge[0])
                to_remove.append(edge[1])
            SG.remove_nodes_from(to_remove)
            depth_nodes = {E.node[x]['seq']:[] for x in E.nodes()}
            if len(depth_nodes)==1:
                E = remove_intersecting_edges(E,SG)
                E = prune_graph_level_plus(E,min_depth)
                EP.add_nodes_from(E.nodes(data=True))
                EP.add_edges_from(E.edges())
                SG.add_nodes_from(C.nodes(data=True))
                SG.add_edges_from(C.edges())

            else:
                #Select the most conserved option - prevents mix up of order
                for x in E.nodes():
                    depth_nodes[E.node[x]['seq']].append(E.node[x]['layer'])
                for s in depth_nodes:
                    depth_nodes[s].sort()
                    depth_nodes[s] = depth_nodes[s][-1]
                 
                kmer = max(depth_nodes,key=depth_nodes.get)

                to_remove = [x for x in E.nodes if E.node[x]['seq']!=kmer]
                to_removeC = [str(E.node[x]['layer'])+'_'+str(E.node[x]['location']+1)+'_'+kmer[1:] for x in to_remove]
              
                E.remove_nodes_from(to_remove)
                C.remove_nodes_from(to_removeC)

                E = remove_intersecting_edges(E,SG)
                E = prune_graph_level_plus(E,min_depth)

                EP.add_nodes_from(E.nodes(data=True))
                EP.add_edges_from(E.edges())
                SG.add_nodes_from(C.nodes(data=True))
                SG.add_edges_from(C.edges())
                

    if len(EP.nodes)>0: 

        complex_paths = getComplexPaths(EP)
        all_complex_paths = []
        for kmer in complex_paths:
            for path in complex_paths[kmer]:
                all_complex_paths.append(path)
        all_complex_paths.sort(key=lambda x: EP.node[x[0][0]]['location'])
    
        for node in EP.nodes():
            layer = EP.node[node]['layer']
            start = EP.node[node]['location']
            if start+1 in connected_bases_start[layer]:
                connected_bases_start[layer].append(start)
    
    
        local_connect = False
        #process complex paths
        for path in all_complex_paths:
            E = nx.DiGraph()
            C = nx.DiGraph()
            nodes_per_layer = {}
            for edge in [e for e in path if EP.node[e[0]]['layer']==1]:
                node = edge[0]
                node2 = edge[1]
                layer1 = EP.node[edge[0]]['layer']
                layer2 = EP.node[edge[1]]['layer']
                start1 = EP.node[node]['location']
                end1 = EP.node[node]['end']
                if end1!=len(sequences[layer1-1]):
                    base1 = sequences[layer1-1][end1]          
                    if end1+1 not in connected_bases_end[layer1]:
                        seq1 = EP.node[node]['seq']+base1
                        posx1 = EP.node[node]['posx']
    
                        start2 = EP.node[node2]['location']
                        end2 = EP.node[node2]['end']
    
                        Allowed = True
                        extenstion_barriers = extenstion_barriers_start[(layer1,layer2)]
    
                        if end1+1 in extenstion_barriers[0] and end2+1 not in extenstion_barriers[1]:
                            Allowed = False
                        elif end1+1 not in extenstion_barriers[0] and end2+1 in extenstion_barriers[1]:
                            Allowed = False   
                  
                        if end2!=len(sequences[layer2-1]) and Allowed:                        
                            base2 = sequences[layer2-1][end2]          
                            if end2+1 not in connected_bases_start[layer2] and base1==base2:
                                seq2 = EP.node[node2]['seq']+base2
                                posx2 = EP.node[node2]['posx']
    
                                E.add_node(str(layer1)+'_'+str(start1)+'_'+seq1,layer=layer1,location = start1,seq=seq1,end=end1+1,posx=float(end1+1)/len(sequences[layer1-1]))
                                E.add_node(str(layer2)+'_'+str(start2)+'_'+seq2,layer=layer2,location = start2,seq=seq2,end=end2+1,posx=float(end2+1)/len(sequences[layer2-1]))
                                E.add_edge(str(layer1)+'_'+str(start1)+'_'+seq1,str(layer2)+'_'+str(start2)+'_'+seq2)
    
                                C.add_node(node,layer=layer1,location=start1,seq=EP.node[node]['seq'],end=end1,posx=posx1)
                                C.add_node(node2,layer=layer2,location=start2,seq=EP.node[node2]['seq'],end=end2,posx=posx2)
                                C.add_edge(node,node2)      
                                local_connect = True
    
    
            i = 2
            while i<last_layer and local_connect:
                local_connect = False
                nodes1 = [x for x in C.nodes() if C.node[x]['layer']==i]
                edges2 = [x for x in path if x[0] in nodes1]
                layer1 = i
                layer2 = i+1
                for edge in edges2:
                    node1 = edge[0]
                    node2 = edge[1]
                    start1 = EP.node[node1]['location']
                    end1 = EP.node[node1]['end']
                    base1 = sequences[layer1-1][end1]          
                    seq1 = EP.node[node]['seq']+base1
                    posx1 = EP.node[node]['posx']       
                    start2 = EP.node[node2]['location']
                    end2 = EP.node[node2]['end']
    
                    Allowed = True
                    extenstion_barriers = extenstion_barriers_start[(layer1,layer2)]
    
                    if end1+1 in extenstion_barriers[0] and end2+1 not in extenstion_barriers[1]:
                        Allowed = False
                    elif end1+1 not in extenstion_barriers[0] and end2+1 in extenstion_barriers[1]:
                        Allowed = False
    
    
    
                    if end2!=len(sequences[layer2-1]) and Allowed:
                        base2 = sequences[layer2-1][end2]
                        seq2 = EP.node[node2]['seq']+base2          
                        if end2+1 not in connected_bases_start[layer2] and seq1==seq2:
                            posx2 = EP.node[node2]['posx']
                            E.add_node(str(layer1)+'_'+str(start1)+'_'+seq1,layer=layer1,location = start1,seq=seq1,end=end1+1,posx=float(end1+1)/len(sequences[layer1-1]))
                            E.add_node(str(layer2)+'_'+str(start2)+'_'+seq2,layer=layer2,location = start2,seq=seq2,end=end2+1,posx=float(end2+1)/len(sequences[layer2-1]))
                            E.add_edge(str(layer1)+'_'+str(start1)+'_'+seq1,str(layer2)+'_'+str(start2)+'_'+seq2)
    
                            C.add_node(node2,layer=layer2,location = start2,seq=EP.node[node2]['seq'],end=end2,posx=posx2)
                            C.add_edge(node1,node2)      
                            local_connect = True
                if not local_connect:
                    #print 'done'
                    #add remaining layers to C
                    for j in range(i,last_layer):
                        nodes1 = [x for x in C.nodes() if C.node[x]['layer']==j]
                        edges2 = [x for x in path if x[0] in nodes1]
                        for edge in edges2:
                            node = edge[1]
                            pn = edge[0]
                            C.add_node(node,layer=j+1,location=EP.node[node]['location'],seq=EP.node[node]['seq'],end=EP.node[node]['end'],posx=EP.node[node]['posx'])
                            C.add_edge(pn,node)
                i+=1
    
            E = prune_graph_level_plus(E,min_depth)
            #print C.nodes()
            #print E.nodes()
            #print EP.nodes()
    
            if len(E.nodes())>0:
                to_remove = []
                for edge in path:
                    to_remove.append(edge[0])
                    to_remove.append(edge[1])
                EP.remove_nodes_from(to_remove)
                depth_nodes = {E.node[x]['seq']:[] for x in E.nodes()}
                if len(depth_nodes)==1:
                    E = remove_intersecting_edges(E,SG)
                    E = prune_graph_level_plus(E,min_depth)
                    #to_removeEP = [x for x in EP.nodes() if str(EP.node[x]['layer'])+'_'+str(EP.node[x]['location']+1)+'_'+EP.node[x]['seq'][1:] not in C.nodes()]
                    #EP.remove_nodes_from(to_removeEP)
                    #EP.add_nodes_from(E.nodes(data=True))
                    #EP.add_edges_from(E.edges())
                    SG.add_nodes_from(E.nodes(data=True))
                    SG.add_edges_from(E.edges())                    
                    SG.add_nodes_from(C.nodes(data=True))
                    SG.add_edges_from(C.edges())
    
                else:
                    #Select the most conserved option - prevents mix up of order
                    for x in E.nodes():
                        depth_nodes[E.node[x]['seq']].append(E.node[x]['layer'])
                    for s in depth_nodes:
                        depth_nodes[s].sort()
                        depth_nodes[s] = depth_nodes[s][-1]
                     
                    kmer = max(depth_nodes,key=depth_nodes.get)
    
                    to_remove = [x for x in E.nodes if E.node[x]['seq']!=kmer]
                    to_removeC = [str(E.node[x]['layer'])+'_'+str(E.node[x]['location']+1)+'_'+kmer[0:-1] for x in to_remove]
                  
                    E.remove_nodes_from(to_remove)
                    C.remove_nodes_from(to_removeC)
                    
                    #to_removeEP = [x for x in EP.nodes() if str(EP.node[x]['layer'])+'_'+str(EP.node[x]['location']+1)+'_'+kmer[1:] not in C.nodes()]
                    #EP.remove_nodes_from(to_removeEP)
    
                    E = remove_intersecting_edges(E,SG)
                    E = prune_graph_level_plus(E,min_depth)
    
                    #EP.add_nodes_from(E.nodes(data=True))
                    #EP.add_edges_from(E.edges())
                    SG.add_nodes_from(E.nodes(data=True))
                    SG.add_edges_from(E.edges())                    
                    SG.add_nodes_from(C.nodes(data=True))
                    SG.add_edges_from(C.edges())
    else:
        
        complex_paths = getComplexPaths(SG)
        all_complex_paths = []
        for kmer in complex_paths:
            for path in complex_paths[kmer]:
                all_complex_paths.append(path)
        all_complex_paths.sort(key=lambda x: SG.node[x[0][0]]['location'])

        local_connect = False
        #process complex paths
        for path in all_complex_paths:
            E = nx.DiGraph()
            C = nx.DiGraph()
            nodes_per_layer = {}
            for edge in [e for e in path if SG.node[e[0]]['layer']==1]:
                node = edge[0]
                node2 = edge[1]
                layer1 = SG.node[edge[0]]['layer']
                layer2 = SG.node[edge[1]]['layer']
                start1 = SG.node[node]['location']
                end1 = SG.node[node]['end']
                if end1!=len(sequences[layer1-1]):
                    base1 = sequences[layer1-1][end1]          
                    if end1+1 not in connected_bases_end[layer1]:
                        seq1 = SG.node[node]['seq']+base1
                        posx1 = SG.node[node]['posx']
    
                        start2 = SG.node[node2]['location']
                        end2 = SG.node[node2]['end']
    
                        Allowed = True
                        extenstion_barriers = extenstion_barriers_start[(layer1,layer2)]
    
                        if end1+1 in extenstion_barriers[0] and end2+1 not in extenstion_barriers[1]:
                            Allowed = False
                        elif end1+1 not in extenstion_barriers[0] and end2+1 in extenstion_barriers[1]:
                            Allowed = False   
                  
                        if end2!=len(sequences[layer2-1]) and Allowed:                        
                            base2 = sequences[layer2-1][end2]          
                            if end2+1 not in connected_bases_start[layer2] and base1==base2:
                                seq2 = SG.node[node2]['seq']+base2
                                posx2 = SG.node[node2]['posx']
    
                                E.add_node(str(layer1)+'_'+str(start1)+'_'+seq1,layer=layer1,location = start1,seq=seq1,end=end1+1,posx=float(end1+1)/len(sequences[layer1-1]))
                                E.add_node(str(layer2)+'_'+str(start2)+'_'+seq2,layer=layer2,location = start2,seq=seq2,end=end2+1,posx=float(end2+1)/len(sequences[layer2-1]))
                                E.add_edge(str(layer1)+'_'+str(start1)+'_'+seq1,str(layer2)+'_'+str(start2)+'_'+seq2)
    
                                C.add_node(node,layer=layer1,location=start1,seq=SG.node[node]['seq'],end=end1,posx=posx1)
                                C.add_node(node2,layer=layer2,location=start2,seq=SG.node[node2]['seq'],end=end2,posx=posx2)
                                C.add_edge(node,node2)      
                                local_connect = True
    
    
            i = 2
            while i<last_layer and local_connect:
                local_connect = False
                nodes1 = [x for x in C.nodes() if C.node[x]['layer']==i]
                edges2 = [x for x in path if x[0] in nodes1]
                layer1 = i
                layer2 = i+1
                for edge in edges2:
                    node1 = edge[0]
                    node2 = edge[1]
                    start1 = SG.node[node1]['location']
                    end1 = SG.node[node1]['end']
                    base1 = sequences[layer1-1][end1]          
                    seq1 = SG.node[node]['seq']+base1
                    posx1 = SG.node[node]['posx']       
                    start2 = SG.node[node2]['location']
                    end2 = SG.node[node2]['end']
    
                    Allowed = True
                    extenstion_barriers = extenstion_barriers_start[(layer1,layer2)]
    
                    if end1+1 in extenstion_barriers[0] and end2+1 not in extenstion_barriers[1]:
                        Allowed = False
                    elif end1+1 not in extenstion_barriers[0] and end2+1 in extenstion_barriers[1]:
                        Allowed = False
    
    
    
                    if end2!=len(sequences[layer2-1]) and Allowed:
                        base2 = sequences[layer2-1][end2]
                        seq2 = SG.node[node2]['seq']+base2          
                        if end2+1 not in connected_bases_start[layer2] and seq1==seq2:
                            posx2 = SG.node[node2]['posx']
                            E.add_node(str(layer1)+'_'+str(start1)+'_'+seq1,layer=layer1,location = start1,seq=seq1,end=end1+1,posx=float(end1+1)/len(sequences[layer1-1]))
                            E.add_node(str(layer2)+'_'+str(start2)+'_'+seq2,layer=layer2,location = start2,seq=seq2,end=end2+1,posx=float(end2+1)/len(sequences[layer2-1]))
                            E.add_edge(str(layer1)+'_'+str(start1)+'_'+seq1,str(layer2)+'_'+str(start2)+'_'+seq2)
    
                            C.add_node(node2,layer=layer2,location = start2,seq=SG.node[node2]['seq'],end=end2,posx=posx2)
                            C.add_edge(node1,node2)      
                            local_connect = True
                if not local_connect:
                    #print 'done'
                    #add remaining layers to C
                    for j in range(i,last_layer):
                        nodes1 = [x for x in C.nodes() if C.node[x]['layer']==j]
                        edges2 = [x for x in path if x[0] in nodes1]
                        for edge in edges2:
                            node = edge[1]
                            pn = edge[0]
                            C.add_node(node,layer=j+1,location=SG.node[node]['location'],seq=SG.node[node]['seq'],end=SG.node[node]['end'],posx=SG.node[node]['posx'])
                            C.add_edge(pn,node)
                i+=1
    
            E = prune_graph_level_plus(E,min_depth)
            #print C.nodes()
            #print E.nodes()
            #print EP.nodes()
    
            if len(E.nodes())>0:
                to_remove = []
                for edge in path:
                    to_remove.append(edge[0])
                    to_remove.append(edge[1])
                SG.remove_nodes_from(to_remove)
                depth_nodes = {E.node[x]['seq']:[] for x in E.nodes()}
                if len(depth_nodes)==1:
                    E = remove_intersecting_edges(E,SG)
                    E = prune_graph_level_plus(E,min_depth)
                    #to_removeEP = [x for x in EP.nodes() if str(EP.node[x]['layer'])+'_'+str(EP.node[x]['location']+1)+'_'+EP.node[x]['seq'][1:] not in C.nodes()]
                    #EP.remove_nodes_from(to_removeEP)
                    #EP.add_nodes_from(E.nodes(data=True))
                    #EP.add_edges_from(E.edges())
                    SG.add_nodes_from(E.nodes(data=True))
                    SG.add_edges_from(E.edges())                    
                    SG.add_nodes_from(C.nodes(data=True))
                    SG.add_edges_from(C.edges())
    
                else:
                    #Select the most conserved option - prevents mix up of order
                    for x in E.nodes():
                        depth_nodes[E.node[x]['seq']].append(E.node[x]['layer'])
                    for s in depth_nodes:
                        depth_nodes[s].sort()
                        depth_nodes[s] = depth_nodes[s][-1]
                     
                    kmer = max(depth_nodes,key=depth_nodes.get)
    
                    to_remove = [x for x in E.nodes if E.node[x]['seq']!=kmer]
                    to_removeC = [str(E.node[x]['layer'])+'_'+str(E.node[x]['location']+1)+'_'+kmer[0:-1] for x in to_remove]
                  
                    E.remove_nodes_from(to_remove)
                    C.remove_nodes_from(to_removeC)
                    
                    #to_removeEP = [x for x in EP.nodes() if str(EP.node[x]['layer'])+'_'+str(EP.node[x]['location']+1)+'_'+kmer[1:] not in C.nodes()]
                    #EP.remove_nodes_from(to_removeEP)
    
                    E = remove_intersecting_edges(E,SG)
                    E = prune_graph_level_plus(E,min_depth)
    
                    #EP.add_nodes_from(E.nodes(data=True))
                    #EP.add_edges_from(E.edges())
                    SG.add_nodes_from(E.nodes(data=True))
                    SG.add_edges_from(E.edges())                      
                    SG.add_nodes_from(C.nodes(data=True))
                    SG.add_edges_from(C.edges())



    #SG.add_nodes_from(EP.nodes(data=True))
    #SG.add_edges_from(EP.edges())       
    return SG
 
def start_lncLOOM(sequences,seq_lengths,kmers_len,stopw,hsps,solver,prune,min_depth,shorttandem,maxedges,noconstraints,tolerance):

    LOOM_at_Levels = {}
    LOOM = build_LOOM(sequences,seq_lengths,kmers_len,stopw,hsps,solver,prune,min_depth,shorttandem,maxedges,noconstraints)

    #prune LOOM to exclude nodes that do not have min_depth
    LOOM = prune_graph_level_plus(LOOM,min_depth)

    e = [x for x in LOOM.edges if LOOM.node[x[0]]['layer']==2]

    if len(LOOM.nodes())>0:
        #for i in range(5):
        LOOM = refine_LOOM(LOOM,sequences,min_depth)

  


    LOOM_at_Levels = get_graphs_per_level(LOOM,min_depth)
    #Main_LOOM_Level = LOOM_at_Levels.copy()
    Main_LOOM_Level = get_graphs_per_level(LOOM,min_depth)

    #Main_LOOM = LOOM.copy()
    #Main_LOOM = prune_graph_level_plus(Main_LOOM,min_depth)    

    LOOM5_Levels = {}
    LOOM3_Levels = {}
    for level in LOOM_at_Levels:
        LOOM5_Levels[level] = nx.DiGraph()
        LOOM3_Levels[level] = nx.DiGraph()

    Corrected5 = nx.DiGraph()
    Corrected3 = nx.DiGraph()

    details5 = []
    details3 = []

    total_seq = len(sequences)

    #next build a series of graphs that skip layers


    #Next build graphs for sequences that have longer 5' or 3' end (only do this if the top layer has longer 5' or 3' end)
    layer_ids = {}
    if len(LOOM.nodes())>0:

        #Determine iteration of graph building for 5' and 3' caps based on first and last nodes    
        remove5prime = calculateMissing5Prime(LOOM,total_seq,tolerance)
        remove3prime = calculateMissing3Prime(LOOM,total_seq,tolerance)

        if 1 not in remove5prime and len(remove5prime)>0:
            print "Calculating 5 Prime Graph"
            marked5_sequences = cutFor5PrimeIteration(LOOM,sequences)
            five_prime_seqs = []
            five_prime_hsps = []
            five_prime_lens = []
            layer_ids = {}
            count_5prime = 0
            for i in range(1,total_seq+1):
                if not i in remove5prime:
                    five_prime_seqs.append(marked5_sequences[i-1])
                    five_prime_lens.append(len(marked5_sequences[i-1]))
                    five_prime_hsps.append([])
                    count_5prime+=1
                    layer_ids[count_5prime]=i
                    details5.append((i,len(marked5_sequences[i-1])))


             
            LOOM_5Prime = build_LOOM(five_prime_seqs,five_prime_lens,kmers_len,stopw,five_prime_hsps,solver,prune,2,shorttandem,maxedges,noconstraints)
            if len(LOOM_5Prime.nodes())>0:
                LOOM_5Prime = refine_LOOM(LOOM_5Prime,five_prime_seqs,2)
            if len(LOOM_5Prime.nodes())>0:
                #add nodes to main LOOM: but first edit attributes
                for edge in LOOM_5Prime.edges():
                    node1 = edge[0]
                    node2 = edge[1]
                    layer1 = layer_ids[LOOM_5Prime.node[node1]['layer']]
                    layer2 = layer_ids[LOOM_5Prime.node[node2]['layer']]
                    label1 = str(layer1)+'_'+str(LOOM_5Prime.node[node1]['location'])+'_'+LOOM_5Prime.node[node1]['seq']
                    label2 = str(layer2)+'_'+str(LOOM_5Prime.node[node2]['location'])+'_'+LOOM_5Prime.node[node2]['seq']
                    
                    Corrected5.add_node(label1,layer=layer1,location=LOOM_5Prime.node[node1]['location'],end=LOOM_5Prime.node[node1]['end'],seq=LOOM_5Prime.node[node1]['seq'],posx=(LOOM_5Prime.node[node1]['end'])/float(len(sequences[layer1-1])))
                    Corrected5.add_node(label2,layer=layer2,location=LOOM_5Prime.node[node2]['location'],end=LOOM_5Prime.node[node2]['end'],seq=LOOM_5Prime.node[node2]['seq'],posx=(LOOM_5Prime.node[node2]['end'])/float(len(sequences[layer2-1])))
                    Corrected5.add_edge(label1,label2)
                paths5 = getNonIntersectingPaths(Corrected5)
                #LOOM = add_paths_to_graph(LOOM,Corrected,paths5)

                #Add 3 prime nodes to graph at each level
                to_remove = []
                for kmer in paths5:
                    for path in paths5[kmer]:
                        path.sort(key=lambda x: Corrected5.node[x[0]]['layer'])
                        depth = Corrected5.node[path[-1][1]]['layer']
                        if depth<min_depth:
                            for edge in path:
                                to_remove.append(edge[0])
                                to_remove.append(edge[1])
                        for l in range(depth,min_depth-1,-1):
                            LOOM_at_Levels[l] = add_paths_to_graph(LOOM_at_Levels[l],Corrected5,{kmer:[path]})
                            LOOM5_Levels[l] = add_paths_to_graph(LOOM5_Levels[l],Corrected5,{kmer:[path]})
                            prune_layers = [x for x in LOOM_at_Levels[l].nodes() if LOOM_at_Levels[l].node[x]['layer']>l]
                            LOOM_at_Levels[l].remove_nodes_from(prune_layers)
                            LOOM5_Levels[l].remove_nodes_from(prune_layers)

                Corrected5.remove_nodes_from(to_remove)


        
        if 1 not in remove3prime and len(remove3prime)>0:
            print "Calculating 3 Prime Graph"
            marked3_sequences = cutFor3PrimeIteration(LOOM,sequences)
            three_prime_seqs = []
            three_prime_lens = []
            three_prime_hsps = []
            layer_ids = {}
            count_3prime = 0
            end_of_last_nodes = [x[1] for x in marked3_sequences]
            for i in range(1,total_seq+1):
                if not i in remove3prime:
                    three_prime_seqs.append(marked3_sequences[i-1][0])
                    three_prime_lens.append(len(marked3_sequences[i-1][0]))
                    three_prime_hsps.append([])
                    count_3prime+=1
                    layer_ids[count_3prime]=i
                    details3.append((i,len(marked3_sequences[i-1][0])))

            LOOM_3Prime = build_LOOM(three_prime_seqs,three_prime_lens,kmers_len,stopw,three_prime_hsps,solver,prune,2,shorttandem,maxedges,noconstraints)
            if len(LOOM_3Prime.nodes())>0:
                LOOM_3Prime = refine_LOOM(LOOM_3Prime,three_prime_seqs,2)
            if len(LOOM_3Prime.nodes())>0:
                #add nodes to main LOOM: but first edit attributes
                for edge in LOOM_3Prime.edges():
                    node1 = edge[0]
                    node2 = edge[1]
                    layer1 = layer_ids[LOOM_3Prime.node[node1]['layer']]
                    layer2 = layer_ids[LOOM_3Prime.node[node2]['layer']]
                    label1 = str(layer1)+'_'+str(LOOM_3Prime.node[node1]['location']+end_of_last_nodes[layer1-1])+'_'+LOOM_3Prime.node[node1]['seq']
                    label2 = str(layer2)+'_'+str(LOOM_3Prime.node[node2]['location']+end_of_last_nodes[layer2-1])+'_'+LOOM_3Prime.node[node2]['seq']
                    
                    Corrected3.add_node(label1,layer=layer1,location=LOOM_3Prime.node[node1]['location']+end_of_last_nodes[layer1-1],end=LOOM_3Prime.node[node1]['end']+end_of_last_nodes[layer1-1],seq=LOOM_3Prime.node[node1]['seq'],posx=(LOOM_3Prime.node[node1]['end']+end_of_last_nodes[layer1-1])/float(len(sequences[layer1-1])))
                    Corrected3.add_node(label2,layer=layer2,location=LOOM_3Prime.node[node2]['location']+end_of_last_nodes[layer2-1],end=LOOM_3Prime.node[node2]['end']+end_of_last_nodes[layer2-1],seq=LOOM_3Prime.node[node2]['seq'],posx=(LOOM_3Prime.node[node2]['end']+end_of_last_nodes[layer2-1])/float(len(sequences[layer2-1])))
                    Corrected3.add_edge(label1,label2)
                paths3 = getNonIntersectingPaths(Corrected3)
                #LOOM = add_paths_to_graph(LOOM,Corrected,paths3) 
                #Add 3 prime nodes to graph at each level
                to_remove = []
                for kmer in paths3:
                    for path in paths3[kmer]:
                        path.sort(key=lambda x: Corrected3.node[x[0]]['layer'])
                        depth = Corrected3.node[path[-1][1]]['layer']
                        if depth<min_depth:
                            for edge in path:
                                to_remove.append(edge[0])
                                to_remove.append(edge[1])
                        for l in range(depth,min_depth-1,-1):
                            LOOM_at_Levels[l] = add_paths_to_graph(LOOM_at_Levels[l],Corrected3,{kmer:[path]})
                            LOOM3_Levels[l] = add_paths_to_graph(LOOM3_Levels[l],Corrected3,{kmer:[path]})
                            prune_layers = [x for x in LOOM_at_Levels[l].nodes() if LOOM_at_Levels[l].node[x]['layer']>l]
                            LOOM_at_Levels[l].remove_nodes_from(prune_layers)
                            LOOM3_Levels[l].remove_nodes_from(prune_layers)
                Corrected3.remove_nodes_from(to_remove)

    #LOOM = prune_graph_level_plus(LOOM,min_depth)    
                 
    print "FINAL GRAPH RETURNED"
    print "Number of nodes = "+str(len(LOOM.nodes())+len(Corrected3.nodes())+len(Corrected5.nodes()))
    print "Number of edges = "+str(len(LOOM.edges())+len(Corrected3.edges())+len(Corrected5.nodes()))
    return LOOM_at_Levels,LOOM,Main_LOOM_Level,Corrected5,Corrected3,LOOM5_Levels,LOOM3_Levels,details5,details3






