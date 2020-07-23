import os
import sys
import argparse
import LincMotif as SM
import LincRNA_Blast as LB
import LincRNA_Annotate as LA
import LincRNA_Output as LO
import LincStats as LS
import shutil

              
def wellcome():
    print '==================================================================================================='
    print '                                            LncLOOM                                                ' 
    print '=========================== Ulitsky Lab - Weizmann Institute of Science ==========================='
    print ' Graph-based framework to mine ordered combinations of short conserved kmers in lincRNA sequences  '
    print '===================================================================================================\n'

def get_options():
    parser = argparse.ArgumentParser()
    parser.add_argument("-f","--fasta", dest="fasta", required=True, help="Input a fasta file of lincRNA sequences for motif discovery", default = None)
    parser.add_argument("-s", "--solver", dest="solver", help="Select linear programming solver: GUROBI or CBC: (Default = CBC)", type=str, default = "CBC")
    parser.add_argument("-m", "--multiprocess", dest="multiprocess",help="Split statistical iterations over multiple cores", type=int,default = 4)
    parser.add_argument("-w", "--startw", dest="startw",help="lncLOOM will begin scanning for motifs of length=startw and then iteratively decrease to a length of 6 (default) bases pairs", type=int,default = 15)
    parser.add_argument("-k", "--stopw", dest="stopw",help="lncLOOM will begin scanning for motifs of length=startw and then iteratively decrease to length=stopw", type=int,default = 6)
    parser.add_argument("-d", "--mindepth", dest="mindepth",help="Minimum depth(d) required for motif identifcation. lncLOOM will only identify motifs that appear it at least d layers of sequences", type=int,default = 2)
    parser.add_argument("-r", "--iterations", dest="iterations",help="Define the number of random iterations for statistical analysis", type=int,default = 0)
    parser.add_argument("-p", "--prune", dest="prune",help="Set the maximum number(n) of repeated k-mers. K-mers present more than n times will be removed from graph", type=int,default = 8)  
    parser.add_argument("-b", "--hspblast", dest="hspblast",help="Imposes BLAST HSPs constraints on graph construction: (Default = False)",action="store_true")
    parser.add_argument("-nt", "--shorttandem", dest="shorttandem",help="Excludes tandem repeats (complex paths) found in iterations of larger k: (Default = False)",action="store_true")
    parser.add_argument("-nc", "--noconstraints", dest="noconstraints",help="If set to true, then simple paths calculated from longer kmers will not be used as constraints on edge construction (Default = False)",action="store_true")
    parser.add_argument("-i", "--inputorder", dest="inputorder",help="Graph will be calculated based on the order of sequences in the input fasta file: (Default = False): Sequences reordered by homology",action="store_true")
    parser.add_argument("-sl", "--select", dest="select",help="Enter the name of the sequence to which you like to find conserved motifs", type=str, default = "DEEPEST")
    parser.add_argument("-n", "--pname", dest="projectname",help="Enter a name for the project: (Default = LincMotif)", type=str, default = "LincMotif")
    parser.add_argument("-o", "--outdir", dest="outdir",help="Directory for output (Default: current working directory)", type=str, default = "./")
    parser.add_argument("-x", "--similarity", dest="similarity",help="Excludes sequences that have greater than specified similarity", type=float,default = 100.00)
    parser.add_argument("-t5", "--tol5", dest="tolerance5prime",help="Tolerance step from median postion of first and last nodes to determine exclusion from 5' extension graphs", type=float,default = 0.50)
    parser.add_argument("-t3", "--tol3", dest="tolerance3prime",help="Tolerance step from median postion of first and last nodes to determine exclusion from 3' extension graphs", type=float,default = 0.50)
    parser.add_argument("-e", "--maxedges", dest="maxedges",help="Sets the maximum number of edges allowed on graph. If the graph is too large, a solution is not solved", type=int,default = 1200)
    
   #Add parameters associated with annotation - default is to not annotate
    parser.add_argument("-tscan", "--targetscan", dest="targetscan",help="Annotates with TargetScan: (Default = False)", action="store_true")
    parser.add_argument("-eclip", "--eclip", dest="eclip",help="Upload text file formatted for eCLIP annotation", default = None)
    parser.add_argument("-track", "--track", dest="track",help="Outputs a track bedfile of conserved motifs in query sequence (Default = False)", action="store_true")
    parser.add_argument("-ncol", "--newcolours", dest="newcolours",help="Generates a new random colour set for output: (Default = False)", action="store_true")
    parser.add_argument("-u", "--rna", dest="rna",help="converts DNA sequence to RNA", action="store_true")

    args = parser.parse_args()
    return args


def check_args(args):
    Valid = True
    if args.fasta == None:
        print "ERROR! No fasta file uploaded"
        Valid = False
    elif not os.path.isfile(str(args.fasta)):
        print "ERROR! File not found:"+str(args.fasta)
        Valid = False


    if str((args.solver)).upper()!="GUROBI" and str((args.solver)).upper()!="CBC":
        print "ERROR! Invalid solver:"+args.solver
        Valid = False

    if args.similarity<0 or args.similarity>100:
        print "ERROR! Similarity threshold but be between 0 and 100"
        Valid = False

    if args.tolerance5prime<0:
        print "Tolerance step must be postive"
        Valid = False
    if args.tolerance3prime<0:
        print "Tolerance step must be postive"
        Valid = False


    if args.maxedges<50:
        print "ERROR! Maxedges must be greater than 50"
        Valid = False



    if args.startw < 6:
        print "ERROR! startw cannot be less than 6"
        Valid = False

    if args.stopw < 6 or args.stopw>args.startw:
        print "ERROR! startw cannot be less than 6 and must be less than or equal to startw "
        Valid = False

    if args.mindepth < 2:
        print "ERROR! mindepth cannot be less than 2"
        Valid = False

    if args.prune < 2:
        print "ERROR! prune cannot be less than 2"
        Valid = False

    if args.multiprocess < 1:
        print "ERROR! multiprocess cannot be less than 1"
        Valid = False

    if args.iterations < 0:
        print "ERROR! Number of iterations must be positive"
        Valid = False


    if args.eclip!= None:
        if not os.path.exists(str(args.eclip)):
            print "ERROR! File not found:"+str(args.eclip)
            Valid = False

    if args.track == True:
        #check bedfile_output.txt is found
        if not os.path.exists('./src/for_track_output.txt'):
            print "ERROR! File not found: ./src/for_track_output.txt"
            Valid = False


    return Valid



def main():
    wellcome()
    args=get_options()
    Valid = check_args(args)
    if Valid:
        print "\n==============================================================\nINPUT ARGUMENTS VALID...STARTING PROGRAM...\n"


        #Define variables based on input arguments
        filename = str(args.fasta)
        hspblast = args.hspblast
        shorttandem = args.shorttandem
        inputorder = args.inputorder
        kmers_len = args.startw
        stopw = args.stopw
        multiprocess = args.multiprocess
        min_depth = args.mindepth
        prune = args.prune
        outdir = args.outdir
        project_name = args.projectname
        solver = (args.solver).upper()
        targetScan =args.targetscan
        eCLIP = str(args.eclip)
        select = args.select.upper()
        iterations = args.iterations
        similarity = args.similarity
        newcolours = args.newcolours
        noconstraints = args.noconstraints
        maxedges = args.maxedges
        uracil = args.rna
        tolerancefive = args.tolerance5prime
        tolerancethree = args.tolerance3prime
        tolerance = (tolerancefive,tolerancethree)
        track = args.track

        if outdir!='./':
            if not os.path.exists(outdir):
                os.mkdir(outdir)
     
        

        #Make directory to contain output and results for the current project run
        if not os.path.exists(outdir+"/"+project_name):
            os.mkdir(outdir+"/"+project_name)
        else:
            shutil.rmtree(outdir+"/"+project_name)
            os.mkdir(outdir+"/"+project_name)

 
        #Make a sub-directory to store run/log files ect
        run_files_path = outdir+"/"+project_name+"/Run_Files"
        if not os.path.exists(run_files_path):
            os.mkdir(run_files_path)

        #Make a sub-directory to store html output file
        html_files_path = outdir+"/"+project_name+"/Html_Files"
        if not os.path.exists(html_files_path):
            os.mkdir(html_files_path)
        

        #CALL CLASSES AND METHODS AS NEEDED 

        #First step: Parse FASTA files
        headers,sequences,seq_lengths,intron_indices,uracil = SM.parse_fasta_file(outdir,project_name,filename,uracil)
        if len(sequences)==0:
            print "\n==============================================================\nNO SEQUENCES FOUND IN INPUT FILE...\n=============================================================="
            print "\n===================================================================================================\n\nLncLOOM COMPLETE\n\n==================================================================================================="
            sys.exit()
        #Check if fasta file was parsed correctly
        if len(sequences)!=len(headers) or len(sequences)!=len(intron_indices) or len(sequences)!=len(seq_lengths):
            print "\n==============================================================\nERROR! FORMAT OF FASTA FILE NOT RECOGNISED...\n=============================================================="
            print "\n===================================================================================================\n\nLncLOOM COMPLETE\n\n==================================================================================================="
            sys.exit()
        



        #Then filter ti remove sequences with high percentage ID
        #Run a mafft alignment - only if needed (if similarity threshold<100 or iterations>0)
        mafft = False
        if similarity<100:
            mafft = SM.run_mafft_msa(outdir,project_name)

        headers,sequences,seq_lengths,intron_indices =SM.exlcudeSimilarSequences(outdir,project_name,sequences,headers,seq_lengths,intron_indices,mafft,similarity)
            
        number_of_layers = len(sequences)
        
        if min_depth >number_of_layers:
            min_depth = number_of_layers

        hsps=[[] for s in sequences]

        #Second step: Reorder sequences and calculate HSPs, as specified by user input
        if inputorder:
            #BLAST is not run, order of sequences in maintained according to fasta file input
            print "\n==============================================================\nBLAST NOT IMPLEMENTED, INPUT ORDER RETAINED...\n==============================================================\n"

        else:
            #Run BLAST
            blastFile = outdir+"/"+project_name+"/Run_Files/Sequences.fasta"
            successful_blast = LB.run_blast(outdir,project_name,blastFile,"main")
            if successful_blast:
                #Now order sequences, headers, lengths and hsps according to BLAST results
                print "\n==============================================================\nBLAST SUCCESSFUL - SEQUENCES ORDERD BY HOMOLOGY...\n"
                headers,sequences,seq_lengths,intron_indices,hsps = LB.parse_blast(outdir,project_name,blastFile,headers,sequences,seq_lengths,intron_indices,"main")


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
                #Impose HSP contraints if blast == True
                if hspblast:
                    print "HSPs CONSTRAINTS IMPOSED ON GRAPH CONSTRUCTION...\n==============================================================\n"
                else:
                    hsps=[[] for s in sequences]
                    print "NO CONSTRAINTS IMPOSED ON GRAPH CONSTRUCTION...\n==============================================================\n"
            else:
                print "\n==============================================================\nBLAST NOT SUCCESSFUL - BLAST RESULTS NOT IMPLEMENTED - INPUT ORDER RETAINED...\n" 

        #Identify consecutive layers of identical sequences:
        '''identical_sequences = [] 
        for i in range(len(sequences)-1):
            if sequences[i]==sequences[i+1]:
                identical_sequences.append((i+1,i+2))'''

        #If bedout or eclip parameters selected and user chosen to run blat: then run blat"
        track_query = 0
        blat_successful = False
        track_method = ''
        genome = ''
        track_bed = ''
        format_file = True
        if track ==True:
            b = open('./src/for_track_output.txt','r')
            blines = b.readlines()
            b.close()
            query = blines[0].split(':')
            if (query[0].strip()).upper()=="QUERY LAYER":
                if int(query[1].strip())>0 and int(query[1].strip())<=len(sequences):
                    track_query = int(query[1].strip())
                else:
                    print "ERROR! Query layer specified for track output is not in graph"
                    format_file = False
            else:
                print "ERROR! for_track_output.txt format is incorrect"
                format_file = False

            query = blines[1].split(':')
            if (query[0].strip()).upper()=="BLAT":
                if os.path.isfile(query[1].strip()):
                    track_method = "BLAT"
                    genome = query[1].strip()
                else:
                    print "ERROR! File not found: "+query[1].strip()
                    format_file = False
            elif (query[0].strip()).upper()=="BED":
                if os.path.isfile(query[1].strip()):
                    track_method = "BED"
                    track_bed = query[1].strip()
                else:
                    print "ERROR! File not found: "+query[1].strip()
                    format_file = False
            else:
                print "ERROR! for_track_output.txt format is incorrect"
                format_file = False

        psl_file = ''
        if track_method == "BLAT" and format_file:
            print "Running BLAT"
            query_header = headers[track_query-1]
            query_seq = sequences[track_query-1]
            aligned,psl_file =LA.run_blat(outdir,project_name,query_header,query_seq,genome)
            blat_successful = aligned

        
        #Get a tuple containing 0)the main LncLOOM graph and 1) a dictionary of graphs where nodes have conservation to each layer 
        Motif_Graph_Levels,MainGraph,Main_LOOM_Level,LOOM5,LOOM3,LOOM5_Levels,LOOM3_Levels,details5,details3 = SM.start_lncLOOM(sequences,seq_lengths,kmers_len,stopw,hsps,solver,prune,min_depth,shorttandem,maxedges,noconstraints,tolerance)
        Modules = SM.get_modules(Main_LOOM_Level, LOOM5_Levels, LOOM3_Levels,number_of_layers)

        #Motif_Graph_Levels = Motif_Graph_All[1]
        print "\n==============================================================\nMOTIFS DISCOVERY SUCCESSFUL\n"

        #Fourth step: Extract information from Motif_Graph
        ## - firstly, obtain dictionaries of locations and seqs of nodes/kmers per layer (sequence)
        
        kmers_dict = SM.get_kmers_per_layer(Motif_Graph_Levels,number_of_layers)



        LevelSolutions = ()
        msg=''
        select = select.strip('>')
        if select!='DEEPEST': 
            try:
                depth = headers.index('>'+select)+1
                if depth<min_depth:
                    msg= "Depth of "+select+" is greater than mindepth parameter: "+str(min_depth)+"<br><br>Specified conservation not calculated"
                else:
                    if depth in Motif_Graph_Levels:
                        all_deep_levels = [x for x in Motif_Graph_Levels if x>= depth]
                        current_graphs = {x:Motif_Graph_Levels[x] for x in all_deep_levels}
                        kmers_dict_depth = SM.get_kmers_per_layer(current_graphs,number_of_layers)
                        #Motif_Graph_Depth = Motif_Graph_Levels[depth]
                        #kmers_dict_depth = SM.get_kmers_per_layer({depth:Motif_Graph_Depth},number_of_layers)
                        LevelSolutions = (kmers_dict_depth,depth)
                    else:
                    
                        msg= "No conserved nodes found in sequence: "+select

            except ValueError:
                msg = "Name of sequence not found: Specified conservation not calculated"
                
        else:
           #Find the deepest path
           level_keys = Motif_Graph_Levels.keys()
           if level_keys:
               level_keys.sort(reverse=True)
               #lg = Motif_Graph_Levels[level_keys[0]]
               all_deep_levels = [x for x in Motif_Graph_Levels if x>= level_keys[0]]
               current_graphs = {x:Motif_Graph_Levels[x] for x in all_deep_levels}
               kmers_dict_depth = SM.get_kmers_per_layer(current_graphs,number_of_layers)
               #kmers_dict_depth = SM.get_kmers_per_layer({level_keys[0]:lg},number_of_layers)
               LevelSolutions = (kmers_dict_depth,level_keys[0])
           else:
               msg = "No conserved motifs found! Graph is empty"
        
        if msg!='':
            LO.write_error_page(outdir,project_name,msg,'kmers_in_seqs_level.html')
            LO.write_error_page(outdir,project_name,msg,'kmers_in_blocks_level.html')
            LO.write_error_page(outdir,project_name,msg,'kmers_in_seqs_level_graded.html')

        #layers = Motif_Graph_Levels.keys()
        graded_colours = LO.define_colour_grade(number_of_layers)
        kmer_graded_colours = {}
        
        #Calculate the depth of each kmer
        Depth_of_kmers = LS.get_depth_of_kmers(number_of_layers,kmers_dict)
        #Add species name to Depth_of_kmers
        for kmer in Depth_of_kmers:
            d = Depth_of_kmers[kmer]
            Depth_of_kmers[kmer] = [d,headers[d-1].strip('>')]
            kmer_graded_colours[kmer] = graded_colours[d]

        #Define colours for output
        kmer_colour_dict = LO.define_kmer_colours(outdir,project_name,kmers_dict,newcolours,uracil)

        #Annotate the nodes in the current solution
        Annotations = {i+1:{} for i in range(number_of_layers)}
        Annotations_Conservation = {i:{} for i in range(1,number_of_layers)}
        #AnnotationsLevel = {i+1:{} for i in range(number_of_layers)}
        stats_dict = {}

        #Write Homepage of results
        top = headers[0].strip('>')[0:10]
        LO.write_homepage(outdir,project_name,select,top+'.')
        

        #Annotate eCLIP: here we are annotating and creating tables for all combined and sub-kmers in the selected query layer.
        #This only needs to be done once, as the annotations for all k-mers will be stored - loop not needed here
        if eCLIP!='None':
            e = open(eCLIP,'r')
            elines = e.readlines()
            e.close()
            #Set default values
            queryID =1
            genome =None
            bedfile  = None 
            eCLIP_Paths = []
            eCLIP_Path_names = []
            #Obtain values from file:
            for line in elines:
                case = ((line.split(':')[0]).strip()).upper()
                value = (line.split(':')[1]).strip()
                if case =='BLAT':
                    if os.path.isfile(value):
                        genome = value
                    else:
                        print "ERROR! Genome file not found: "+value
                        Valid=False

                elif case == 'QUERY LAYER':
                    try:
                        queryID = int(value)
                    except ValueError:
                        print "ERROR! Invalid query sequence in eCLIP Annotation File: "+info[1].strip()+'\nDefault Layer = 1 has been used for annotation mapping'

                elif case =='ECLIP':
                    eCLIP_Path_names.append(value)
                    try:
                        e_path = (line.split(':')[2]).strip()
                        if os.path.exists(e_path):
                            eCLIP_Paths.append(e_path)
                        else:
                            print "ERROR! Path to eCLIP data not found: "+e_path
                            Valid = False
                    except IndexError:
                        print "ERROR! Invalid format in eCLIP Annotation File"
                        Valid = False
                elif case == 'BED':
                    if os.path.isfile(value):
                        bedfile = value
                    else:
                        print "ERROR! Bed File not found: "+value
                        Valid=False
 
                else:
                    print "ERROR! Format of eCLIP Annotation File Invalid"
                    Valid = False
                       
                   
            if Valid:
                Annotations[queryID]['eCLIP']={}
                #AnnotationsLevel[queryID]['eCLIP']={}
                query_seq = sequences[queryID-1]
                query_seq = query_seq.replace('U','T')
                query_header = headers[queryID-1]
                nodes_of_query = []

                
                for x in kmers_dict[queryID]:
                    for node in x[1]: 
                       nodes_of_query.append((node[0],node[1],node[2]))     
                if eCLIP_Paths and genome!=None:
                    try:
                    
                    #In this case - run blat to obtain alignment of query sequence to genome
                        if track_query==queryID and blat_successful:
                            aligned,chrm,strand,num_of_exons,size_of_exons,seq_exons,chromosome_coords = LA.parse_psl(psl_file,query_seq)
                        else:
                            print "Running BLAT"
                            aligned,psl_file = LA.run_blat(outdir,project_name,query_header,query_seq,genome)
                            aligned,chrm,strand,num_of_exons,size_of_exons,seq_exons,chromosome_coords = LA.parse_psl(psl_file,query_seq)

                        if aligned:
                            start_min = chromosome_coords[0] #start position of the first exon on the chromosome
                            start_max = chromosome_coords[-1]+size_of_exons[-1] #eCLIP data will only be extracted if it starts before the end of the last exon
                            seq_enrichment = LA.seqEnrichment(eCLIP_Paths,size_of_exons,chromosome_coords,len(query_seq),strand,chrm)
                            eCLIP_data = LA.extract_eCLIP_data(eCLIP_Paths,chrm,start_min,start_max,strand)
                            eCLIP_data = LA.format_eCLIP(eCLIP_data,seq_enrichment)
                            eCLIP_Annotations,nodes_in_chromosome = LA.map_eCLIP(num_of_exons,size_of_exons,seq_exons,chromosome_coords,nodes_of_query,eCLIP_data,strand,query_seq)
                            LO.write_eCLIP(eCLIP_Annotations,nodes_in_chromosome,outdir,project_name,'BLAT',query_header,strand,chrm,kmer_colour_dict)
                            Annotations[queryID]['eCLIP']['BLAT']=eCLIP_Annotations
                            #AnnotationsLevel[queryID]['eCLIP']['BED']=eCLIP_Annotations
                            for layer in Annotations_Conservation: 
                                Annotations_Conservation[layer]['eCLIP']={'BLAT':eCLIP_Annotations} 
                    
                    except Exception as e:
                        print "eCLIP ANNOTATION ERROR!"
                        print str(e)

                if eCLIP_Paths and bedfile:
                    try:
                        #In this case obtain alignment info from the bedfile uploaded bu user
                        parsed,chrm,strand,num_of_exons,size_of_exons,seq_exons,chromosome_coords = LA.parse_bedFile(bedfile,query_seq)
                        if parsed:
                            start_min = chromosome_coords[0] #start position of the first exon on the chromosome
                            start_max = chromosome_coords[-1]+size_of_exons[-1] #eCLIP data will only be extracted if it starts before the end of the last exon
                            seq_enrichment = LA.seqEnrichment(eCLIP_Paths,size_of_exons,chromosome_coords,len(query_seq),strand,chrm)
                            eCLIP_data = LA.extract_eCLIP_data(eCLIP_Paths,chrm,start_min,start_max,strand)
                            eCLIP_data = LA.format_eCLIP(eCLIP_data,seq_enrichment)   
                            eCLIP_Annotations,nodes_in_chromosome = LA.map_eCLIP(num_of_exons,size_of_exons,seq_exons,chromosome_coords,nodes_of_query,eCLIP_data,strand,query_seq)
                            LO.write_eCLIP(eCLIP_Annotations,nodes_in_chromosome,outdir,project_name,'BED',query_header,strand,chrm,kmer_colour_dict)
                            Annotations[queryID]['eCLIP']['BED']=eCLIP_Annotations
                            #AnnotationsLevel[queryID]['eCLIP']['BED']=eCLIP_Annotations
                            for layer in Annotations_Conservation: 
                                Annotations_Conservation[layer]['eCLIP']={'BED':eCLIP_Annotations} 
                    except Exception as e:
                        print "eCLIP ANNOTATION ERROR!"
                        print str(e)


        #Print html files for main graph
        if targetScan:
            targetScanAnnotations,miRNA_Tables = LA.annotate_TargetScan(targetScan,project_name,kmers_dict,sequences,uracil)
            LO.write_miRNA_Tables(outdir,project_name,headers,miRNA_Tables,kmer_colour_dict)
            for layer in targetScanAnnotations:
                Annotations[layer]['TargetScan']=targetScanAnnotations[layer]

        print "Writing Pre-HTML Result (prior to statistical analysis)..."
        #LO.write_kmers_in_seq_html(headers,sequences,intron_indices,kmers_dict,project_name,Annotations,stats_dict,Depth_of_kmers,kmer_colour_dict,'kmers_in_seqs.html')
        LO.write_lncLOOM_txt(kmers_dict,MainGraph,Main_LOOM_Level,number_of_layers,Depth_of_kmers,LOOM5,LOOM3,LOOM5_Levels,LOOM3_Levels,details5,details3,sequences,headers,outdir,project_name,stats_dict,Annotations)
        LO.write_kmers_in_seq_overlap_html(headers,sequences,intron_indices,kmers_dict,outdir,project_name,Annotations,stats_dict,Depth_of_kmers,kmer_colour_dict,'kmers_in_seqs.html',number_of_layers)
        LO.write_kmers_graded_html(headers,[x[0:25].strip('>') for x in headers],sequences,intron_indices,kmers_dict,outdir,project_name,Annotations,stats_dict,Depth_of_kmers,kmer_graded_colours,graded_colours,'kmers_in_seqs_graded.html',number_of_layers)
        LO.write_kmers_in_blocks_overlap_html(headers,sequences,intron_indices,kmers_dict,outdir,project_name,Annotations,stats_dict,Depth_of_kmers,kmer_colour_dict,'kmers_in_blocks.html',number_of_layers)
        #LO.write_kmers_in_seq_overlap_html(headers,sequences,intron_indices,kmers_dict,project_name,Annotations,stats_dict,Depth_of_kmers,kmer_colour_dict,'kmers_in_seqs_overlap.html')
        #LO.write_kmers_in_blocks_overlap_html(headers,sequences,intron_indices,kmers_dict,project_name,Annotations,stats_dict,Depth_of_kmers,kmer_colour_dict,'kmers_in_blocks_overlap.html')
        #Write modules

        LO.write_modules(Modules,headers,kmer_colour_dict,Annotations,outdir,project_name,kmer_graded_colours,graded_colours,stats_dict,Depth_of_kmers)

        
        #Print html files for specified level kmers
        if LevelSolutions:
            kmers_dict_level = LevelSolutions[0]
            depth = LevelSolutions[1]
            '''if targetScan:
                #targetScanAnnotations,miRNA_Tables = LA.annotate_TargetScan(targetScan,project_name,kmers_dict_level,sequences,uracil)
                for layer in AnnotationsLevel:
                    AnnotationsLevel[layer]['TargetScan']=targetScanAnnotations[layer]

            LO.write_kmers_in_seq_overlap_html(headers[0:depth],sequences[0:depth],intron_indices[0:depth],kmers_dict_level,project_name,Annotations,stats_dict,Depth_of_kmers,kmer_colour_dict,'kmers_in_seqs_level.html')
            LO.write_kmers_graded_html(headers[0:depth],[x[0:25].strip('>') for x in headers[0:depth]],sequences[0:depth],intron_indices[0:depth],kmers_dict_level,project_name,Annotations,stats_dict,Depth_of_kmers,kmer_graded_colours,graded_colours,'kmers_in_seqs_level_graded.html')
            LO.write_kmers_in_blocks_overlap_html(headers[0:depth],sequences[0:depth],intron_indices[0:depth],kmers_dict_level,project_name,Annotations,stats_dict,Depth_of_kmers,kmer_colour_dict,'kmers_in_blocks_level.html')'''


            LO.write_kmers_in_seq_overlap_html(headers,sequences,intron_indices,kmers_dict_level,outdir,project_name,Annotations,stats_dict,Depth_of_kmers,kmer_colour_dict,'kmers_in_seqs_level.html',depth)
            LO.write_kmers_graded_html(headers,[x[0:25].strip('>') for x in headers],sequences,intron_indices,kmers_dict_level,outdir,project_name,Annotations,stats_dict,Depth_of_kmers,kmer_graded_colours,graded_colours,'kmers_in_seqs_level_graded.html',depth)
            LO.write_kmers_in_blocks_overlap_html(headers,sequences,intron_indices,kmers_dict_level,outdir,project_name,Annotations,stats_dict,Depth_of_kmers,kmer_colour_dict,'kmers_in_blocks_level.html',depth)




        #Print html files for conserved kmers per layer
        topNodes = []
        for x in kmers_dict[1]:
            for node in x[1]: 
                topNodes.append(node)     
        kmers_dict_conservation = SM.kmers_by_conservation(topNodes,number_of_layers-1)

        #kmers_dict_conservation = {}
        headers_conservation = []
        sequences_conservation = []
        intron_indices_conservation = []

        for i in range(1,number_of_layers):
            #kmers_dict_conservation[i]={}
            headers_conservation.append(headers[0]+' TO '+headers[i].strip('>'))
            sequences_conservation.append(sequences[0])
            intron_indices_conservation.append(intron_indices[0])

        '''for layer in range(1,min_depth-1):
            kmers_dict_depth = SM.get_kmers_per_layer({layer:Motif_Graph_Levels[min_depth]},number_of_layers)
            kmers_dict_conservation[layer] = kmers_dict_depth[1]


        for layer in Motif_Graph_Levels:
            kmers_dict_depth = SM.get_kmers_per_layer({layer:Motif_Graph_Levels[layer]},number_of_layers)
            kmers_dict_conservation[layer-1] = kmers_dict_depth[1]'''



        if targetScan:
            #targetScanAnnotations,miRNA_Tables = LA.annotate_TargetScan(targetScan,project_name,kmers_dict_conservation,sequences_conservation,uracil) 
            for layer in Annotations_Conservation:
                Annotations_Conservation[layer]['TargetScan']=targetScanAnnotations[1]
           

        LO.write_kmers_graded_html(headers_conservation,[x[0:25].strip('>') for x in headers[1:]],sequences_conservation,intron_indices_conservation,kmers_dict_conservation,outdir,project_name,Annotations_Conservation,stats_dict,Depth_of_kmers,kmer_graded_colours,graded_colours,'kmers_in_seqs_layer_conservation.html',number_of_layers)
        LO.write_graded_blocks_html(headers_conservation,[x[0:25].strip('>') for x in headers[1:]],sequences_conservation,intron_indices_conservation,kmers_dict_conservation,outdir,project_name,Annotations_Conservation,stats_dict,Depth_of_kmers,kmer_graded_colours,graded_colours,'kmers_in_blocks_layer_conservation.html')


        #write track file
        if track == True:
            query_seq = sequences[track_query-1]
            query_head = headers[track_query-1]
            parsed = False
            if blat_successful and track_method == "BLAT":
                parsed,chrm,strand,num_of_exons,size_of_exons,seq_exons,chromosome_coords = LA.parse_psl(psl_file,query_seq)

            elif format_file and track_method == "BED":
                parsed,chrm,strand,num_of_exons,size_of_exons,seq_exons,chromosome_coords = LA.parse_bedFile(bedfile,query_seq)

            if parsed:
                print "Printing track output..."
                kmers_in_query = kmers_dict[track_query]
                #LO.print_trackRGB(kmers_in_query,kmer_graded_colours,chrm,strand,num_of_exons,size_of_exons,seq_exons,chromosome_coords,outdir,project_name,query_head,query_seq )
                LO.print_trackScore(kmers_in_query,kmer_graded_colours,chrm,strand,num_of_exons,size_of_exons,seq_exons,chromosome_coords,outdir,project_name,query_head,headers,Depth_of_kmers,query_seq )

            else:
                print "ERROR! Track output not successful (psl or bedfile not recognised)"












        #Run stats
        if iterations>0:
            mafft = SM.run_mafft_msa(outdir,project_name)
            print "Starting Stats Analysis..."
            p1,p2 = LS.run_stats_msa(iterations,number_of_layers,kmers_dict,outdir,project_name,kmers_len,stopw,seq_lengths,prune,solver,min_depth,shorttandem,maxedges,noconstraints,multiprocess,mafft,hspblast,tolerance)
            print "Stats: Percentage ID Iterations Complete"
            #Random generated sequence sets with same dinucleotide frequencies
            p3,p4 = LS.run_stats_random(iterations,number_of_layers,kmers_dict,outdir,project_name,kmers_len,stopw,sequences,seq_lengths,prune,solver,min_depth,shorttandem,maxedges,noconstraints,multiprocess,hspblast,tolerance)
            print "Stats: Random Iterations Complete"
            for kmer in p1:
                stats_dict[kmer] = {1:p1[kmer],2:p2[kmer],3:p3[kmer],4:p4[kmer]}


            print "Printing Final HTML Results (Stats included)..."
            #LO.write_kmers_in_seq_html(headers,sequences,intron_indices,kmers_dict,project_name,Annotations,stats_dict,Depth_of_kmers,kmer_colour_dict,'kmers_in_seqs.html')
            LO.write_kmers_in_seq_overlap_html(headers,sequences,intron_indices,kmers_dict,outdir,project_name,Annotations,stats_dict,Depth_of_kmers,kmer_colour_dict,'kmers_in_seqs.html',number_of_layers)
            LO.write_kmers_graded_html(headers,[x[0:25].strip('>') for x in headers],sequences,intron_indices,kmers_dict,outdir,project_name,Annotations,stats_dict,Depth_of_kmers,kmer_graded_colours,graded_colours,'kmers_in_seqs_graded.html',number_of_layers)
            LO.write_kmers_in_blocks_overlap_html(headers,sequences,intron_indices,kmers_dict,outdir,project_name,Annotations,stats_dict,Depth_of_kmers,kmer_colour_dict,'kmers_in_blocks.html',number_of_layers)
            #LO.write_kmers_in_seq_overlap_html(headers,sequences,intron_indices,kmers_dict,project_name,Annotations,stats_dict,Depth_of_kmers,kmer_colour_dict,'kmers_in_seqs_overlap.html')
            #LO.write_kmers_in_blocks_overlap_html(headers,sequences,intron_indices,kmers_dict,project_name,Annotations,stats_dict,Depth_of_kmers,kmer_colour_dict,'kmers_in_blocks_overlap.html')
            #Write modules
            LO.write_modules(Modules,headers,kmer_colour_dict,Annotations,outdir,project_name,kmer_graded_colours,graded_colours,stats_dict,Depth_of_kmers)

            #Print html files for specified level kmers
            if LevelSolutions:
                LO.write_kmers_in_seq_overlap_html(headers,sequences,intron_indices,kmers_dict_level,outdir,project_name,Annotations,stats_dict,Depth_of_kmers,kmer_colour_dict,'kmers_in_seqs_level.html',depth)
                LO.write_kmers_graded_html(headers,[x[0:25].strip('>') for x in headers],sequences,intron_indices,kmers_dict_level,outdir,project_name,Annotations,stats_dict,Depth_of_kmers,kmer_graded_colours,graded_colours,'kmers_in_seqs_level_graded.html',depth)
                LO.write_kmers_in_blocks_overlap_html(headers,sequences,intron_indices,kmers_dict_level,outdir,project_name,Annotations,stats_dict,Depth_of_kmers,kmer_colour_dict,'kmers_in_blocks_level.html',depth)

            #Print html files for conserved kmers per layer
            LO.write_kmers_graded_html(headers_conservation,[x[0:25].strip('>') for x in headers[1:]],sequences_conservation,intron_indices_conservation,kmers_dict_conservation,outdir,project_name,Annotations_Conservation,stats_dict,Depth_of_kmers,kmer_graded_colours,graded_colours,'kmers_in_seqs_layer_conservation.html',number_of_layers)
            LO.write_graded_blocks_html(headers_conservation,[x[0:25].strip('>') for x in headers[1:]],sequences_conservation,intron_indices_conservation,kmers_dict_conservation,outdir,project_name,Annotations_Conservation,stats_dict,Depth_of_kmers,kmer_graded_colours,graded_colours,'kmers_in_blocks_layer_conservation.html')

            #LO.write_lncLOOM_txt(Main_LOOM_Level,number_of_layers,Depth_of_kmers,LOOM5_Levels,LOOM3_Levels,details5,details3,sequences,headers,project_name,stats_dict,Annotations)
            LO.write_lncLOOM_txt(kmers_dict,MainGraph,Main_LOOM_Level,number_of_layers,Depth_of_kmers,LOOM5,LOOM3,LOOM5_Levels,LOOM3_Levels,details5,details3,sequences,headers,outdir,project_name,stats_dict,Annotations)
    else:
        print "\n==============================================================\nINPUT ARGUMENTS ARE INVALID. EXITING PROGRAM...\n=============================================================="
        print "\n===================================================================================================\n\nLncLOOM COMPLETE\n\n==================================================================================================="
        sys.exit()
    
    print "\n===================================================================================================\n\nLncLOOM COMPLETE\n\n==================================================================================================="
    


if __name__ == '__main__':
    main()

