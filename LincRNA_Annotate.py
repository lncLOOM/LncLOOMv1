import subprocess
import os
import sys
import re
import pyBigWig
from operator import itemgetter


def get_reverse_complement(RNA,uracil):
    reverse_complement=''
    if uracil:
        ComplimentBases = {'A':'U','C':'G','G':'C','U':'A'}
        for b in RNA[::-1]:
            reverse_complement+=ComplimentBases[b]
    else:
        ComplimentBases = {'A':'T','C':'G','G':'C','U':'A'}
        for b in RNA[::-1]:
            reverse_complement+=ComplimentBases[b]

    return reverse_complement


def parse_mir_species_table():

    SpeciesFile = 'src/TargetScan/species100.html'
    Species = {}
    #Get ID of species
    try:
        f = open(SpeciesFile,'r')
        speciesHTML = f.readlines()
        f.close()
        start = speciesHTML.index("<TR><TD><B>Symbol</TD><TD><B>Common name</TD><TD><B>Species</TD><TD><B>Species ID</TD></TR>\n")+1
        end = speciesHTML.index("</TABLE>\n")
        all_species = speciesHTML[start:end]
        
        for spec in all_species:
            tabs = spec.split('<TD>')[1:]
            
            ID = tabs[3].rstrip('</TD></TR>\n')
            Name = tabs[1].rstrip('</TD>').strip()
            ScientificName = tabs[2].rstrip('</TD>')
            Species[ID] = Name+" ("+ScientificName+")"
    except:
        print "ERROR! Species file is not compatable"

    return Species





def annotate_TargetScan(filename,project_name,combined_kmers_dict,sequences,uracil):

    #For now TargetScan annotations are done according to src files that are set - user can update these files if needed???
    #TargetScan annotations are performed if user sets --targetScan True (default = False)

    TargetScanMatches={}
    miRNA_Tables = {}
    
    FamilyFile = 'src/TargetScan/miR_Family_Info.txt'
    Species = parse_mir_species_table()

    SeedFamily = {}
    SeedConservation = {}
    SeedSpecies = {}
    SeedMatches = {}
    Conservation = {'2':'Broadly Conserved','1':'Conserved','0':'Poorly Conserved','-1':'Poorly Conserved'}   


    #Parses the text file and extracts summary details for each miRNA family     
    try:
        #Read in miR-Family according to seed
        f = open(FamilyFile,'r')
        mirFamily = f.readlines()[1:]
        f.close()
        #Remove last blankline
        if mirFamily[-1]=='\n':
            mirFamily = mirFamily[0:-1]
        for record in mirFamily:
            info = record.split()
            seed = info[1].strip()
            family = info[0].strip()
            species = info[2].strip()
            conservation = info[5].strip()
            #if conservation =='2' or conservation =='1':
            if conservation >='1' and species=='9606':
                if not seed in SeedFamily:
                    SeedFamily[seed] = family
                    SeedConservation[seed] = conservation
                    SeedSpecies[seed] = [species]

        for record in mirFamily[0:-1]:
            info = record.split()
            seed = info[1].strip()
            family = info[0].strip()
            species = info[2].strip()
            if seed in SeedFamily:
                if family == SeedFamily[seed]:
                    if species not in SeedSpecies[seed]:
                        SeedSpecies[seed].append(species)

    except:
        print "ERROR! miRNA Family file is not compatable"


    #Compare node kmers to seeds of families
    #All seeds are 7-kmers in lenght (6mer+m8 site)
    

    for layer in combined_kmers_dict:
        TargetScanMatches[layer] = {}
        combined_nodes_of_layer = []
        for x in combined_kmers_dict[layer]:
            for km in x[1]:
                combined_nodes_of_layer.append((km[0],km[1],km[2]))
        #combined_nodes_of_layer = [x[0] for x in combined_kmers_dict[layer]]
        sequence = sequences[layer-1]
        for node in combined_nodes_of_layer:
            matches = []
            kmer = node[2]
            klen = len(kmer)
            kmer_7M8 = kmer
            kmer_7A1 = kmer
            kmer_8 = kmer
            

            start = node[0]
            if start!=1:
                kmer_7M8 = sequence[start-2]+kmer_7M8

            if start+klen-1<len(sequence):
                kmer_7A1 = kmer_7A1+sequence[start+klen-1]


            if start>1 and start+klen-1<len(sequence):
                kmer_8 = sequence[start-2]+kmer_8+sequence[start+klen-1]


            for seed in SeedFamily:
                family = SeedFamily[seed]
                M7mer_m8 = get_reverse_complement(seed,uracil)
                M8mer = M7mer_m8+'A'
                M6mer = M7mer_m8[1:]
                M7mer_A1 = M6mer+'A'
                
                if M8mer in kmer_8:
                    try:
                        SeedMatches[seed].append((layer,kmer,'8mer'))
                    except KeyError:
                        SeedMatches[seed]=[(layer,kmer,'8mer')]
                    matches.append(family+':'+seed)

                elif M7mer_m8 in kmer_7M8:
                    try:
                        SeedMatches[seed].append((layer,kmer,'7mer-m8'))
                    except KeyError:
                        SeedMatches[seed]=[(layer,kmer,'7mer-m8')]
                    matches.append(family+':'+seed)

                elif M7mer_A1 in kmer_7A1:
                    try:
                        SeedMatches[seed].append((layer,kmer,'7mer-A1'))
                    except KeyError:
                        SeedMatches[seed]=[(layer,kmer,'7mer-A1')]
                    matches.append(family+':'+seed)
                elif M6mer in kmer:
                    try:
                        SeedMatches[seed].append((layer,kmer,'6mer'))
                    except KeyError:
                        SeedMatches[seed]=[(layer,kmer,'6mer')]
                    matches.append(family+':'+seed)
                       
           
            TargetScanMatches[layer][node]=matches

    #Get a dictionary that defines info for miRNA families that match the lincRNA motifs - this is used later in output
   
    for seed in SeedMatches:
        
        conservation = SeedConservation[seed]
        conservation = Conservation[conservation]
        species = SeedSpecies[seed]
        for i,s in enumerate(species):
            try:
                species[i] = Species[s]
            except KeyError:
                species[i] = 'Unknown Species ID'
        matches = SeedMatches[seed]
        family = SeedFamily[seed]
        table_info = [family,conservation,species,matches]
        miRNA_Tables[seed] = table_info
               
    return TargetScanMatches,miRNA_Tables



def run_blat(outdir,project_name,default_head,default_seq,genome):
    num_of_exons = 0
    exon_sizes = []
    exon_pos_seq = []
    exon_pos_chrom = []
    chrm = ''
    strand = ''
    aligned = True
    try:
        default_head=default_head.strip('>')

        #create fasta file with single sequence
        w= open(outdir+'/'+project_name+'/Run_Files/query_'+default_head+'.fa','w',0) #set buffereing to false so that file is written immediately 
        w.write('>'+default_head+'\n'+default_seq)
        w.close
    
        #delay code until query.fa exists
        while not os.path.exists(outdir+'/'+project_name+'/Run_Files/query_'+default_head+'.fa'):
            time.sleep(1)


        subprocess.call(["./src/blat",genome,outdir+'/'+project_name+'/Run_Files/query_'+default_head+'.fa', outdir+'/'+project_name+'/Run_Files/'+default_head+".psl"])
        while not os.path.exists(outdir+'/'+project_name+'/Run_Files/'+default_head+".psl"):
            time.sleep(1)

    except Exception as e:
        print "\nERROR! BLAT UNSUCCESSFUL!\n"
        print str(e)+'\n'
        aligned = False

    return aligned,outdir+'/'+project_name+'/Run_Files/'+default_head+".psl"



def parse_psl(psl_file,default_seq):

    num_of_exons = 0
    exon_sizes = []
    exon_pos_seq = []
    exon_pos_chrom = []
    chrm = ''
    strand = ''
    aligned = True

    try:

        #Parse psl output from blat and obtain exon and chromosome positions
        f = open(psl_file,'r')
        psl_output = f.readlines()
        f.close
    
        
        #Check if a match was found:
        if len(psl_output)==5:
            print "SEQUENCE HAD NO SIGNIFICANT ALIGNMENT TO GENOME"
            aligned = False
        else: #get the output for the top hit
            hits = psl_output[5:]
            for h in range(len(hits)):
                hits[h] = hits[h].split()
                hits[h][0]=int(hits[h][0])
            hits = sorted(hits,key=itemgetter(0))
            hit_found = False
            i = len(hits)-1
            #find the top hit where exon sizes match the query sequence length
            while not hit_found and i>=0:
                top_hit = hits[i]
                strand = top_hit[8].strip() #If given a reverse compliment BLAT returns co-ords matches relative to sense strand of chromosome
                chrm = top_hit[13].strip()
                num_of_exons = int(top_hit[17].strip())
                exon_sizes = map(int,top_hit[18].strip().strip(',').split(','))
                exon_pos_seq = map(int,top_hit[19].strip().strip(',').split(','))
                exon_pos_chrom = map(int,top_hit[20].strip().strip(',\n').split(','))
                len_exons = sum(exon_sizes)
                if len_exons==len(default_seq):
                    hit_found=True
                i-=1
            if not hit_found:
                print "ERROR!!! No matches from BLAT obtained...\n"
                aligned = False

    except Exception as e:
        print "\nERROR! BLAT UNSUCCESSFUL!\n"
        print str(e)+'\n'
        aligned = False

    return aligned,chrm,strand,num_of_exons,exon_sizes,exon_pos_seq,exon_pos_chrom


def parse_bedFile(bedfile,sequence):
    num_of_exons = 0
    exon_sizes = []
    exon_pos_seq = []
    exon_pos_chrom = []
    chrm = ''
    strand = ''
    parsed = True

    try:
        f = open(bedfile,'r')
        bed_lines = f.readlines()
        f.close()
        exon_sizes = []
        #loop through lines until the first line containing data is found (sometimes the files can have headers ect)
        regex = re.compile('chr\d')#regular expression to check if first line contains alignment data
        for line in bed_lines:
            top_hit = line.split()
            if regex.match(top_hit[0]):
                break
            
        chrm = top_hit[0].strip()
        strand = top_hit[5].strip()
        start_on_chromosome = int(top_hit[1])
        num_of_exons = top_hit[9]
        exon_sizes = map(int,top_hit[10].strip().strip(',').split(','))
        relative_exon_pos = map(int,top_hit[11].strip().strip(',\n').split(','))
        exon_pos_seq = []
        exon_pos_chrom = []
        start_in_seq = 0
        for i,exon in enumerate(relative_exon_pos):
            exon_pos_chrom.append(exon+start_on_chromosome)
            #To calcualte position of exons in query sequence you must assume that the seq consists only of consecutive exons - no alignment info is given in bedfile
            exon_pos_seq.append(start_in_seq)
            start_in_seq+=exon_sizes[i]
    except:
        print "\nERROR! BED FILE FORMAT IS NOT COMPATABLE\nRECOMMENDATIONS:\n1)CHECK COLUMN SPACING IN FILE\n2)RUN WITH BLAT!!!"
        parsed = False

  
    len_exons = sum(exon_sizes)
    if len(sequence)!=len_exons:
        print "\nERROR! BED FILE IS NOT COMPATABLE WITH LENGTH OF QUERY SEQUENCE\nRECOMMENDATIONS:\n1)RUN WITH BLAT!!!\n"
        parsed = False
        num_of_exons = 0
        exon_sizes = []
        exon_pos_seq = []
        exon_pos_chrom = []
        chrm = ''
        strand = ''

    
    return parsed,chrm,strand,num_of_exons,exon_sizes,exon_pos_seq,exon_pos_chrom


def pyBigWig_eCLIP_data(eCLIP_Paths,chrm,start_min,start_max,strand):

    extracted_eCLIP = []
    for path in eCLIP_Paths:
        big_bed_files = [f for f in os.listdir(path) if f.endswith('.bb')]
    
        for bbf in big_bed_files:
            bb = pyBigWig.open(path+'/'+bbf)
            extract = bb.entries(chrm, start_min, start_max)
            bb.close()
            if extract:
                for e in extract:
                    strand_e = e[2].split('\t')[2]
                    if strand_e == strand:
                        extracted_eCLIP.append(e)
             
    return extracted_eCLIP

def seqEnrichment(eCLIP_paths,size_of_exons,chromosome_coords,length,strand,chrm):

    protein_enrichment ={}
    exon_eCLIP_dict ={}
    for i in range(len(size_of_exons)):
        exon_start = chromosome_coords[i]
        exon_end = exon_start+size_of_exons[i]
        eCLIP_data= extract_eCLIP_data(eCLIP_paths,chrm,exon_start,exon_end,strand)
        exon_eCLIP_dict[i+1] = eCLIP_data
        for eCLIP in eCLIP_data:
            protein_enrichment[eCLIP[3]] = 0.00

    for i in range(len(size_of_exons)):
        exon_eCLIP = exon_eCLIP_dict[i+1]
        exon_start = chromosome_coords[i]
        exon_end = exon_start+size_of_exons[i]
        exon_eCLIP.sort(key=lambda x: x[1])

        protein_exon = {x[3]:[] for x in exon_eCLIP}
        for eCLIP in exon_eCLIP:
            protein_exon[eCLIP[3]].append(eCLIP)

        #print protein_exon

        
        for protein in protein_exon:
            peclip = protein_exon[protein]
            peclip.sort(key=lambda x:x[1])
            if peclip[-1][2]>exon_end:
                peclip[-1][2]=exon_end


            j = 0
            BoundBases=peclip[j][2]-peclip[j][1]

            for j in range(1,len(peclip)):
                si = peclip[j][1]
                se = peclip[j][2]
                
                if si<peclip[j-1][2] and se>peclip[j-1][2]:
                    BoundBases+= se-peclip[j-1][2]
                elif si>=peclip[j-1][2]:
                    BoundBases+=se-si

            protein_enrichment[protein]+=BoundBases

    for protein in protein_enrichment:
        protein_enrichment[protein] = round((protein_enrichment[protein]/float(length))*100,2)

    return protein_enrichment
    

def format_eCLIP(eCLIP_data,protein_enrichment):
    formatted = []
    for eCLIP in eCLIP_data:
        if eCLIP[3] in protein_enrichment:
            enrich = protein_enrichment[eCLIP[3]]
            if enrich >10.00:
                eCLIP[3] = eCLIP[3].lower()
            else:
                eCLIP[3] = eCLIP[3].upper()
            eCLIP.append(enrich)
            formatted.append(eCLIP)

    return formatted  

def extract_eCLIP_data(eCLIP_paths,chrm,start_min,start_max,strand):
    all_eCLIP = []
    extracted_eCLIP = pyBigWig_eCLIP_data(eCLIP_paths,chrm,start_min,start_max,strand)

    extracted_eCLIP.sort()
    enrichmentCutOff = 2.0
    significanceCutOff = 2.0
 
    for e in extracted_eCLIP:#formats the extract data to make mapping and output easier later

        start = e[0]
        end = e[1]
        annotations = e[2].split('\t')
        protein = annotations[0].split('_')[0]
        cell_line = annotations[0].split('_')[1]
        strand_e = annotations[2]
        enrichment = float(annotations[3].strip())
        significance = float(annotations[4].strip())
        if enrichment>=enrichmentCutOff and significance>= significanceCutOff:
            hit = [chrm,start,end,protein,cell_line,strand_e,enrichment,significance] 
            all_eCLIP.append(hit)
    #Sorting eCLIP data by enrichment
    #all_eCLIP = sorted(all_eCLIP,key=itemgetter(6),reverse=True)
    #Sorting eCLIP data alphabetically
    all_eCLIP = sorted(all_eCLIP,key=itemgetter(3))
    return all_eCLIP 


def map_eCLIP(num_of_exons,size_of_exons,seq_exons,chromosome_coords,nodes,all_eCLIP,strand,sequence):
   
    eCLIP_Annotations = {}

    #num_of_exons: holds the total number of exons in seq
    #size_of exons: holds the length of each exon
    #seq_exons: holds the start postion of each exon relative to the query sequence - index 0 based
    #chromosome_coords: holds the start position of each exon relative to the chromosome - index 0 based
    #nodes: holds an ordered list of tuples for each node (start,kmer) - index here is 1 based
    #all_eCLIP holds a list of lists with info on eCLIP data [start, end, binding protein, cell line]
    
    #Map eCLIP data to node postions relative to the chromosome

    #1) Obtain positions of each node relative to chromosome co-ordinates (remember nodes are index 1 based)

    #Determine which nodes are located in which exon and if any nodes span two exons
  
    nodes_in_chromosome = {}

    exon_ends_in_seq =[]
    for i,start in enumerate(seq_exons):
        exon_ends_in_seq.append(start+size_of_exons[i]) # Index for exon ends in sequence is 1 based


    nodes_starts = [] #holds start of nodes - which can be reversed if neg strand
    #If bedfile/blat indicate the lncRNA is on - strand, reverse the position of each node, so the last node (at 3' end) is the first node, appearing first in chrm relative to reference genome
    if strand =='-':
        len_seq = len(sequence)
        for node in nodes:
            node_end = node[1] #end of node is 1 indexed
            #len of seq is 1 indexed, node end is 1 indexed len_seq-node_end = 0 indexed start of node in opposite direction  
            nodes_starts.append(len_seq-node_end) #The start of the node must be presented counting from start of sequence at 3' end. So the start of each node is actually the end of each node
    elif strand == '+':
        for node in nodes:
            nodes_starts.append(node[0]-1)
    



    for n,node in enumerate(nodes): #nodes are ordered from start of sequence to end
        node_start = nodes_starts[n] #node_start is 0 indexed: same as eclip and blat
        node_end = node_start+len(node[2]) #node_end is 1 indexed: same as eclip and blat
        
        i=0
        still_to_be_mapped = True       
        while i<num_of_exons and still_to_be_mapped:
            exon_seq_start = seq_exons[i]
            exon_seq_end = exon_ends_in_seq[i]
            
            if node_start>=exon_seq_start and node_end<=exon_seq_end:
                #In this case the node is contained in a single exon (the current exon)
                chromosome_start_exon = chromosome_coords[i]
                relative_node_start_exon = node_start-exon_seq_start
                node_chrom_start = chromosome_start_exon+relative_node_start_exon
                node_chrom_end = node_chrom_start+len(node[2])
                if strand =='+':
                    exon = i+1
                elif strand =='-':
                    exon = num_of_exons-i
                nodes_in_chromosome[node]=[(exon,node_chrom_start,node_chrom_end)] #here the key of nodes_in_chromosome holds original start in lncRNA seq (5-3 orientation)
          
                #once chromosome positions of node is obtained, break from inner for loop for efficiency
                still_to_be_mapped = False

            elif node_start<exon_seq_end and node_end>exon_seq_end:
                
                #In this case the node spans over into the next exon
                #Determine how much of node is in the current exon 
                len_node_exon1= exon_seq_end-node_start
                
                #Chromosome co-ords of segment in this exon 
                chromosome_start_exon = chromosome_coords[i]
                relative_node_start_exon = node_start-exon_seq_start
                node_chrom_start = chromosome_start_exon+relative_node_start_exon
                node_chrom_end = node_chrom_start+len_node_exon1
                if strand =='+':
                    exon = i+1
                elif strand =='-':
                    exon = num_of_exons-i
                nodes_in_chromosome[node]=[(exon,node_chrom_start,node_chrom_end)]

                #loop through the next exons until node no longer maps
                stop = False
                while i+1<num_of_exons and stop==False:
                    #Determine how much of node is in next exon
                    len_node_exon2= len(node[2])-len_node_exon1
                    #Chromosome co-ords of segment in next exon   
                    node_chrom_start = chromosome_coords[i+1]
                    if node_end>exon_ends_in_seq[i+1]:
                        node_chrom_end = node_chrom_start+size_of_exons[i+1]
                        len_node_exon1+= size_of_exons[i+1]
                    else:
                        node_chrom_end = node_chrom_start+len_node_exon2
                        stop = True
                    if strand =='+':
                        exon = i+2
                    elif strand =='-':
                        exon = num_of_exons-i-1
                    nodes_in_chromosome[node].append((exon,node_chrom_start,node_chrom_end))
                    i+=1
                    

                #once chromosome positions of node is obtained, break from inner for loop for efficiency
                still_to_be_mapped = False
            else:
                nodes_in_chromosome[node]=[]
                i+=1


                
    #2) Determine overlaps between nodes (based on chromosome positions) and eCLIP data
    for node in nodes_in_chromosome:
        chromosome_pos = nodes_in_chromosome[node]
        eCLIP_matches = []
        for eCLIP in all_eCLIP:
            e_start = eCLIP[1]
            e_end = eCLIP[2]
            for segment in chromosome_pos: #loop through multiple exon spans
                seg_start = segment[1]
                seg_end = segment[2]
                if seg_start>=e_start and seg_end<=e_end:
                    #Complete overlap
                    eCLIP_matches.append((eCLIP,1)) #boolean added to mark partial overlap(0) or complete overlap(1) to eCLIP data
                elif e_start>=seg_start and e_end<=seg_end:
                    #Complete overlap
                    eCLIP_matches.append((eCLIP,1))
                elif e_end>=seg_start and e_end<seg_end:
                    eCLIP_matches.append((eCLIP,0))
                elif seg_end>=e_start and seg_start<e_start:
                    eCLIP_matches.append((eCLIP,0))
                
        eCLIP_Annotations[node] = eCLIP_matches
        
    return eCLIP_Annotations,nodes_in_chromosome
    



                





