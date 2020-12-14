import subprocess
import random
import LincMotif as SM
from multiprocessing import Pool
import LincRNA_Blast as LB


###############################
#Statistics

#Set global variables with default values
number_of_layers = 0
kmers_dict = {}
project_name = 'LincMotif'
outdir = './'
kmers_len = 15
stopw = 6
aligned_seqs = []
sequences = []
seq_lengths = []
prune = 8
solver = 'CBC'
min_depth = 2
shorttandem = False
maxedges = 1200
noconstraints = False
hspblast = False
tolerance = 0.2



def run_mafft_msa():
    mafft = True
    in_file = outdir+"/"+project_name+"/Run_Files/Sequences.fasta"
    out_file = open(outdir+"/"+project_name+"/Run_Files/Mafft.fasta","w",0)
    try:
        subprocess.call(["mafft","--maxiterate","1000","--inputorder","--quiet","--6merpair",in_file],stdout=out_file)        
    except:
        mafft = False
    out_file.close()
    return mafft


def extract_msa():
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


def remove_gaps(shuffled_msa,mp,iteration):
    NoGaps = []
    for sequence in shuffled_msa:
        seq = ''
        for base in sequence:
            if base!='-':
                seq+=base
        NoGaps.append(seq)

    w = open(outdir+"/"+project_name+"/Run_Files/ShuffledMSA_"+str(mp)+'_'+str(iteration+1)+".fasta",'w',0)
    for i in range(len(NoGaps)):
        w.write('>'+str(i+1)+'\n')
        w.write(NoGaps[i]+'\n')
    w.close()

    return NoGaps


def shuffle_alignment():
    shuffled_msa = list(aligned_seqs)
    msa_length = len(shuffled_msa[0])
    random_indexes = []
    while len(random_indexes)<msa_length*10:
        index = random.randint(0,msa_length-1)
        random_indexes.append(index)
    for i in range(0,len(random_indexes)-1,2):
        j = random_indexes[i]
        insert = random_indexes[i+1]
        for s,sequence in enumerate(shuffled_msa):
            if insert==j:
                continue
            elif insert<j:
                
                seq = sequence[0:insert]+sequence[j]+sequence[insert+1:j]+sequence[insert]+sequence[j+1:]
                shuffled_msa[s] = seq
            else:
                seq = sequence[0:j]+sequence[insert]+sequence[j+1:insert]+sequence[j]+sequence[insert+1:]
                shuffled_msa[s] = seq
    return shuffled_msa




'''def shuffle_alignment(aligned_seqs,kmers_len):
    shuffled_msa = list(aligned_seqs)
    msa_length = len(shuffled_msa[0])
    random_indexes = []
    while len(random_indexes)<msa_length*100:
        index = random.randint(0,msa_length-1)
        random_indexes.append(index)
    for i in range(0,len(random_indexes)-1,2):
        j = random_indexes[i]
        insert = random_indexes[i+1]
        for s,sequence in enumerate(shuffled_msa):
            extract = sequence[j:j+kmers_len]
            if insert==j:
                continue
            elif insert<j:
                seq = sequence[0:insert]+extract+sequence[insert:j]+sequence[j+kmers_len:]
                shuffled_msa[s] = seq
            elif insert>j+kmers_len:
                seq = sequence[0:j]+sequence[j+kmers_len:insert]+extract+sequence[insert:]
                shuffled_msa[s] = seq
            else:
                seq = sequence[0:j]+sequence[j+kmers_len:insert]+extract+sequence[insert+(j+kmers_len-insert):]
                shuffled_msa[s] = seq
 

    return shuffled_msa'''
            
def build_random_sequence_set(mp,iteration):
    #calculate nucleotide frequencies of first sequence
    nucleotides = {'A':0,'T':0,'C':0,'G':0}
    dinucleotides = {'AA':0,'AT':0,'AC':0,'AG':0,'CA':0,'CT':0,'CC':0,'CG':0,'GA':0,'GT':0,'GC':0,'GG':0,'TA':0,'TT':0,'TC':0,'TG':0}

    random_sequences = []         
    eligable_sequences = [seq for seq in sequences if len(seq)!=0]
    for seq in eligable_sequences:
        nucleotides = {'A':0,'T':0,'C':0,'G':0}
        dinucleotides = {'AA':0,'AT':0,'AC':0,'AG':0,'CA':0,'CT':0,'CC':0,'CG':0,'GA':0,'GT':0,'GC':0,'GG':0,'TA':0,'TT':0,'TC':0,'TG':0}
        seq_length=len(seq)
        for i in range(seq_length-1):
            base1 = seq[i]
            base2 = seq[i+1]
            if base1 in ['A','T','C','G'] and base2 in ['A','T','C','G'] :
                nucleotides[base1]+=1
                nucleotides[base2]+=1
                dinucleotides[base1+base2]+=1
            elif base1 =='N' and base2 =='N':
                for n in nucleotides:
                    nucleotides[n]+=2
                for dn in dinucleotides:
                    dinucleotides[dn]+=1
            elif base1 =='N' and base2 in ['A','T','C','G']:
                for n in nucleotides:
                    nucleotides[n]+=1
                    dinucleotides[n+base2]+=1
                    dinucleotides[n+base2]+=1
                    dinucleotides[n+base2]+=1
                    dinucleotides[n+base2]+=1
                nucleotides[base2]+=1
            elif base2 =='N' and base1 in ['A','T','C','G']:
                for n in nucleotides:
                    nucleotides[n]+=1
                    dinucleotides[base1+n]+=1
                    dinucleotides[base1+n]+=1
                    dinucleotides[base1+n]+=1
                    dinucleotides[base1+n]+=1
                nucleotides[base1]+=1

        total_bases = float(sum(nucleotides.values()))
        total_dibases = float(sum(dinucleotides.values()))
        if total_bases!=0:
            for n in nucleotides:
                nucleotides[n]=(nucleotides[n]/total_bases)*100
        if total_dibases !=0:
            for dn in dinucleotides:
                dinucleotides[dn]=dinucleotides[dn]/total_dibases

        transition_matrix = {'A':{'A':0,'T':0,'C':0,'G':0},'C':{'A':0,'T':0,'C':0,'G':0},'G':{'A':0,'T':0,'C':0,'G':0},'T':{'A':0,'T':0,'C':0,'G':0}}
        if (dinucleotides['AA']+dinucleotides['AC']+dinucleotides['AG']+dinucleotides['AT'])!=0:

            transition_matrix['A']['A'] = (dinucleotides['AA']/(dinucleotides['AA']+dinucleotides['AC']+dinucleotides['AG']+dinucleotides['AT']))*100
            transition_matrix['A']['T'] = (dinucleotides['AT']/(dinucleotides['AA']+dinucleotides['AC']+dinucleotides['AG']+dinucleotides['AT']))*100
            transition_matrix['A']['C'] = (dinucleotides['AC']/(dinucleotides['AA']+dinucleotides['AC']+dinucleotides['AG']+dinucleotides['AT']))*100
            transition_matrix['A']['G'] = (dinucleotides['AG']/(dinucleotides['AA']+dinucleotides['AC']+dinucleotides['AG']+dinucleotides['AT']))*100

        if (dinucleotides['TA']+dinucleotides['TC']+dinucleotides['TG']+dinucleotides['TT'])!=0:
            transition_matrix['T']['A'] = (dinucleotides['TA']/(dinucleotides['TA']+dinucleotides['TC']+dinucleotides['TG']+dinucleotides['TT']))*100      
            transition_matrix['T']['C'] = (dinucleotides['TC']/(dinucleotides['TA']+dinucleotides['TC']+dinucleotides['TG']+dinucleotides['TT']))*100 
            transition_matrix['T']['G'] = (dinucleotides['TG']/(dinucleotides['TA']+dinucleotides['TC']+dinucleotides['TG']+dinucleotides['TT']))*100
            transition_matrix['T']['T'] = (dinucleotides['TT']/(dinucleotides['TA']+dinucleotides['TC']+dinucleotides['TG']+dinucleotides['TT']))*100

        if (dinucleotides['CA']+dinucleotides['CC']+dinucleotides['CG']+dinucleotides['CT'])!=0:
            transition_matrix['C']['A'] = (dinucleotides['CA']/(dinucleotides['CA']+dinucleotides['CC']+dinucleotides['CG']+dinucleotides['CT']))*100       
            transition_matrix['C']['T'] = (dinucleotides['CT']/(dinucleotides['CA']+dinucleotides['CC']+dinucleotides['CG']+dinucleotides['CT']))*100 
            transition_matrix['C']['C'] = (dinucleotides['CC']/(dinucleotides['CA']+dinucleotides['CC']+dinucleotides['CG']+dinucleotides['CT']))*100 
            transition_matrix['C']['G'] = (dinucleotides['CG']/(dinucleotides['CA']+dinucleotides['CC']+dinucleotides['CG']+dinucleotides['CT']))*100

        if (dinucleotides['GA']+dinucleotides['GC']+dinucleotides['GG']+dinucleotides['GT'])!=0:
            transition_matrix['G']['A'] = (dinucleotides['GA']/(dinucleotides['GA']+dinucleotides['GC']+dinucleotides['GG']+dinucleotides['GT']))*100  
            transition_matrix['G']['T'] = (dinucleotides['GT']/(dinucleotides['GA']+dinucleotides['GC']+dinucleotides['GG']+dinucleotides['GT']))*100 
            transition_matrix['G']['C'] = (dinucleotides['GC']/(dinucleotides['GA']+dinucleotides['GC']+dinucleotides['GG']+dinucleotides['GT']))*100 
            transition_matrix['G']['G'] = (dinucleotides['GG']/(dinucleotides['GA']+dinucleotides['GC']+dinucleotides['GG']+dinucleotides['GT']))*100 

        #Generate strings for random base selection, according to nucleotide frequencies
        nucleotides_choice = ''
        for n in nucleotides:
            nucleotides_choice+=n*int(round(nucleotides[n]))

        if nucleotides_choice=='': #extreme case if sequence has no reconizable bases
            nucleotides_choice='ACTG'

        

        dinucleotides_choice = {'A':'','C':'','G':'','T':''}
        for n in dinucleotides_choice:
            for dn in transition_matrix[n]:
                dinucleotides_choice[n]+= dn*int(round(transition_matrix[n][dn]))

              
       

        #Generate randome sequence
        random_seq = random.choice(nucleotides_choice)
        i=1
        while i<seq_length:
            choice = dinucleotides_choice[random_seq[i-1]]
            if choice=='':
                choice = 'ATCG' # occurs in extreme case where a single base occurs only once and at end of sequence
            random_seq+= random.choice(choice)
            i+=1

        random_sequences.append(random_seq)

    w = open(outdir+"/"+project_name+"/Run_Files/Random_"+str(mp)+'_'+str(iteration+1)+".fasta",'w',0)
    for i in range(len(random_sequences)):
        w.write('>'+str(i+1)+'\n')
        w.write(random_sequences[i]+'\n')
    w.close()
            
    return random_sequences

def calculate_pValue(project_name):
    print "Calculating p-value"

def calculate_pValue6mer(project_name):
    print "Calculating p-value of 6mer at level"


def msa_subprocess(arguments):
    global seq_lengths
    iterations = arguments[0]
    mp = arguments[1]

    
    #dictionaries for iteration counting:
    count_similar_kmers = {l:{} for l in range(1,number_of_layers+1)}
    count_exact_kmers = {l:{} for l in range(1,number_of_layers+1)}

    similar_combinations = {l:{} for l in range(1,number_of_layers+1)}
    count_similar_combinations = {l:{} for l in range(1,number_of_layers+1)}
    for x in kmers_dict[1]:
        for kmer in x[1]:
            if len(kmer[2]) in similar_combinations[kmer[3]]:
                similar_combinations[kmer[3]][len(kmer[2])].append(kmer)
            else:
                similar_combinations[kmer[3]][len(kmer[2])] = [kmer]

    total_similar_combinations = {l:{} for l in range(1,number_of_layers+1)}
    for level in similar_combinations:
        for length in similar_combinations[level]:
            total = sum([len(similar_combinations[level][x]) for x in similar_combinations[level] if x>=length])
            total_similar_combinations[level][length] = total
            count_similar_combinations[level][length] = 0

    for level in range(1,number_of_layers+1):
        inner_kmers = [x[1] for x in kmers_dict[level]]
        kmers = []
        for x in inner_kmers:
            kmers.extend([a[2] for a in x])
        kmers = list(dict.fromkeys(kmers)) #removes repeat elements
        #kmers = [x[0][2] for x in kmers_dict[level]]
        #kmers = list(dict.fromkeys(kmers)) #removes repeat elements
        for kmer in kmers:
            count_exact_kmers[level][kmer]=0.0


    for i in range(iterations):
        
        print '...Iteration: '+str(i+1) +' of '+str(iterations)+' of MSA subprocess '+str(mp)      
        iter_aligned_seqs = shuffle_alignment()    
        iter_aligned_seqs = remove_gaps(iter_aligned_seqs,mp,i)
        iter_lens = [len(x) for x in iter_aligned_seqs]
        hsps = [ [] for s in iter_aligned_seqs ]
        
        if hspblast:
            blastFile = outdir+"/"+project_name+"/Run_Files/ShuffledMSA_"+str(mp)+'_'+str(i+1)+".fasta"
            headers = ['>'+str(i+1) for i in range(len(iter_aligned_seqs))]
            intron_indices = [[] for i in range(len(iter_aligned_seqs))]
            successful_blast = LB.run_blast(outdir,project_name,blastFile,str(mp))
            if successful_blast:
                headers,iter_aligned_seqs,iter_lens,intron_indices,hsps = LB.parse_blast(outdir,project_name,blastFile,headers,iter_aligned_seqs,iter_lens,intron_indices,str(mp))

        IterationGraph = SM.start_lncLOOM(iter_aligned_seqs,iter_lens,kmers_len,stopw,hsps,solver,prune,min_depth,shorttandem,maxedges,noconstraints,tolerance)
        IterationKmers = SM.get_kmers_per_layer(IterationGraph[0],number_of_layers)

        '''iteration_combinations = {l:[] for l in range(1,number_of_layers+1)}
        for x in IterationKmers[1]:
            for kmer in x[1]:
                iteration_combinations[kmer[3]].append(kmer)
        for level in iteration_combinations:
            if len(iteration_combinations[level])>=similar_combinations[level]:
                count_similar_combinations[level]+=1'''


        iteration_combinations = {l:{} for l in range(1,number_of_layers+1)}

        for x in IterationKmers[1]:
            for kmer in x[1]:
                if len(kmer[2]) in iteration_combinations[kmer[3]]:
                    iteration_combinations[kmer[3]][len(kmer[2])].append(kmer)
                else:
                    iteration_combinations[kmer[3]][len(kmer[2])] = [kmer]

        for level in count_similar_combinations:
            for length in count_similar_combinations[level]:
                total = sum([len(iteration_combinations[level][x]) for x in iteration_combinations[level] if x>=length])
                if total >= total_similar_combinations[level][length] and total!=0:
                    count_similar_combinations[level][length]+=1



        for level in kmers_dict:
            inner_kmers = [x[1] for x in IterationKmers[level]]
            kmers = []
            for x in inner_kmers:
                kmers.extend([a[2] for a in x])
            kmers = list(dict.fromkeys(kmers)) #removes repeat elements
            #kmers = [x[0][2] for x in IterationKmers[level]]
            #kmers = list(dict.fromkeys(kmers)) #removes repeat elements - calculating chance of finding it at least once. 
            for kmer in kmers:
                if kmer in count_exact_kmers[level]:
                    count_exact_kmers[level][kmer]+=1
                '''if len(kmer) in count_similar_kmers[level]: 
                    count_similar_kmers[level][len(kmer)]+=1
                else:
                    count_similar_kmers[level][len(kmer)]=1.0'''

  
    #return (count_similar_kmers,count_exact_kmers)
    return (count_similar_combinations,count_exact_kmers)



def run_stats_msa(iterations,nol,kdict,outd,pname,klen,sw,seq_len,pr,sol,mdepth,tan,medges,noconst,multiprocess,mafft,hsp,tol):

    global number_of_layers
    global kmers_dict
    global project_name
    global outdir
    global kmers_len
    global stopw
    global seq_lengths
    global prune
    global solver
    global min_depth
    global shorttandem
    global aligned_seqs
    global maxedges
    global noconstraints
    global hspblast
    global tolerance

    number_of_layers = nol
    kmers_dict = kdict
    project_name = pname
    outdir = outd
    kmers_len = klen
    stopw = sw
    seq_lengths = seq_len
    prune = pr
    solver = sol
    min_depth = mdepth
    shorttandem = tan
    maxedges = medges
    noconstraints = noconst
    hspblast = hsp
    tolerance = tol

   
    #Run stats
    print '#############################################'
    print 'RUNNING STATS: Shuffling MSA'
    print '#############################################'

    Pvalue1 = {} #similar motif of same length at same depth
    Pvalue2 = {} #chance of finding exact motif to a depth x

    Depth_of_kmers = {}#dictionary that specifies the depth of every kmer found in graph

    for level in range(1,number_of_layers+1):
        inner_kmers = [x[1] for x in kmers_dict[level]]
        kmers = []
        for x in inner_kmers:
            kmers.extend([a[2] for a in x])
        #kmers = [x[0][2] for x in kmers_dict[level]]
        kmers = list(dict.fromkeys(kmers)) #removes repeat elements
        for kmer in kmers:
            Depth_of_kmers[kmer] = level
    #if 1==1:
    try:
        
        if mafft:
            aligned_seqs = extract_msa()
            #split total iterations into  subprocesses for multicore processing
            
            pool = Pool(processes=multiprocess)
            if iterations%multiprocess==0:
                sub_iterations = [iterations/multiprocess for i in range(multiprocess)]
            else:
                sub_iterations = [iterations/multiprocess for i in range(multiprocess)]
                remainder = iterations%multiprocess
                for i in range(remainder):
                    sub_iterations[i]+=1

            paramaters = [[sub_iterations[i],i+1] for i in range(multiprocess)]
            results = pool.map(msa_subprocess,paramaters)
            count_similar_kmers,count_exact_kmers = sum_iteration_subprocess(results,number_of_layers,kmers_dict) 
            

            #Calculate first set of P-values for each kmer

            for kmer in Depth_of_kmers:
                depth = Depth_of_kmers[kmer]
                k_len = len(kmer)
                #Total occurances of a motif with k_len or greater at level=depth
                #Pvalue1[kmer] = "{:.3f}".format(float(sum([count_similar_kmers[depth][x] for x in count_similar_kmers[depth] if x>=k_len]))/iterations)
                Pvalue1[kmer] = "{:.3f}".format(float(count_similar_kmers[depth][len(kmer)])/iterations)
                #Total occurances of a motif with k_len or greater at level=depth
                Pvalue2[kmer] = "{:.3f}".format((count_exact_kmers[depth][kmer])/iterations)


        else: #if mafft fails
            for kmer in Depth_of_kmers:
                Pvalue1[kmer] = "Undefined"
                Pvalue2[kmer] = "Undefined"
    #else:
    except:
        for kmer in Depth_of_kmers:
            Pvalue1[kmer] = "Undefined"
            Pvalue2[kmer] = "Undefined"

    return Pvalue1,Pvalue2



def get_depth_of_kmers(nol,kdict):
    Depth_of_kmers = {}#dictionary that specifies the depth of every kmer found in graph
    for level in range(1,nol+1):
        #kmers = [x[0][2] for x in kdict[level]]    
        inner_kmers = [x[1] for x in kdict[level]]
        kmers = []
        for x in inner_kmers:
            kmers.extend([a[2] for a in x])
        kmers = list(dict.fromkeys(kmers)) #removes repeat elements
        for kmer in kmers:
            Depth_of_kmers[kmer] = level
    return Depth_of_kmers


def random_subprocess(arguments):

    global seq_lengths   
    iterations = arguments[0]
    mp = arguments[1]

    #dictionaries for iteration counting:
    '''count_similar_kmers = {l:{} for l in range(1,number_of_layers+1)}
    count_exact_kmers = {l:{} for l in range(1,number_of_layers+1)}

    similar_combinations = {l:[] for l in range(1,number_of_layers+1)}
    count_similar_combinations = {l:0 for l in range(1,number_of_layers+1)}
    for x in kmers_dict[1]:
        for kmer in x[1]:
            similar_combinations[kmer[3]].append(kmer)
    for level in similar_combinations:
        similar_combinations[level] = len(similar_combinations[level])'''



    #dictionaries for iteration counting:
    count_similar_kmers = {l:{} for l in range(1,number_of_layers+1)}
    count_exact_kmers = {l:{} for l in range(1,number_of_layers+1)}

    similar_combinations = {l:{} for l in range(1,number_of_layers+1)}
    count_similar_combinations = {l:{} for l in range(1,number_of_layers+1)}
    for x in kmers_dict[1]:
        for kmer in x[1]:
            if len(kmer[2]) in similar_combinations[kmer[3]]:
                similar_combinations[kmer[3]][len(kmer[2])].append(kmer)
            else:
                similar_combinations[kmer[3]][len(kmer[2])] = [kmer]

    total_similar_combinations = {l:{} for l in range(1,number_of_layers+1)}
    for level in similar_combinations:
        for length in similar_combinations[level]:
            total = sum([len(similar_combinations[level][x]) for x in similar_combinations[level] if x>=length])
            total_similar_combinations[level][length] = total
            count_similar_combinations[level][length] = 0


    for level in range(1,number_of_layers+1):
        #kmers = [x[0][2] for x in kmers_dict[level]]
        inner_kmers = [x[1] for x in kmers_dict[level]]
        kmers = []
        for x in inner_kmers:
            kmers.extend([a[2] for a in x])
        kmers = list(dict.fromkeys(kmers)) #removes repeat elements
        for kmer in kmers:
            count_exact_kmers[level][kmer]=0.0


    for i in range(iterations):
        
        print '...Iteration: '+str(i+1) +' of '+str(iterations)+' of random subprocess '+str(mp)      
        iter_aligned_seqs = build_random_sequence_set(mp,i)
        hsps = [ [] for s in iter_aligned_seqs ]
        iter_lens = [len(x) for x in iter_aligned_seqs]
        if hspblast:
            blastFile = outdir+"/"+project_name+"/Run_Files/Random_"+str(mp)+'_'+str(i+1)+".fasta"
            headers = ['>'+str(i+1) for i in range(len(iter_aligned_seqs))]
            intron_indices = [[] for i in range(len(iter_aligned_seqs))]
            successful_blast = LB.run_blast(outdir,project_name,blastFile,str(mp))
            if successful_blast:
                headers,iter_aligned_seqs,iter_lens,intron_indices,hsps = LB.parse_blast(outdir,project_name,blastFile,headers,iter_aligned_seqs,iter_lens,intron_indices,str(mp))
                


        IterationGraph = SM.start_lncLOOM(iter_aligned_seqs,iter_lens,kmers_len,stopw,hsps,solver,prune,min_depth,shorttandem,maxedges,noconstraints,tolerance)
        IterationKmers = SM.get_kmers_per_layer(IterationGraph[0],number_of_layers)

        '''iteration_combinations = {l:[] for l in range(1,number_of_layers+1)}
        for x in IterationKmers[1]:
            for kmer in x[1]:
                iteration_combinations[kmer[3]].append(kmer)
        for level in iteration_combinations:
            if len(iteration_combinations[level])>=similar_combinations[level]:
                count_similar_combinations[level]+=1'''

        iteration_combinations = {l:{} for l in range(1,number_of_layers+1)}

        for x in IterationKmers[1]:
            for kmer in x[1]:
                if len(kmer[2]) in iteration_combinations[kmer[3]]:
                    iteration_combinations[kmer[3]][len(kmer[2])].append(kmer)
                else:
                    iteration_combinations[kmer[3]][len(kmer[2])] = [kmer]

        for level in count_similar_combinations:
            for length in count_similar_combinations[level]:
                total = sum([len(iteration_combinations[level][x]) for x in iteration_combinations[level] if x>=length])
                if total >= total_similar_combinations[level][length] and total!=0:
                    count_similar_combinations[level][length]+=1


        

        for level in kmers_dict:
            inner_kmers = [x[1] for x in IterationKmers[level]]
            kmers = []
            for x in inner_kmers:
                kmers.extend([a[2] for a in x])
            kmers = list(dict.fromkeys(kmers)) #removes repeat elements

            #kmers = [x[0][2] for x in IterationKmers[level]]
            #kmers = list(dict.fromkeys(kmers)) #removes repeat elements - calculating chance of finding it at least once. 
            for kmer in kmers:
                if kmer in count_exact_kmers[level]:
                    count_exact_kmers[level][kmer]+=1
                if len(kmer) in count_similar_kmers[level]: 
                    count_similar_kmers[level][len(kmer)]+=1
                else:
                    count_similar_kmers[level][len(kmer)]=1.0

    #return (count_similar_kmers,count_exact_kmers)
    return (count_similar_combinations,count_exact_kmers)


def sum_iteration_subprocess(list_of_results_tuples,nol,kdict):


    #Define dictionaries with keys
    total_similar_kmers = {l:{} for l in range(1,nol+1)}
    total_exact_kmers = {l:{} for l in range(1,nol+1)}

    for level in range(1,nol+1):
        inner_kmers = [x[1] for x in kdict[level]]
        kmers = []
        for x in inner_kmers:
            kmers.extend([a[2] for a in x])
        #kmers = [x[0][2] for x in kdict[level]]
        for kmer in kmers:
            total_exact_kmers[level][kmer]=0.0

    for result in list_of_results_tuples:
        similar = result[0]
        exact = result[1]
        for level in exact:
            for kmer in exact[level]:
                total_exact_kmers[level][kmer]+=exact[level][kmer]

        for level in similar:
            for kmer_len in similar[level]:
                if kmer_len in total_similar_kmers[level]:
                    total_similar_kmers[level][kmer_len]+=similar[level][kmer_len]
                else:
                    total_similar_kmers[level][kmer_len]=similar[level][kmer_len]


    return total_similar_kmers,total_exact_kmers



def run_stats_random(iterations,nol,kdict,outd,pname,klen,sw,seq,seq_len,pr,sol,mdepth,tan,medges,noconst,multiprocess,hsp,tol):

    global number_of_layers
    global kmers_dict
    global project_name
    global outdir
    global kmers_len
    global stopw
    global seq_lengths
    global prune
    global solver
    global min_depth
    global shorttandem
    global sequences
    global maxedges
    global noconstraints
    global hspblast
    global tolerance
    

    number_of_layers = nol
    kmers_dict = kdict
    project_name = pname
    outdir = outd
    kmers_len = klen
    stopw = sw
    seq_lengths = seq_len
    prune = pr
    solver = sol
    min_depth = mdepth
    shorttandem = tan
    sequences = seq
    maxedges = medges
    noconstraints = noconst
    hspblast = hsp
    tolerance =tol




    #Run stats
    print '#############################################'
    print 'RUNNING STATS: Generating Random Sequence Set'
    print '#############################################'

    Pvalue1 = {} #similar motif of same length at same depth
    Pvalue2 = {} #chance of finding exact motif to a depth x

    Depth_of_kmers = {}#dictionary that specifies the depth of every kmer found in graph
    for level in range(1,number_of_layers+1):
        inner_kmers = [x[1] for x in kmers_dict[level]]
        kmers = []
        for x in inner_kmers:
            kmers.extend([a[2] for a in x])
        kmers = list(dict.fromkeys(kmers)) #removes repeat elements
        #kmers = [x[0][2] for x in kmers_dict[level]]
        #kmers = list(dict.fromkeys(kmers)) #removes repeat elements
        for kmer in kmers:
            Depth_of_kmers[kmer] = level
    #if 1==1:       
    try:
        #split total iterations into subprocesses for multicore processing
        
        pool = Pool(processes=multiprocess)
        if iterations%multiprocess==0:
            sub_iterations = [iterations/multiprocess for i in range(multiprocess)]
        else:
            sub_iterations = [iterations/multiprocess for i in range(multiprocess)]
            remainder = iterations%multiprocess
            for i in range(remainder):
                sub_iterations[i]+=1

        paramaters = [[sub_iterations[i],i+1] for i in range(multiprocess)]
        results = pool.map(random_subprocess,paramaters)

        count_similar_kmers,count_exact_kmers = sum_iteration_subprocess(results,number_of_layers,kmers_dict) 
        #Calculate first set of P-values for each kmer

        for kmer in Depth_of_kmers:
            depth = Depth_of_kmers[kmer]
            k_len = len(kmer)
            #Total occurances of a motif with k_len or greater at level=depth
            #Pvalue1[kmer] = "{:.3f}".format(float(sum([count_similar_kmers[depth][x] for x in count_similar_kmers[depth] if x>=k_len]))/iterations)
            #Pvalue1[kmer] = "{:.3f}".format(float(count_similar_kmers[depth])/iterations)
            Pvalue1[kmer] = "{:.3f}".format(float(count_similar_kmers[depth][len(kmer)])/iterations)
            #Total occurances of a motif with k_len or greater at level=depth
            Pvalue2[kmer] = "{:.3f}".format((count_exact_kmers[depth][kmer])/iterations)
    #else:
    except: #if error occurs
        for kmer in Depth_of_kmers:
            Pvalue1[kmer] = "Undefined"
            Pvalue2[kmer] = "Undefined"

    return Pvalue1,Pvalue2


 
