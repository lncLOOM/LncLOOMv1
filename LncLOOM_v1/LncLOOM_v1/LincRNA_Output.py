import os, sys
import random
import LincMotif as SM
import networkx as nx
from operator import itemgetter

def write_homepage(outdir,project_name,level,top):

    logo = os.path.join(os.path.dirname(__file__), 'src', 'logo.png')
    html_output ='<!DOCTYPE html><html><head><style>header{background-color:white; outline: 7px solid purple;text-align:center;line-height:1}th,td {padding: 20px; border: 6px solid white;}th{font-size:35px;} td{font-size:20px; text-align:center}table{table-layout : fixed;} </style>'
    html_output+='\n<title>LncLOOM Results</title></head><body>'
    html_output+='\n<header><h1 style="font-size:45px;">'+project_name+' Results</h1><h2><img src="'+logo+'" style="vertical-align: middle;"><br>Ulitsky Lab <br>Weizmann Institute of Science</h2></header>'
    html_output+='\n<br><table style="width:100%">'
    html_output+='\n<tr><th style = "background-color:purple;color:white">All K-mers</th><th style = "background-color:DodgerBlue " >Conservation by Species</th><th style = "background-color:orangered">Selected Species</th></tr>'
    html_output+='\n<tr><td><a href = ./Html_Files/kmers_in_seqs.html style="text-decoration:none;color:purple;font-size:25px;font-weight:bold">&#9654; KMERS IN SEQUENCES</a></td>'
    html_output+='\n<td><a href = ./Html_Files/kmers_in_seqs_layer_conservation.html style="text-decoration:none;color:DodgerBlue ;font-size:25px;font-weight:bold">&#9654; KMERS MAPPED TO ANCHOR SEQ.</a></td>'
    html_output+='\n<td><a href = ./Html_Files/kmers_in_seqs_level.html style="text-decoration:none;color:orangered;font-size:25px;font-weight:bold">&#9654; '+level[0:15]+'. LEVEL SPECIFIC KMERS</a></td></tr>'
    html_output+='\n<tr><td><a href = ./Html_Files/kmers_in_blocks.html style="text-decoration:none;color:purple;font-size:25px;font-weight:bold">&#9654; BLOCK DIAGRAMS</a></td>'
    html_output+='\n<td><a href = ./Html_Files/kmers_in_blocks_layer_conservation.html style="text-decoration:none;color:DodgerBlue ;font-size:25px;font-weight:bold">&#9654; BLOCK DIAGRAMS MAPPED TO ANCHOR SEQ.</a></td>'
    html_output+='\n<td><a href = ./Html_Files/kmers_in_blocks_level.html style="text-decoration:none;color:orangered;font-size:25px;font-weight:bold">&#9654; '+level[0:15]+'. KMERS IN BLOCKS</a></td></tr>'
    html_output+='\n<tr><td></td>'
    html_output+='\n<td><a href = ./Html_Files/Modules.html style="text-decoration:none;color:DodgerBlue;font-size:25px;font-weight:bold">&#9654; MODULES</a></td></tr>'
    #html_output+='\n<tr><td><a href = ./Html_Files/kmers_in_blocks_overlap.html style="text-decoration:none;color:purple;font-size:25px;font-weight:bold">&#9654; OVERLAPPING KMERS IN BLOCKS</a></td></tr>'
    html_output+='</table></body></html>' 
    w= open(outdir+'/'+project_name+'/'+project_name+'_RESULTS.html','w',0)
    w.write(html_output)
    w.close()



def write_error_page(outdir,project_name,msg,fname):
    html_output ='<!DOCTYPE html><html><head><style>header{background-color:SkyBlue;text-align:center;}table, th, td {padding: 15px;  border: 1px solid black;border-collapse: collapse;}th {text-align: left;font-size:20px}</style><title>SPECIFIC DEPTH KMERS</title></head><body>'
    html_output+='<header><h1 style="font-size:30px">Selected Motif Conservation</h1></header>'
    html_output += '<span style = "color:red">'+msg+'</span><br></body></html>'
    w= open(outdir+"/"+project_name+"/Html_Files/"+fname,'w',0)
    w.write(html_output)
    w.close()



def write_annotation_form(outdir,project_name):
    print "Annotation Form"

def write_eCLIP(eCLIP_dict,nodes_in_chromosome,outdir,project_name,source,header,strand,chrm,kmer_colour_dict):

    #Make a sub-directory to store html output file
    html_files_path = outdir+"/"+project_name+"/Html_Files"
    if not os.path.exists(html_files_path):
        os.mkdir(html_files_path)

    

    
    html_output ='<!DOCTYPE html><html><head><style>header{background-color:SkyBlue;text-align:center;}table, th, td {padding: 15px;  border: 1px solid black;border-collapse: collapse;}th {text-align: left;font-size:20px}</style><title>eCLIP Matches ('+source+')</title></head><body>'
    html_output+='<header><h1 style="font-size:30px">eCLIP Annotation Results<br><br>('+source+')</h1><br></header><br>'

    header = header.strip('>')
    nodes = eCLIP_dict.keys()
    nodes = sorted(nodes,key=itemgetter(0))
    
    count_matches = 0
    for node in nodes:
        eCLIP = eCLIP_dict[node]
        if eCLIP:
            count_matches+=1
            colour = kmer_colour_dict[node[2]][0]
            fc = kmer_colour_dict[node[2]][1]
            node_seg_chrm = nodes_in_chromosome[node]
            nSeg = len(node_seg_chrm)
            html_output+='<table id="'+str(node)+'" style="width:100%"><tr><th  colspan="8">Match '+str(count_matches)+' in '+header+'</th></tr>'
            
            html_output+='<tr><th>Motif</th><th>Start in Seq (1 Indexed)</th><th>End in Seq (1 Indexed)</th><th>Strand</th><th>Chrm</th><th>Exon</th><th>Start in Chrm (0 Indexed)</th><th>End in Chrm (1 Indexed)</th></tr>'
            html_output+='<tr><td rowspan="'+str(nSeg)+'"><mark style="background-color:'+colour+';color:'+fc+'">'+node[2]+'</mark></td><td rowspan="'+str(nSeg)+'">'+str(node[0])+'</td><td rowspan="'+str(nSeg)+'">'+str(node[1])+'</td><td rowspan="'+str(nSeg)+'">'+strand+'</td><td rowspan="'+str(nSeg)+'">'+chrm+'</td>'
            seg = node_seg_chrm[0]
            html_output+='<td>'+str(seg[0])+'</td><td>'+str(seg[1])+'</td><td>'+str(seg[2])+'</td></tr>'
            for seg in node_seg_chrm[1:]:
                html_output+='<tr><td>'+str(seg[0])+'</td><td>'+str(seg[1])+'</td><td>'+str(seg[2])+'</td></tr>'
            html_output+='<tr><th>eCLIP Fold-Enrichment</th><th>Binding Protein</th><th>Cell Line</th><th>Strand</th><th colspan="2">Chrm</th><th>Start in Chrm (0 Indexed)</th><th>End in Chrm (1 Indexed)</th></tr>'
            for ec in eCLIP:
                e = ec[0]
                overlap = ec[1]
                if overlap==1:
                    row_style = 'style="color:red;font-weight: bold"'
                else:
                    row_style = 'style="color:black"'
                html_output+='<tr '+row_style+'><td>'+str(e[6])+'</td><td>'+e[3]+' (bg='+str(e[8])+'%)</td><td>'+e[4]+'</td><td>'+e[5]+'</td><td colspan="2">'+e[0]+'</td><td>'+str(e[1])+'</td><td>'+str(e[2])+'</td></tr>'


            html_output+='</table><br><br>'

                
    html_output+= '</body></html>'

    w= open(outdir+"/"+project_name+"/Html_Files/eCLIP_results_"+source+".html",'w',0)
    w.write(html_output)
    w.close()


def generate_colour_set(number_of_nodes):
    colour_set = []
    count = 0
    while count<number_of_nodes:
   
       r = int(random.random()*256)
       g = int(random.random()*256)
       b = int(random.random()*256)
       c = 'rgb('+str(r)+','+str(g)+','+str(b)+')'          
       if not c in colour_set:
           #Determine white or black font colour
           fc = "black"
           if (r+g+b)/3<127:
               fc = "white"
           c_set=(c,fc)
           colour_set.append(c_set)
           count+=1
    return colour_set

##################
def colour_kmers_in_seq_html_2(sequence,overlapping_nodes,colours_dict,overlap_flag,block_flag,Annotations,introns,project_name,stats_dict,Depth_of_kmers):
    #If overlap_flag is true, overlapping kmers will be individually coloured and indicated else only the combined kmer will be coloured

    kmers_in_seq_html = ' '

    if overlapping_nodes:
        node_tuples = [node[0] for node in overlapping_nodes]
        full_kmer_starts = [node[0][0] for node in overlapping_nodes]
        full_kmer_ends = [node[0][1] for node in overlapping_nodes]
        full_kmer_seqs = [node[0][2] for node in overlapping_nodes]

        spliced_sequence = []
        seg_start = 0
        for i in range(len(full_kmer_seqs)):
            kmer_start = full_kmer_starts[i]-1 # -1 for index of start
            seg = sequence[seg_start:kmer_start]
            spliced_sequence.append(seg)
            spliced_sequence.append(str(i))
            seg_start = full_kmer_ends[i]
        spliced_sequence.append(sequence[seg_start:])


        #Mark colours of each kmer segment in spliced_sequence
     
        count_bases = 0
        count_with_introns = 0

        #set submenu if no annotation was performed by user (ADD THIS FUNCTION LATER -LINK TO PAGE TO RUN ANNOTATION ANALYSIS)
        #if not Annotations:
            #button = '<a style="color:DarkRed;line-height:1.0;padding: 12px 16px;font-size:18px" href="#">&#9654; RUN ANNOTATION ANALYSIS</a>'
            
        
        for segment in spliced_sequence:
            if segment.isdigit():
                full_kmer = full_kmer_seqs[int(segment)]
                if overlap_flag:
                    inner_kmers = overlapping_nodes[int(segment)][1]
                else:
                    inner_kmers = [overlapping_nodes[int(segment)][0]]

                full_start = full_kmer_starts[int(segment)]
                #create dictionary that specifies colour of each character in the kmer
                kmer_colours = {}
                drop_down = {}
                sub_list = []
                significance = [True for i in inner_kmers]
                for b,inner in enumerate(inner_kmers):
                    relative_start = inner[0]-full_start
                    inner_colours = colours_dict[inner[2]]
                    for index in range(len(inner[2])):
                        relative_index = index+relative_start
                        kmer_colours[relative_index]=inner_colours
                        drop_down[relative_index]=b
                        
                    nodeID = (inner[0],inner[1],inner[2])


                    submenu='<br><span style="font-size:14px;line-height:0.5;padding: 12px 16px;color:black">Depth:'+str(Depth_of_kmers[nodeID[2]][0])+' ('+Depth_of_kmers[nodeID[2]][1]+')</span><br>'
                    significant = True
                    if stats_dict:
                        if float(stats_dict[nodeID[2]][1])>0.05 or float(stats_dict[nodeID[2]][3])>0.05 or float(stats_dict[nodeID[2]][2])>0.05 or float(stats_dict[nodeID[2]][4])>0.05:
                            significant = False
                        submenu+='<span style="font-size:14px;line-height:0.5;padding: 12px 16px;color:black">E<sub>i</sub>-value:'+stats_dict[nodeID[2]][1]+', P<sub>i</sub>-value:'+stats_dict[nodeID[2]][2]+'</span>'
                        submenu+='<br><span style="font-size:14px;line-height:0.5;padding: 12px 16px;color:black">E<sub>r</sub>-value:'+stats_dict[nodeID[2]][3]+', P<sub>r</sub>-value:'+stats_dict[nodeID[2]][4]+'</span><br>'
                    significance[b]=significant
      
                    if 'eCLIP' in Annotations: #Expand this
                        matches_to_eCLIP = False
                        try:
                            if 'BLAT' in Annotations['eCLIP']:
                                if Annotations['eCLIP']['BLAT']:
                                    kmer_annotations = Annotations['eCLIP']['BLAT'][nodeID]
                                    if kmer_annotations:
                                        included=[]
                                        submenu+='<span style="font-size:18px;line-height:1.0;padding: 12px 16px;color:DarkBlue">eCLIP MATCHES</span>'
                                        matches_to_eCLIP = True
                                        for annotation in kmer_annotations:
                                            anno = annotation[0]
                                            if anno[3] not in included:
                                                included.append(anno[3])
                                                submenu+= '<a style="color:DarkBlue" href="eCLIP_results_BLAT.html#'+str(nodeID)+'">&#9654;'+anno[3]+' (bg='+str(anno[8])+'%)</a>'
                     
                            if 'BED' in Annotations['eCLIP']:
                                if Annotations['eCLIP']['BED']:
                                    kmer_annotations = Annotations['eCLIP']['BED'][nodeID]
                                    if kmer_annotations:
                                        included=[]
                                        submenu+='<span style="font-size:18px;line-height:1.0;padding: 12px 16px;color:DarkBlue">eCLIP MATCHES</span>'
                                        matches_to_eCLIP = True
                                        for annotation in kmer_annotations:
                                            anno = annotation[0]
                                            if anno[3] not in included:
                                                included.append(anno[3])                                        
                                                submenu+= '<a style="color:DarkBlue" href="eCLIP_results_BED.html#'+str(nodeID)+'">&#9654;'+anno[3]+' (bg='+str(anno[8])+'%)</a>'
                        except Exception as e:
                            print str(e)

                        if not matches_to_eCLIP:
                            submenu+='<a style="color:DarkBlue">No matches to eCLIP Data</a>'
                    
                    if 'TargetScan' in Annotations:
                        kmer_annotations = Annotations['TargetScan'][nodeID]
                        if kmer_annotations:
                            submenu+='<span style="font-size:18px;line-height:1.0;padding: 12px 16px;color:DarkRed">MATCHES To TargetScan</span>'
                            for anno in kmer_annotations:
                                seedID = anno.split(':')[1].strip()
                                submenu+='<a style="color:DarkRed" href="miRNA_Matches.html#'+seedID+'">&#9654; '+anno+'</a>'
                        else:
                            submenu+='<a style="color:DarkRed">No matches to TargetScan</a>'
                    sub_list.append(submenu)


                seg_output = '<div class="dropdown"><span class="dropbtn">'
                marked_kmer = ''

                for index,base in enumerate(full_kmer):
                    kcolor = kmer_colours[index][0]
                    if block_flag:
                        fc = kmer_colours[index][0]
                    else:
                        fc = kmer_colours[index][1]
                    si = drop_down[index]
                    submenu = sub_list[si]
                    i_kmer = inner_kmers[si][2]
                    significant = significance[si]
                    if significant:
                        i_kmer = i_kmer.upper()
                        base = base.upper() 
                    else:
                        i_kmer = i_kmer.lower()
                        base = base.lower() 
                    marked_kmer+='<mark style="background-color:'+kcolor+';color:'+fc+'">'+base+'</mark>'

                    count_bases+=1
                    count_with_introns+=1

                    if count_with_introns%120==0:
                        marked_kmer+='</span><div class="dropdown-content"><span style="font-size:18px;line-height:1.0;padding: 12px 16px;color:black">'+i_kmer+'</span>'+submenu+'</div></div><span style ="font-size:15px"> '+str(count_bases)+'</span><br> <div class="dropdown"><span class="dropbtn">'
                    if count_bases in introns:
                        for i in range(1,3):
                            marked_kmer+= '<span style = "color:red"><b>|</b></span>'
                            count_with_introns+=1
                            if count_with_introns%120==0:
                                marked_kmer+='</span><div class="dropdown-content"><span style="font-size:18px;line-height:1.0;padding: 12px 16px;color:black">'+i_kmer+'</span>'+submenu+'</div></div><span style ="font-size:15px"> '+str(count_bases)+'</span><br> <div class="dropdown"><span class="dropbtn">'
                   
                    if index<len(full_kmer)-1:
                        if si!=drop_down[index+1]:
                            seg_output+=marked_kmer+'</span><div class="dropdown-content"><span style="font-size:18px;line-height:1.0;padding: 12px 16px;color:black">'+i_kmer+'</span>'+submenu+'</div></div><div class="dropdown"><span class="dropbtn">'
                            marked_kmer = ''
                seg_output+=marked_kmer+'</span><div class="dropdown-content"><span style="font-size:18px;line-height:1.0;padding: 12px 16px;color:black">'+i_kmer+'</span>'+submenu+'</div></div>'
                kmers_in_seq_html+=seg_output
            else:
                seg_with_break = ''
                for base in segment:
                    seg_with_break+=base
                    count_bases+=1
                    count_with_introns+=1                            
                    if count_with_introns%120==0:
                        seg_with_break+='<span style ="font-size:15px"> '+str(count_bases)+'</span><br> '

                    if count_bases in introns:
                        for i in range(1,3):
                            seg_with_break+='<span style = "color:red"><b>|</b></span>'
                            count_with_introns+=1
                            if count_with_introns%120==0:
                                seg_with_break+='<span style ="font-size:15px"> '+str(count_bases)+'</span><br> '
                                 
                kmers_in_seq_html+=seg_with_break

      

    else:
        kmers_in_seq_html+= '<span style = "color:red">NO CONSERVED NODES FOUND</span><br>'
        count_bases=0
        count_with_introns = 0
        for s in sequence:
            kmers_in_seq_html+=s
            count_bases+=1
            count_with_introns+=1           
            if count_with_introns%120==0:
                kmers_in_seq_html+='<span style ="font-size:15px"> '+str(count_bases)+'</span><br> '
            if count_bases in introns:
                for i in range(1,3):
                    kmers_in_seq_html+='<span style = "color:red"><b>|</b></span>'
                    count_with_introns+=1
                    if count_with_introns%120==0:
                        kmers_in_seq_html+='<span style ="font-size:15px"> '+str(count_bases)+'</span><br> '          
        
    kmers_in_seq_html+=' '*(120-(count_with_introns%120))+' <span style ="font-size:15px">'+str(count_bases)+'</span>'      
    return kmers_in_seq_html













###################








def colour_kmers_in_seq_html(sequence,overlapping_nodes,colours_dict,overlap_flag,block_flag,Annotations,introns,project_name,stats_dict,Depth_of_kmers):
    #If overlap_flag is true, overlapping kmers will be individually coloured and indicated else only the combined kmer will be coloured

    kmers_in_seq_html = ' '

    if overlapping_nodes:
        node_tuples = [node[0] for node in overlapping_nodes]
        full_kmer_starts = [node[0][0] for node in overlapping_nodes]
        full_kmer_ends = [node[0][1] for node in overlapping_nodes]
        full_kmer_seqs = [node[0][2] for node in overlapping_nodes]

        spliced_sequence = []
        seg_start = 0
        for i in range(len(full_kmer_seqs)):
            kmer_start = full_kmer_starts[i]-1 # -1 for index of start
            seg = sequence[seg_start:kmer_start]
            spliced_sequence.append(seg)
            spliced_sequence.append(str(i))
            seg_start = full_kmer_ends[i]
        spliced_sequence.append(sequence[seg_start:])


        #Mark colours of each kmer segment in spliced_sequence
     
        count_bases = 0
        count_with_introns = 0

        #set submenu if no annotation was performed by user (ADD THIS FUNCTION LATER -LINK TO PAGE TO RUN ANNOTATION ANALYSIS)
        #if not Annotations:
            #button = '<a style="color:DarkRed;line-height:1.0;padding: 12px 16px;font-size:18px" href="#">&#9654; RUN ANNOTATION ANALYSIS</a>'
            
        
        for segment in spliced_sequence:
            if segment.isdigit():
                full_kmer = full_kmer_seqs[int(segment)]
                nodeID = node_tuples[int(segment)]
                if overlap_flag:
                    inner_kmers = overlapping_nodes[int(segment)][1]
                else:
                    inner_kmers = [overlapping_nodes[int(segment)][0]]

                full_start = full_kmer_starts[int(segment)]
                #create dictionary that specifies colour of each character in the kmer
                kmer_colours = {}
                for inner in inner_kmers:
                    relative_start = inner[0]-full_start
                    inner_colours = colours_dict[inner[2]]
                    for index in range(len(inner[2])):
                        relative_index = index+relative_start
                        kmer_colours[relative_index]=inner_colours

                #print overlapping_nodes
                #Handle annotations for the current kmer
                submenu='<br><span style="font-size:14px;line-height:0.5;padding: 12px 16px;color:black">Depth:'+str(Depth_of_kmers[nodeID[2]][0])+' ('+Depth_of_kmers[nodeID[2]][1]+')</span><br>'
                significant = True
                if stats_dict:
                    if float(stats_dict[nodeID[2]][1])>0.05 or float(stats_dict[nodeID[2]][3])>0.05 or float(stats_dict[nodeID[2]][2])>0.05 or float(stats_dict[nodeID[2]][4])>0.05:
                        significant = False
                    submenu+='<span style="font-size:14px;line-height:0.5;padding: 12px 16px;color:black">E<sub>i</sub>-value:'+stats_dict[nodeID[2]][1]+', P<sub>i</sub>-value:'+stats_dict[nodeID[2]][2]+'</span>'
                    submenu+='<br><span style="font-size:14px;line-height:0.5;padding: 12px 16px;color:black">E<sub>r</sub>-value:'+stats_dict[nodeID[2]][3]+', P<sub>r</sub>-value:'+stats_dict[nodeID[2]][4]+'</span><br>'

      
                if 'eCLIP' in Annotations: #Expand this
                    matches_to_eCLIP = False
                    try:
                        if 'BLAT' in Annotations['eCLIP']:
                            if Annotations['eCLIP']['BLAT']:
                                kmer_annotations = Annotations['eCLIP']['BLAT'][nodeID]
                                if kmer_annotations:
                                    included=[]
                                    submenu+='<span style="font-size:18px;line-height:1.0;padding: 12px 16px;color:DarkBlue">eCLIP MATCHES</span>'
                                    matches_to_eCLIP = True
                                    for annotation in kmer_annotations:
                                        anno = annotation[0]
                                        if anno[3] not in included:
                                            included.append(anno[3])
                                            submenu+= '<a style="color:DarkBlue" href="eCLIP_results_BLAT.html#'+str(nodeID)+'">&#9654;'+anno[3]+' (bg='+str(anno[8])+'%)</a>'

                        
                        if 'BED' in Annotations['eCLIP']:
                            if Annotations['eCLIP']['BED']:
                                kmer_annotations = Annotations['eCLIP']['BED'][nodeID]
                                if kmer_annotations:
                                    included=[]
                                    submenu+='<span style="font-size:18px;line-height:1.0;padding: 12px 16px;color:DarkBlue">eCLIP MATCHES</span>'
                                    matches_to_eCLIP = True
                                    for annotation in kmer_annotations:
                                        anno = annotation[0]
                                        if anno[3] not in included:
                                           included.append(anno[3])                                        
                                           submenu+= '<a style="color:DarkBlue" href="eCLIP_results_BED.html#'+str(nodeID)+'">&#9654;'+anno[3]+' (bg='+str(anno[8])+'%)</a>'
                    except Exception as e:
                        print str(e)

                    if not matches_to_eCLIP:
                        submenu+='<a style="color:DarkBlue">No matches to eCLIP Data</a>'
                    
                    
                
                if 'TargetScan' in Annotations:
                    kmer_annotations = Annotations['TargetScan'][nodeID]
                    if kmer_annotations:
                        submenu+='<span style="font-size:18px;line-height:1.0;padding: 12px 16px;color:DarkRed">MATCHES To TargetScan</span>'
                        for anno in kmer_annotations:
                            seedID = anno.split(':')[1].strip()
                            submenu+='<a style="color:DarkRed" href="miRNA_Matches.html#'+seedID+'">&#9654; '+anno+'</a>'
                    else:
                        submenu+='<a style="color:DarkRed">No matches to TargetScan</a>'

                marked_kmer = ''

                if not significant:
                    full_kmer = full_kmer.lower() 
                for index,base in enumerate(full_kmer):
                    kcolor = kmer_colours[index][0]
                    if block_flag:
                        fc = kmer_colours[index][0]
                    else:
                        fc = kmer_colours[index][1]
                    marked_kmer+='<mark style="background-color:'+kcolor+';color:'+fc+'">'+base+'</mark>'
                    count_bases+=1
                    count_with_introns+=1
                    if count_with_introns%120==0:
                        marked_kmer+='</span><div class="dropdown-content"><span style="font-size:18px;line-height:1.0;padding: 12px 16px;color:black">'+full_kmer+'</span>'+submenu+'</div></div><span style ="font-size:15px"> '+str(count_bases)+'</span><br> <div class="dropdown"><span class="dropbtn">'
                    if count_bases in introns:
                        for i in range(1,3):
                            marked_kmer+= '<span style = "color:red"><b>|</b></span>'
                            count_with_introns+=1
                            if count_with_introns%120==0:
                                marked_kmer+='</span><div class="dropdown-content"><span style="font-size:18px;line-height:1.0;padding: 12px 16px;color:black">'+full_kmer+'</span>'+submenu+'</div></div><span style ="font-size:15px"> '+str(count_bases)+'</span><br> <div class="dropdown"><span class="dropbtn">'
                   


                kmers_in_seq_html+='<div class="dropdown"><span class="dropbtn">'+marked_kmer+'</span><div class="dropdown-content"><span style="font-size:18px;line-height:1.0;padding: 12px 16px;color:black">'+full_kmer+'</span>'+submenu+'</div></div>'
                
            else:
                seg_with_break = ''
                for base in segment:
                    seg_with_break+=base
                    count_bases+=1
                    count_with_introns+=1                            
                    if count_with_introns%120==0:
                        seg_with_break+='<span style ="font-size:15px"> '+str(count_bases)+'</span><br> '

                    if count_bases in introns:
                        for i in range(1,3):
                            seg_with_break+='<span style = "color:red"><b>|</b></span>'
                            count_with_introns+=1
                            if count_with_introns%120==0:
                                seg_with_break+='<span style ="font-size:15px"> '+str(count_bases)+'</span><br> '
                                 
                kmers_in_seq_html+=seg_with_break

      

    else:
        kmers_in_seq_html+= '<span style = "color:red">NO CONSERVED NODES FOUND</span><br>'
        count_bases=0
        count_with_introns = 0
        for s in sequence:
            kmers_in_seq_html+=s
            count_bases+=1
            count_with_introns+=1           
            if count_with_introns%120==0:
                kmers_in_seq_html+='<span style ="font-size:15px"> '+str(count_bases)+'</span><br> '
            if count_bases in introns:
                for i in range(1,3):
                    kmers_in_seq_html+='<span style = "color:red"><b>|</b></span>'
                    count_with_introns+=1
                    if count_with_introns%120==0:
                        kmers_in_seq_html+='<span style ="font-size:15px"> '+str(count_bases)+'</span><br> '          
        
    kmers_in_seq_html+=' '*(120-(count_with_introns%120))+' <span style ="font-size:15px">'+str(count_bases)+'</span>'      
    return kmers_in_seq_html



def write_modules(Modules,headers,colours_dict,AnnotationsDict,outdir,project_name,graded_dict,graded_colours,stats_dict,Depth_of_kmers):

    #establish a dictionary of combined kmers and their segments (for colouring purpose)
    '''segments_dict = {}
    for layer in kmers_dict:
        nodes_of_layer = kmers_dict[layer]
        layer_seg = {}
        for node in nodes_of_layer:
            combined_node = node[0]
            segments = node[1]
            layer_seg[combined_node]= segments
        segments_dict[layer] = layer_seg'''

    html_key='<table style ="width:70%;margin-left:15%;margin-right:15%;border-collapse: collapse;table-layout: fixed"><tr>'
    levels = graded_colours.keys()
    #levels.sort(reverse=True)
    levels.sort()
    for level in levels:
        html_key+='<td style="background-color:'+graded_colours[level][0]+';color:'+graded_colours[level][0]+';text-align:center">|</td>'
    html_key+='</tr>'
    html_key+='<tr>'
    for level in levels:
        if level==levels[0]:
            html_key+='<td style="color:black;text-align:left">'+str(level)+'</td>'
        elif level==levels[-1]:
            html_key+='<td style="color:black;text-align:right">'+str(level)+'</td>'
        elif level==levels[len(levels)/2]:
            html_key+='<td style="color:black;text-align:center">'+str(level)+'</td>'
        else:
            html_key+='<td></td>'
    html_key+='</tr><tr><td style="color:black;text-align:center" colspan="'+str(len(levels))+'">Depth of motif conservation (number of species)</td></tr></table>'

    #<a href="#'+str(level)+'" style="cursor:pointer;text-decoration:none;color:black;font-size: 20px">&#9660</a>

    mod_levels = Modules.keys()
    mod_levels.sort(reverse=True)

    submenu = ''
    for depth in mod_levels:
        submenu+='<a href="#'+str(depth)+'" style="color:black;text-align:left">&#9654;'+headers[depth-1][0:25].strip('>')+' (depth:'+str(depth)+')</a>'

    #Write modules.html
    output ='<!DOCTYPE html><html><head><style>header{background-color:SkyBlue;text-align:center;} .dropbtn {font-size: 22px;border: none;cursor: pointer;}.dropdown {position: relative;display: inline-block;}.dropdown-content {display: none;position: absolute;background-color: rgb(255, 153, 153);white-space: nowrap; overflow:hidden;box-shadow: 0px 8px 50px 0px rgba(0,0,0,0.2);z-index: 1;}.dropdown-content a {line-height:1.0;font-size: 15px;padding: 12px 16px;text-decoration: none;display: block;}.dropdown:hover .dropdown-content {display: block;}.dropdown-content a:hover {background-color: yellow;} </style><title>MODULES</title></head><body>\n'
    output+='<header><h1 style="font-size:30px">MODULES</h1><h2>start-<span style="color:red">MOTIF</span><span style="color:black">-end</span><span style="color:grey">--(number of bases in between)--</span><span style="color:black">start-</span><span style="color:red">MOTIF</span>-end</h2>'
    output+='<br><div class="dropdown"><span class="dropbtn"><mark style="background-color:LightGrey;color:black;font-size:20px;padding: 15px 15px 15px 15px;text-align: center;border: 3px solid black;">NAVIGATE &#9660</mark></span><div class="dropdown-content">'+submenu+'</div></div><br><br><br></header>\n'
    output+='<br><br><br>'+html_key+'<br><br>'


    
    


    for depth in mod_levels:
        modules_at_depth = Modules[depth]
        main_modules = Modules[depth]['Main']
        prime5 = Modules[depth]['5prime']
        prime3 = Modules[depth]['3prime']

        output+= '<h1 id="'+str(depth)+'">Modules conserved to  '+headers[depth-1].strip('>')+' (Depth: '+str(depth)+')</h1>'
        if main_modules:
            output+= '<h2>Modules in Main Graph (All sequences considered):</h2>'
        mod_keys = main_modules.keys()
        mod_keys.sort()
        for mod_id in mod_keys: #Output longest module first, with shorter fragments of module following
            mod = main_modules[mod_id]
            #Set the fontsize, linespacing ect for html display                   
            output+='<pre style="font-size:20px;line-height:1.5"><br>'
            reg_exp = mod['RE']            

            #Display all modules clustered into the current group 
            mod_layers = mod.keys() # extracts a dictionary of the loctaion of the current module in respective sequences(layers)
            mod_layers.remove('RE')

            for layer in mod_layers:
                head = headers[layer-1]
                output+=head+'<br>    '
                seq_nodes_of_mod = mod[layer][0] # extracts a list of tuples for combined kmers eg[(5,12,'ATTTCCCC'),(20,28,'ATTTCTCTC')] that constitute the particular module in the current layer
                inner_kmers_layer = mod[layer][1]
                Annotations = AnnotationsDict[layer]
                for n,node in enumerate(seq_nodes_of_mod):
                    start = node[0]
                    kmer = node[2]
                    inner_kmers = inner_kmers_layer[n]

                    #create dictionary that specifies colour of each character in the kmer
                    kmer_colours = {}
                    drop_down = {}
                    sub_list = []
                    significance = [True for i in inner_kmers]

                    for b,inner in enumerate(inner_kmers):
                        relative_start = inner[0]-start
                        inner_colours = graded_dict[inner[2]]
                        for index in range(len(inner[2])):
                            relative_index = index+relative_start
                            kmer_colours[relative_index]=inner_colours
                            drop_down[relative_index]=b
                        
                        nodeID = (inner[0],inner[1],inner[2])


                        submenu='<br><span style="font-size:14px;line-height:0.5;padding: 12px 16px;color:black">Depth:'+str(Depth_of_kmers[nodeID[2]][0])+' ('+Depth_of_kmers[nodeID[2]][1]+')</span><br>'
                        significant = True
                        if stats_dict:
                            if float(stats_dict[nodeID[2]][1])>0.05 or float(stats_dict[nodeID[2]][3])>0.05 or float(stats_dict[nodeID[2]][2])>0.05 or float(stats_dict[nodeID[2]][4])>0.05:
                                significant = False
                            submenu+='<span style="font-size:14px;line-height:0.5;padding: 12px 16px;color:black">E<sub>i</sub>-value:'+stats_dict[nodeID[2]][1]+', P<sub>i</sub>-value:'+stats_dict[nodeID[2]][2]+'</span>'
                            submenu+='<br><span style="font-size:14px;line-height:0.5;padding: 12px 16px;color:black">E<sub>r</sub>-value:'+stats_dict[nodeID[2]][3]+', P<sub>r</sub>-value:'+stats_dict[nodeID[2]][4]+'</span><br>'
                        significance[b]=significant
      
                        if 'eCLIP' in Annotations: #Expand this
                            matches_to_eCLIP = False
                            try:
                                if 'BLAT' in Annotations['eCLIP']:
                                    if Annotations['eCLIP']['BLAT']:
                                        kmer_annotations = Annotations['eCLIP']['BLAT'][nodeID]
                                        if kmer_annotations:
                                            included=[]
                                            submenu+='<span style="font-size:18px;line-height:1.0;padding: 12px 16px;color:DarkBlue">eCLIP MATCHES</span>'
                                            matches_to_eCLIP = True
                                            for annotation in kmer_annotations:
                                                anno = annotation[0]
                                                if anno[3] not in included:
                                                    included.append(anno[3])
                                                    submenu+= '<a style="color:DarkBlue" href="eCLIP_results_BLAT.html#'+str(nodeID)+'">&#9654;'+anno[3]+' (bg='+str(anno[8])+'%)</a>'
                     
                                if 'BED' in Annotations['eCLIP']:
                                    if Annotations['eCLIP']['BED']:
                                        kmer_annotations = Annotations['eCLIP']['BED'][nodeID]
                                        if kmer_annotations:
                                            included=[]
                                            submenu+='<span style="font-size:18px;line-height:1.0;padding: 12px 16px;color:DarkBlue">eCLIP MATCHES</span>'
                                            matches_to_eCLIP = True
                                            for annotation in kmer_annotations:
                                                anno = annotation[0]
                                                if anno[3] not in included:
                                                    included.append(anno[3])                                        
                                                    submenu+= '<a style="color:DarkBlue" href="eCLIP_results_BED.html#'+str(nodeID)+'">&#9654;'+anno[3]+' (bg='+str(anno[8])+'%)</a>'
                            except Exception as e:
                                print str(e)

                            if not matches_to_eCLIP:
                                submenu+='<a style="color:DarkBlue">No matches to eCLIP Data</a>'
                    
                        if 'TargetScan' in Annotations:

                            kmer_annotations = Annotations['TargetScan'][nodeID]
                            if kmer_annotations:
                                submenu+='<span style="font-size:18px;line-height:1.0;padding: 12px 16px;color:DarkRed">MATCHES To TargetScan</span>'
                                for anno in kmer_annotations:
                                    seedID = anno.split(':')[1].strip()
                                    submenu+='<a style="color:DarkRed" href="miRNA_Matches.html#'+seedID+'">&#9654; '+anno+'</a>'
                            else:
                                submenu+='<a style="color:DarkRed">No matches to TargetScan</a>'
                        sub_list.append(submenu)


                    '''for inner in inner_kmers:

                        relative_start = inner[0]-start
                        inner_colours = graded_dict[inner[2]]
                        for index in range(len(inner[2])):
                            relative_index = index+relative_start
                            kmer_colours[relative_index]=inner_colours
                    marked_kmer = ''
                    marked_kmer_graded = ''
                    for index,base in enumerate(kmer):
                        kcolor = kmer_colours[index][0]
                        fc = kmer_colours[index][1]
                        marked_kmer_graded+='<mark style="background-color:'+kcolor+';color:'+fc+'">'+base+'</mark>'

                    kcolor = colours_dict[kmer][0]
                    fc = colours_dict[kmer][1]
                    marked_kmer+='<mark style="background-color:'+kcolor+';color:'+fc+'">'+kmer+'</mark>'
                    #print inner'''


                    seg_output = '<div class="dropdown"><span class="dropbtn">'
                    marked_kmer = ''

                    for index,base in enumerate(kmer):
                        kcolor = kmer_colours[index][0]
                        fc = kmer_colours[index][1]
                        si = drop_down[index]
                        submenu = sub_list[si]
                        i_kmer = inner_kmers[si][2]
                        significant = significance[si]
                        if significant:
                            i_kmer = i_kmer.upper()
                            base = base.upper() 
                        else:
                            i_kmer = i_kmer.lower()
                            base = base.lower() 
                        marked_kmer+='<mark style="background-color:'+kcolor+';color:'+fc+'">'+base+'</mark>'



                        if index<len(kmer)-1:
                            if si!=drop_down[index+1]:
                                seg_output+=marked_kmer+'</span><div class="dropdown-content"><span style="font-size:18px;line-height:1.0;padding: 12px 16px;color:black">'+i_kmer+'</span>'+submenu+'</div></div><div class="dropdown"><span class="dropbtn">'
                                marked_kmer = ''
                    seg_output+=marked_kmer+'</span><div class="dropdown-content"><span style="font-size:18px;line-height:1.0;padding: 12px 16px;color:black">'+i_kmer+'</span>'+submenu+'</div></div>'

                    if n==0:
                        end = node[1]
                        output+=' '*(6-len(str(start)))+str(start)+'-'+seg_output
                    else:
                        
                        basesBetween = (start-end)-1
                        #output+='-'+str(end)+'<span style="color:red">--('+str(basesBetween)+')--</span>'+str(start)+'-'+marked_kmer

                        output+='-'+str(end)+'<span style="color:grey">'+'--('+str(basesBetween)+')--</span>'+str(start)+'-'+seg_output
                        end = node[1]              
                            
                output+='-'+str(end)+'<br><br>'
            

        if prime5:
            output+='</pre>\n'
            output+= "<h2>Modules in 5' Extended Graph: </h2>"
        mod_keys = prime5.keys()
        mod_keys.sort()
        for mod_id in mod_keys: #Output longest module first, with shorter fragments of module following
            mod = prime5[mod_id]
            #Set the fontsize, linespacing ect for html display                   
            output+='<pre style="font-size:20px;line-height:1.5"><br>'
            reg_exp = mod['RE']            

            #Display all modules clustered into the current group 
            mod_layers = mod.keys() # extracts a dictionary of the loctaion of the current module in respective sequences(layers)
            mod_layers.remove('RE')

            for layer in mod_layers:
                head = headers[layer-1]
                output+=head+'<br>    '
                seq_nodes_of_mod = mod[layer][0] # extracts a list of tuples for combined kmers eg[(5,12,'ATTTCCCC'),(20,28,'ATTTCTCTC')] that constitute the particular module in the current layer
                inner_kmers_layer = mod[layer][1]
                Annotations = AnnotationsDict[layer]
                for n,node in enumerate(seq_nodes_of_mod):
                    start = node[0]
                    kmer = node[2]
                    inner_kmers = inner_kmers_layer[n]

                    #create dictionary that specifies colour of each character in the kmer
                    kmer_colours = {}
                    drop_down = {}
                    sub_list = []
                    significance = [True for i in inner_kmers]

                    for b,inner in enumerate(inner_kmers):
                        relative_start = inner[0]-start
                        inner_colours = graded_dict[inner[2]]
                        for index in range(len(inner[2])):
                            relative_index = index+relative_start
                            kmer_colours[relative_index]=inner_colours
                            drop_down[relative_index]=b
                        
                        nodeID = (inner[0],inner[1],inner[2])


                        submenu='<br><span style="font-size:14px;line-height:0.5;padding: 12px 16px;color:black">Depth:'+str(Depth_of_kmers[nodeID[2]][0])+' ('+Depth_of_kmers[nodeID[2]][1]+')</span><br>'
                        significant = True
                        if stats_dict:
                            if float(stats_dict[nodeID[2]][1])>0.05 or float(stats_dict[nodeID[2]][3])>0.05 or float(stats_dict[nodeID[2]][2])>0.05 or float(stats_dict[nodeID[2]][4])>0.05:
                                significant = False
                            submenu+='<span style="font-size:14px;line-height:0.5;padding: 12px 16px;color:black">E<sub>i</sub>-value:'+stats_dict[nodeID[2]][1]+', P<sub>i</sub>-value:'+stats_dict[nodeID[2]][2]+'</span>'
                            submenu+='<br><span style="font-size:14px;line-height:0.5;padding: 12px 16px;color:black">E<sub>r</sub>-value:'+stats_dict[nodeID[2]][3]+', P<sub>r</sub>-value:'+stats_dict[nodeID[2]][4]+'</span><br>'
                        significance[b]=significant
      
                        if 'eCLIP' in Annotations: #Expand this
                            matches_to_eCLIP = False
                            try:
                                if 'BLAT' in Annotations['eCLIP']:
                                    if Annotations['eCLIP']['BLAT']:
                                        kmer_annotations = Annotations['eCLIP']['BLAT'][nodeID]
                                        if kmer_annotations:
                                            included=[]
                                            submenu+='<span style="font-size:18px;line-height:1.0;padding: 12px 16px;color:DarkBlue">eCLIP MATCHES</span>'
                                            matches_to_eCLIP = True
                                            for annotation in kmer_annotations:
                                                anno = annotation[0]
                                                if anno[3] not in included:
                                                    included.append(anno[3])
                                                    submenu+= '<a style="color:DarkBlue" href="eCLIP_results_BLAT.html#'+str(nodeID)+'">&#9654;'+anno[3]+' (bg='+str(anno[8])+'%)</a>'
                     
                                if 'BED' in Annotations['eCLIP']:
                                    if Annotations['eCLIP']['BED']:
                                        kmer_annotations = Annotations['eCLIP']['BED'][nodeID]
                                        if kmer_annotations:
                                            included=[]
                                            submenu+='<span style="font-size:18px;line-height:1.0;padding: 12px 16px;color:DarkBlue">eCLIP MATCHES</span>'
                                            matches_to_eCLIP = True
                                            for annotation in kmer_annotations:
                                                anno = annotation[0]
                                                if anno[3] not in included:
                                                    included.append(anno[3])                                        
                                                    submenu+= '<a style="color:DarkBlue" href="eCLIP_results_BED.html#'+str(nodeID)+'">&#9654;'+anno[3]+' (bg='+str(anno[8])+'%)</a>'
                            except Exception as e:
                                print str(e)

                            if not matches_to_eCLIP:
                                submenu+='<a style="color:DarkBlue">No matches to eCLIP Data</a>'
                    
                        if 'TargetScan' in Annotations:

                            kmer_annotations = Annotations['TargetScan'][nodeID]
                            if kmer_annotations:
                                submenu+='<span style="font-size:18px;line-height:1.0;padding: 12px 16px;color:DarkRed">MATCHES To TargetScan</span>'
                                for anno in kmer_annotations:
                                    seedID = anno.split(':')[1].strip()
                                    submenu+='<a style="color:DarkRed" href="miRNA_Matches.html#'+seedID+'">&#9654; '+anno+'</a>'
                            else:
                                submenu+='<a style="color:DarkRed">No matches to TargetScan</a>'
                        sub_list.append(submenu)


                    '''for inner in inner_kmers:

                        relative_start = inner[0]-start
                        inner_colours = graded_dict[inner[2]]
                        for index in range(len(inner[2])):
                            relative_index = index+relative_start
                            kmer_colours[relative_index]=inner_colours
                    marked_kmer = ''
                    marked_kmer_graded = ''
                    for index,base in enumerate(kmer):
                        kcolor = kmer_colours[index][0]
                        fc = kmer_colours[index][1]
                        marked_kmer_graded+='<mark style="background-color:'+kcolor+';color:'+fc+'">'+base+'</mark>'

                    kcolor = colours_dict[kmer][0]
                    fc = colours_dict[kmer][1]
                    marked_kmer+='<mark style="background-color:'+kcolor+';color:'+fc+'">'+kmer+'</mark>'
                    #print inner'''


                    seg_output = '<div class="dropdown"><span class="dropbtn">'
                    marked_kmer = ''

                    for index,base in enumerate(kmer):
                        kcolor = kmer_colours[index][0]
                        fc = kmer_colours[index][1]
                        si = drop_down[index]
                        submenu = sub_list[si]
                        i_kmer = inner_kmers[si][2]
                        significant = significance[si]
                        if significant:
                            i_kmer = i_kmer.upper()
                            base = base.upper() 
                        else:
                            i_kmer = i_kmer.lower()
                            base = base.lower() 
                        marked_kmer+='<mark style="background-color:'+kcolor+';color:'+fc+'">'+base+'</mark>'



                        if index<len(kmer)-1:
                            if si!=drop_down[index+1]:
                                seg_output+=marked_kmer+'</span><div class="dropdown-content"><span style="font-size:18px;line-height:1.0;padding: 12px 16px;color:black">'+i_kmer+'</span>'+submenu+'</div></div><div class="dropdown"><span class="dropbtn">'
                                marked_kmer = ''
                    seg_output+=marked_kmer+'</span><div class="dropdown-content"><span style="font-size:18px;line-height:1.0;padding: 12px 16px;color:black">'+i_kmer+'</span>'+submenu+'</div></div>'

                    if n==0:
                        end = node[1]
                        output+=' '*(6-len(str(start)))+str(start)+'-'+seg_output
                    else:
                        
                        basesBetween = (start-end)-1
                        #output+='-'+str(end)+'<span style="color:red">--('+str(basesBetween)+')--</span>'+str(start)+'-'+marked_kmer

                        output+='-'+str(end)+'<span style="color:grey">'+'--('+str(basesBetween)+')--</span>'+str(start)+'-'+seg_output
                        end = node[1]              
                            
                output+='-'+str(end)+'<br><br>'
            

        if prime3:
            output+='</pre>\n'
            output+= "<h2>Modules in 3' Extended Graph: </h2>"
        mod_keys = prime3.keys()
        mod_keys.sort()
        for mod_id in mod_keys: #Output longest module first, with shorter fragments of module following
            mod = prime3[mod_id]
            #Set the fontsize, linespacing ect for html display                   
            output+='<pre style="font-size:20px;line-height:1.5"><br>'
            reg_exp = mod['RE']            

            #Display all modules clustered into the current group 
            mod_layers = mod.keys() # extracts a dictionary of the loctaion of the current module in respective sequences(layers)
            mod_layers.remove('RE')

            for layer in mod_layers:
                head = headers[layer-1]
                output+=head+'<br>    '
                seq_nodes_of_mod = mod[layer][0] # extracts a list of tuples for combined kmers eg[(5,12,'ATTTCCCC'),(20,28,'ATTTCTCTC')] that constitute the particular module in the current layer
                inner_kmers_layer = mod[layer][1]
                Annotations = AnnotationsDict[layer]
                for n,node in enumerate(seq_nodes_of_mod):
                    start = node[0]
                    kmer = node[2]
                    inner_kmers = inner_kmers_layer[n]

                    #create dictionary that specifies colour of each character in the kmer
                    kmer_colours = {}
                    drop_down = {}
                    sub_list = []
                    significance = [True for i in inner_kmers]

                    for b,inner in enumerate(inner_kmers):
                        relative_start = inner[0]-start
                        inner_colours = graded_dict[inner[2]]
                        for index in range(len(inner[2])):
                            relative_index = index+relative_start
                            kmer_colours[relative_index]=inner_colours
                            drop_down[relative_index]=b
                        
                        nodeID = (inner[0],inner[1],inner[2])


                        submenu='<br><span style="font-size:14px;line-height:0.5;padding: 12px 16px;color:black">Depth:'+str(Depth_of_kmers[nodeID[2]][0])+' ('+Depth_of_kmers[nodeID[2]][1]+')</span><br>'
                        significant = True
                        if stats_dict:
                            if float(stats_dict[nodeID[2]][1])>0.05 or float(stats_dict[nodeID[2]][3])>0.05 or float(stats_dict[nodeID[2]][2])>0.05 or float(stats_dict[nodeID[2]][4])>0.05:
                                significant = False
                            submenu+='<span style="font-size:14px;line-height:0.5;padding: 12px 16px;color:black">E<sub>i</sub>-value:'+stats_dict[nodeID[2]][1]+', P<sub>i</sub>-value:'+stats_dict[nodeID[2]][2]+'</span>'
                            submenu+='<br><span style="font-size:14px;line-height:0.5;padding: 12px 16px;color:black">E<sub>r</sub>-value:'+stats_dict[nodeID[2]][3]+', P<sub>r</sub>-value:'+stats_dict[nodeID[2]][4]+'</span><br>'
                        significance[b]=significant
      
                        if 'eCLIP' in Annotations: #Expand this
                            matches_to_eCLIP = False
                            try:
                                if 'BLAT' in Annotations['eCLIP']:
                                    if Annotations['eCLIP']['BLAT']:
                                        kmer_annotations = Annotations['eCLIP']['BLAT'][nodeID]
                                        if kmer_annotations:
                                            included=[]
                                            submenu+='<span style="font-size:18px;line-height:1.0;padding: 12px 16px;color:DarkBlue">eCLIP MATCHES</span>'
                                            matches_to_eCLIP = True
                                            for annotation in kmer_annotations:
                                                anno = annotation[0]
                                                if anno[3] not in included:
                                                    included.append(anno[3])
                                                    submenu+= '<a style="color:DarkBlue" href="eCLIP_results_BLAT.html#'+str(nodeID)+'">&#9654;'+anno[3]+' (bg='+str(anno[8])+'%)</a>'
                     
                                if 'BED' in Annotations['eCLIP']:
                                    if Annotations['eCLIP']['BED']:
                                        kmer_annotations = Annotations['eCLIP']['BED'][nodeID]
                                        if kmer_annotations:
                                            included=[]
                                            submenu+='<span style="font-size:18px;line-height:1.0;padding: 12px 16px;color:DarkBlue">eCLIP MATCHES</span>'
                                            matches_to_eCLIP = True
                                            for annotation in kmer_annotations:
                                                anno = annotation[0]
                                                if anno[3] not in included:
                                                    included.append(anno[3])                                        
                                                    submenu+= '<a style="color:DarkBlue" href="eCLIP_results_BED.html#'+str(nodeID)+'">&#9654;'+anno[3]+' (bg='+str(anno[8])+'%)</a>'
                            except Exception as e:
                                print str(e)

                            if not matches_to_eCLIP:
                                submenu+='<a style="color:DarkBlue">No matches to eCLIP Data</a>'
                    
                        if 'TargetScan' in Annotations:

                            kmer_annotations = Annotations['TargetScan'][nodeID]
                            if kmer_annotations:
                                submenu+='<span style="font-size:18px;line-height:1.0;padding: 12px 16px;color:DarkRed">MATCHES To TargetScan</span>'
                                for anno in kmer_annotations:
                                    seedID = anno.split(':')[1].strip()
                                    submenu+='<a style="color:DarkRed" href="miRNA_Matches.html#'+seedID+'">&#9654; '+anno+'</a>'
                            else:
                                submenu+='<a style="color:DarkRed">No matches to TargetScan</a>'
                        sub_list.append(submenu)


                    '''for inner in inner_kmers:

                        relative_start = inner[0]-start
                        inner_colours = graded_dict[inner[2]]
                        for index in range(len(inner[2])):
                            relative_index = index+relative_start
                            kmer_colours[relative_index]=inner_colours
                    marked_kmer = ''
                    marked_kmer_graded = ''
                    for index,base in enumerate(kmer):
                        kcolor = kmer_colours[index][0]
                        fc = kmer_colours[index][1]
                        marked_kmer_graded+='<mark style="background-color:'+kcolor+';color:'+fc+'">'+base+'</mark>'

                    kcolor = colours_dict[kmer][0]
                    fc = colours_dict[kmer][1]
                    marked_kmer+='<mark style="background-color:'+kcolor+';color:'+fc+'">'+kmer+'</mark>'
                    #print inner'''


                    seg_output = '<div class="dropdown"><span class="dropbtn">'
                    marked_kmer = ''

                    for index,base in enumerate(kmer):
                        kcolor = kmer_colours[index][0]
                        fc = kmer_colours[index][1]
                        si = drop_down[index]
                        submenu = sub_list[si]
                        i_kmer = inner_kmers[si][2]
                        significant = significance[si]
                        if significant:
                            i_kmer = i_kmer.upper()
                            base = base.upper() 
                        else:
                            i_kmer = i_kmer.lower()
                            base = base.lower() 
                        marked_kmer+='<mark style="background-color:'+kcolor+';color:'+fc+'">'+base+'</mark>'



                        if index<len(kmer)-1:
                            if si!=drop_down[index+1]:
                                seg_output+=marked_kmer+'</span><div class="dropdown-content"><span style="font-size:18px;line-height:1.0;padding: 12px 16px;color:black">'+i_kmer+'</span>'+submenu+'</div></div><div class="dropdown"><span class="dropbtn">'
                                marked_kmer = ''
                    seg_output+=marked_kmer+'</span><div class="dropdown-content"><span style="font-size:18px;line-height:1.0;padding: 12px 16px;color:black">'+i_kmer+'</span>'+submenu+'</div></div>'

                    if n==0:
                        end = node[1]
                        output+=' '*(6-len(str(start)))+str(start)+'-'+seg_output
                    else:
                        
                        basesBetween = (start-end)-1
                        #output+='-'+str(end)+'<span style="color:red">--('+str(basesBetween)+')--</span>'+str(start)+'-'+marked_kmer

                        output+='-'+str(end)+'<span style="color:grey">'+'--('+str(basesBetween)+')--</span>'+str(start)+'-'+seg_output
                        end = node[1]              
                            
                output+='-'+str(end)+'<br><br>'
        output+='</pre><hr style="height:15px;background-color:grey">'

    w= open(outdir+"/"+project_name+"/Html_Files/Modules.html",'w',0)
    w.write(output)
    w.close()


    return output


def write_miRNA_Tables(outdir,project_name,headers,miRNA_Tables,kmer_colour_dict):
    #Make a sub-directory to store html output file
    html_files_path = outdir+"/"+project_name+"/Html_Files"
    if not os.path.exists(html_files_path):
        os.mkdir(html_files_path)

    html_output ='<!DOCTYPE html><html><head><style>header{background-color:SkyBlue;text-align:center;}table, th, td {padding: 15px;  border: 1px solid black;border-collapse: collapse;}th {text-align: center;font-size:20px}</style><title>miRNA Matches</title></head><body>'
    html_output+='<header><h1 style="font-size:30px">Matches to miRNA Families<br><br>Retrieved from TargetScan </h1></header>'
    html_output+='<br><p style="color:red;font-size:20px">**SEED defined as positions 2-7 of mature miRNA sequence'
    html_output+='<br><span style="color:blue;font-size:20px">**MOTIF MATCHES that correspond to the reverse complement of miRNA SEEDS are displayed from 5`-3` relative to lincRNA sequence</span></p><br>'


    for seed in miRNA_Tables:
        html_output+='<table id="'+seed+'" style="width:95%"><tr><th style="text-align:left" colspan="6">Seed Matches to the miRNA '+miRNA_Tables[seed][0]+'</th></tr>'
        html_output+='<tr><th rowspan ="2">Seed</th><th rowspan ="2">Conservation</th><th rowspan ="2">Species</th><th colspan ="3">Matches</th></tr>'
        html_output+='<tr><th>Sequence</th><th>Motif</th><th>Type</th></tr>'

        matches = miRNA_Tables[seed][3]
        rspan = str(len(matches))
        html_output+='<tr><td rowspan = "'+rspan+'">'+seed+'</td><td rowspan = "'+rspan+'">'+miRNA_Tables[seed][1]+'</td><td rowspan = "'+rspan+'">'


        for species in miRNA_Tables[seed][2]:
            html_output+=species+'<br><br>'
        html_output+='</td>'
        
        match =matches[0]
        html_output+='<td>'+headers[match[0]-1].strip('>')+'</td>'
        colour = kmer_colour_dict[match[1]][0]
        fc = kmer_colour_dict[match[1]][1]

            
        html_output+='<td><mark style="background-color:'+colour+';color:'+fc+'">'+match[1]+'</mark></td>'
        html_output+='<td>'+match[2]+'</td>'
        html_output+='</tr>'

        for match in matches[1:]:
            html_output+='<tr><td>'+headers[match[0]-1].strip('>')+'</td>'
            colour = kmer_colour_dict[match[1]][0]
            fc = kmer_colour_dict[match[1]][1]

            
            html_output+='<td><mark style="background-color:'+colour+';color:'+fc+'">'+match[1]+'</mark></td>'
            html_output+='<td>'+match[2]+'</td>'
            html_output+='</tr>'
        html_output+='</table><br><br><br><br>\n'
   
    w= open(outdir+"/"+project_name+"/Html_Files/miRNA_Matches.html",'w',0)
    w.write(html_output)
    w.close()


def define_kmer_colours(outdir,project_name,kmers_dict,newcolours,uracil):
    colour_path = os.path.join(os.path.dirname(__file__), 'src', 'Kmer_colour_code.txt')

    kmer_colour_dict = {}
    if not newcolours:
        try:
            f = open(colour_path,'r')
            #f = open('src/Kmer_colour_code.txt','r')
            defined_colours = f.readlines()
            f.close()
            for colour in defined_colours:
                colour = colour.strip().split()
                kmer = colour[0]
                highlight = colour[1]
                font = colour[2]
                kmer_colour_dict[kmer] = (highlight,font)

            #Get all unique kmers (including combined kmers from overlaps) and define colours for each node
            all_kmer_sequences = []
            for layer in kmers_dict:
                for kmer in kmers_dict[layer]:
                    combined_kmer_seq = kmer[0][2]
                    seqments_of_combined_kmer = kmer[1]
                    all_kmer_sequences.append(combined_kmer_seq)
                    #for segment in seqments_of_combined_kmer:
                        #all_kmer_sequences.append(segment[2])

            all_kmer_sequences = list(dict.fromkeys(all_kmer_sequences))
            undefined = [k for k in all_kmer_sequences if k not in kmer_colour_dict]
            total_unique_kmers = len(undefined) 
            colour_set = generate_colour_set(total_unique_kmers)
        
            for i,k in enumerate(undefined):
                kmer_colour_dict[k] = colour_set[i]

            #Update all colours in src file
        
            #Write this colour code to a textfile in project folder
            colour_code_text = ''
            for kmer in kmer_colour_dict:
                highlight = kmer_colour_dict[kmer][0]
                font = kmer_colour_dict[kmer][1]
                colour_code_text+= kmer+' '+highlight+' '+font+'\n'
                if uracil:
                    ukmer = kmer.replace('U','T')
                    colour_code_text+= ukmer+' '+highlight+' '+font+'\n'
                else:
                    ukmer = kmer.replace('T','U')
                    colour_code_text+= ukmer+' '+highlight+' '+font+'\n'
            colour_code_text = colour_code_text.strip()
            w = open(colour_path,'w',0)
            #w = open('./src/Kmer_colour_code.txt','w',0)
            w.write(colour_code_text)
            w.close()

        except:
            print "Error! Colour code file not found or format is corrupted! Random color code generated"


            #Get all unique kmers (including combined kmers from overlaps) and define colours for each node
            all_kmer_sequences = []
            for layer in kmers_dict:
                for kmer in kmers_dict[layer]:
                    combined_kmer_seq = kmer[0][2]
                    seqments_of_combined_kmer = kmer[1]
                    all_kmer_sequences.append(combined_kmer_seq)
                    #for segment in seqments_of_combined_kmer:
                        #all_kmer_sequences.append(segment[2])

            all_kmer_sequences = list(dict.fromkeys(all_kmer_sequences))          
            total_unique_kmers = len(all_kmer_sequences) 
            colour_set = generate_colour_set(total_unique_kmers)
        

            for i,k in enumerate(all_kmer_sequences):
                kmer_colour_dict[k] = colour_set[i]

            #Write this colour code to a textfile in project folder
            colour_code_text = ''
            for kmer in kmer_colour_dict:
                highlight = kmer_colour_dict[kmer][0]
                font = kmer_colour_dict[kmer][1]
                colour_code_text+= kmer+' '+highlight+' '+font+'\n'
                if uracil:
                    ukmer = kmer.replace('U','T')
                    colour_code_text+= ukmer+' '+highlight+' '+font+'\n'
                else:
                    ukmer = kmer.replace('T','U')
                    colour_code_text+= ukmer+' '+highlight+' '+font+'\n'
            colour_code_text = colour_code_text.strip()
            w = open(outdir+"/"+project_name+"/Run_Files/Kmer_colour_code.txt",'w',0)
            w.write(colour_code_text)
            w.close()

    else:
        #Get all unique kmers (including combined kmers from overlaps) and define colours for each node
        all_kmer_sequences = []
        for layer in kmers_dict:
            for kmer in kmers_dict[layer]:
                combined_kmer_seq = kmer[0][2]
                seqments_of_combined_kmer = kmer[1]
                all_kmer_sequences.append(combined_kmer_seq)
                for segment in seqments_of_combined_kmer:
                    all_kmer_sequences.append(segment[2])

        all_kmer_sequences = list(dict.fromkeys(all_kmer_sequences))          
        total_unique_kmers = len(all_kmer_sequences) 
        colour_set = generate_colour_set(total_unique_kmers)
        

        for i,k in enumerate(all_kmer_sequences):
            kmer_colour_dict[k] = colour_set[i]

        #get colours of other kmers on src file - so only overwrite colours of current kmers
        previous_colours = {}
        kmers_to_save = []
        try:
            f = open(colour_path,'r')
            #f = open('src/Kmer_colour_code.txt','r')
            defined_colours = f.readlines()
            f.close()
            for colour in defined_colours:
                colour = colour.strip().split()
                kmer = colour[0]
                highlight = colour[1]
                font = colour[2]
                previous_colours[kmer] = (highlight,font)

            base_swap = []
            if uracil:
                base_swap.extend([k.replace('U','T') for k in kmer_colour_dict])
            else:
                base_swap.extend([k.replace('T','U') for k in kmer_colour_dict])

            kmers_to_save = [k for k in previous_colours if k not in kmer_colour_dict and k not in base_swap]

        except:
            print "src colour code file not found! New file generated"


        #Write this colour code to a textfile in project folder
        colour_code_text = ''
        for kmer in kmer_colour_dict:
            highlight = kmer_colour_dict[kmer][0]
            font = kmer_colour_dict[kmer][1]
            colour_code_text+= kmer+' '+highlight+' '+font+'\n'
            if uracil:
                ukmer = kmer.replace('U','T')
                colour_code_text+= ukmer+' '+highlight+' '+font+'\n'
            else:
                ukmer = kmer.replace('T','U')
                colour_code_text+= ukmer+' '+highlight+' '+font+'\n'

        for kmer in kmers_to_save:
            highlight = previous_colours[kmer][0]
            font = previous_colours[kmer][1]
            colour_code_text+= kmer+' '+highlight+' '+font+'\n'
        colour_code_text = colour_code_text.strip()
        w = open(colour_path,'w',0)
        #w = open('./src/Kmer_colour_code.txt','w',0)
        w.write(colour_code_text)
        w.close()

    return kmer_colour_dict


def define_colour_grade(total):


    grade_colours = {}
    #layers.sort(reverse=True)
    #total = len(layers)

    red = [(100,0,0)]
    yellow = [(255,0,0)]
    green = [(155,155,255)]
    aqua = [(0,0,250)]
    blue = [(0,0,100)]
    total = total-1 #min two layer depth

    if total>5:
        
        num_incr = (total-5)/4
        if (total-5)%4>0:
            num_incr+=1
        increment = 155/(num_incr+1)

        for i in range(1,num_incr+1):
            ri = red[i-1]
            ri = (ri[0]+increment,ri[1],ri[2])
            red.append(ri)

            yi = yellow[i-1]
            yi = (yi[0],yi[1]+increment,yi[2]+increment)
            yellow.append(yi)

            gi = green[i-1]
            gi = (gi[0]-increment,gi[1]-increment,gi[2]+increment)
            green.append(gi)

            ai = aqua[i-1]
            ai = (ai[0],ai[1],ai[2]-increment)
            aqua.append(ai)

    '''red = [(255,0,0)]
    yellow = [(255,255,0)]
    green = [(0,255,0)]
    aqua = [(0,255,255)]
    blue = [(0,0,255)]

    num_incr = (total-5)/4
    if (total-5)%4>0:
        num_incr+=1
    increment = 255/(num_incr+1)

    for i in range(1,num_incr+1):
        ri = red[i-1]
        ri = (ri[0],ri[1]+increment,ri[2])
        red.append(ri)

        yi = yellow[i-1]
        yi = (yi[0]-increment,yi[1],yi[2])
        yellow.append(yi)

        gi = green[i-1]
        gi = (gi[0],gi[1],gi[2]+increment)
        green.append(gi)

        ai = aqua[i-1]
        ai = (ai[0],ai[1]-increment,ai[2])
        aqua.append(ai)'''


    all_colours = []
    all_colours.extend(red)
    all_colours.extend(yellow)
    all_colours.extend(green)
    all_colours.extend(aqua)
    all_colours.extend(blue)

    index = range(total+1,1,-1)

    for i in range(total):
        rgb = all_colours[i]
        #l = layers[i]
        l = index[i]
        c = 'rgb('+str(rgb[0])+','+str(rgb[1])+','+str(rgb[2])+')'          
        #Determine white or black font colour
        fc = "black"
        if (rgb[0]+rgb[1]+rgb[2])/3<127:
            fc = "white"
        grade_colours[l]=(c,fc)
    return grade_colours



def write_kmers_in_seq_html(headers,sequences,intron_indices,kmers_dict,project_name,Annotations,stats_dict,Depth_of_kmers,kmer_colour_dict,fname):
    ###############################################################
    #WRITING HTML PAGES FOR COMBINED KMERS
    ###############################################################
    
    #Write kmers_in_seqs.html
    html_output ='<!DOCTYPE html><html><head><style>header{background-color:SkyBlue;text-align:center;} .dropbtn {font-size: 22px;border: none;cursor: pointer;}.dropdown {position: relative;display: inline-block;}.dropdown-content {display: none;position: absolute;background-color: rgb(255, 153, 153);white-space: nowrap; overflow:hidden;box-shadow: 0px 8px 50px 0px rgba(0,0,0,0.2);z-index: 1;}.dropdown-content a {line-height:1.0;font-size: 15px;padding: 12px 16px;text-decoration: none;display: block;}.dropdown:hover .dropdown-content {display: block;}.dropdown-content a:hover {background-color: yellow;}</style><title>KMERS_IN_SEQS</title></head><body>\n'
    html_output+='<header><h1 style="font-size:30px">KMERS IN SEQUENCES</h1><h2>Overlapping KMERS combined</h2></header>' 

    for i,seq in enumerate(sequences):
        head = headers[i]
        ls = str(len(seq))
        introns = intron_indices[i]
        seq_nodes = kmers_dict[i+1]
        seq_annotations = Annotations[i+1]
        seq = colour_kmers_in_seq_html(seq,seq_nodes,kmer_colour_dict,0,0,seq_annotations,introns,project_name,stats_dict,Depth_of_kmers) 
        html_output+='<h2>'+head+' ('+ls+' bases)</h2><pre style="font-size:22px;line-height:2.0;color:grey">'+seq+'</pre><hr>\n'

    html_output+= '</body></html>'  
                
    w= open("./"+project_name+"/Html_Files/"+fname,'w',0)
    w.write(html_output)
    w.close()
    ###############################################################


def write_kmers_in_blocks_html(headers,sequences,intron_indices,kmers_dict,project_name,Annotations,stats_dict,Depth_of_kmers,kmer_colour_dict,fname):
    #Write kmers_in_blocks.html
    html_output ='<!DOCTYPE html><html><head><style>header{background-color:SkyBlue;text-align:center;} .dropbtn {font-size: 22px;border: none;cursor: pointer;}.dropdown {position: relative;display: inline-block;}.dropdown-content {display: none;position: absolute;background-color: rgb(255, 153, 153);white-space: nowrap; overflow:hidden;box-shadow: 0px 8px 50px 0px rgba(0,0,0,0.2);z-index: 1;}.dropdown-content a {line-height:1.0;font-size: 18px;padding: 12px 16px;text-decoration: none;display: block;}.dropdown:hover .dropdown-content {display: block;}.dropdown-content a:hover {background-color: yellow;}</style><title>KMERS_IN_BLOCKS</title></head><body>'
    html_output+='<header><h1 style="font-size:30px">KMERS IN BLOCKS</h1><h2>Overlapping KMERS combined</header>'
    for i,seq in enumerate(sequences):
        ls = str(len(seq))
        seq='-'*len(seq)
        head = headers[i]
        introns = intron_indices[i]
        seq_nodes = kmers_dict[i+1]
        seq_annotations = Annotations[i+1]
        seq = colour_kmers_in_seq_html(seq,seq_nodes,kmer_colour_dict,0,1,seq_annotations,introns,project_name,stats_dict,Depth_of_kmers)
        html_output+='<h2>'+head+' ('+ls+' bases)</h2><pre style="font-size:22px;line-height:2.0;color:grey">'+seq+'</pre><hr>\n'
                
    html_output+= '</body></html>'        
    

    w= open("./"+project_name+"/Html_Files/"+fname,'w',0)
    w.write(html_output)
    w.close()

def write_graded_blocks_html(headers,headers_original,sequences,intron_indices,kmers_dict,outdir,project_name,Annotations,stats_dict,Depth_of_kmers,kmer_colour_dict,graded_colours,fname):
    html_output ='<!DOCTYPE html><html><head><style>header{background-color:SkyBlue;text-align:center;} .dropbtn {font-size: 22px;border: none;cursor: pointer;}.dropdown {position: relative;display: inline-block;}.dropdown-content {display: none;position: absolute;background-color: rgb(255, 153, 153);white-space: nowrap; overflow:hidden;box-shadow: 0px 8px 50px 0px rgba(0,0,0,0.2);z-index: 1;}.dropdown-content a {line-height:1.0;font-size: 15px;padding: 12px 16px;text-decoration: none;display: block;}.dropdown:hover .dropdown-content {display: block;}.dropdown-content a:hover {background-color: yellow;}\n'
    html_output+='</style><title>BLOCKS_CONSERVATION</title></head><body>'
    html_output+='<header><h1 style="font-size:30px">MOTIF CONSERVATION IN BLOCK DIAGRAMS</h1>'
    if fname == 'kmers_in_blocks_layer_conservation.html':
        html_output+='<h2>Motifs mapped to anchor sequence</h2>'

    html_key='<table style ="width:70%;margin-left:5%;margin-right:25%;border-collapse: collapse;table-layout: fixed"><tr>'
    levels = graded_colours.keys()
    #levels.sort(reverse=True)
    levels.sort()
    for level in levels:
        html_key+='<td style="background-color:'+graded_colours[level][0]+';color:'+graded_colours[level][0]+';text-align:center">|</td>'
    html_key+='</tr>'
    html_key+='<tr>'
    for level in levels:
        if level==levels[0]:
            html_key+='<td style="color:black;text-align:left">'+str(level)+'</td>'
        elif level==levels[-1]:
            html_key+='<td style="color:black;text-align:right">'+str(level)+'</td>'
        elif level==levels[len(levels)/2]:
            html_key+='<td style="color:black;text-align:center">'+str(level)+'</td>'
        else:
            html_key+='<td></td>'
    html_key+='</tr><tr><td style="color:black;text-align:center" colspan="'+str(len(levels))+'">Depth of motif conservation (number of species)</td></tr></table>'

    submenu = ''
    for i,seq in enumerate(sequences):
        submenu+='<a href="#'+str(i+1)+'" style="color:black;text-align:left">&#9654;'+headers_original[i]+'</a>'
    html_output+='<br><div class="dropdown"><span class="dropbtn"><mark style="background-color:LightGrey;color:black;font-size:20px;padding: 15px 15px 15px 15px;text-align: center;border: 3px solid black;">NAVIGATE &#9660</mark></span><div class="dropdown-content">'+submenu+'</div></div><br><br><br></header>\n'




    for i,seq in enumerate(sequences):
        ls = str(len(seq))
        seq='-'*len(seq)
        head = headers[i]
        introns = intron_indices[i]
        seq_nodes = kmers_dict[i+1]
        seq_annotations = Annotations[i+1]
        seq = colour_kmers_in_seq_html_2(seq,seq_nodes,kmer_colour_dict,1,1,seq_annotations,introns,project_name,stats_dict,Depth_of_kmers)
        html_output+='<h2 id="'+str(i+1)+'">'+head+' ('+ls+' bases)</h2><pre style="font-size:22px;line-height:2.0;color:grey">'+seq+'</pre>'
        html_output+=html_key+'<br><br><hr>\n'
                
    html_output+= '</body></html>'        
    

    w= open(outdir+"/"+project_name+"/Html_Files/"+fname,'w',0)
    w.write(html_output)
    w.close()


def write_kmers_graded_html(headers,headers_original,sequences,intron_indices,kmers_dict,outdir,project_name,Annotations,stats_dict,Depth_of_kmers,kmer_colour_dict,graded_colours,fname,select):
    html_output ='<!DOCTYPE html><html><head><style>header{background-color:SkyBlue;text-align:center;} .dropbtn {font-size: 22px;border: none;cursor: pointer;}.dropdown {position: relative;display: inline-block;}.dropdown-content {display: none;position: absolute;background-color: rgb(255, 153, 153);white-space: nowrap; overflow:hidden;box-shadow: 0px 8px 50px 0px rgba(0,0,0,0.2);z-index: 1;}.dropdown-content a {line-height:1.0;font-size: 15px;padding: 12px 16px;text-decoration: none;display: block;}.dropdown:hover .dropdown-content {display: block;}.dropdown-content a:hover {background-color: yellow;}\n'
    html_output+='</style><title>MOTIF CONSERVATION</title></head><body>'
    html_output+='<header><h1 style="font-size:30px">MOTIF CONSERVATION</h1>'
    if fname == 'kmers_in_seqs_layer_conservation.html':   
        html_output+='<h2>Motifs mapped to anchor sequence</h2>'
    elif fname == 'kmers_in_seqs_level_graded.html':
        html_output+='<h2>Motifs conserved to (and beyond) '+headers[select-1][0:25].strip('>')+' (depth:' +str(select)+')</h2>'  

    html_key='<table style ="width:70%;margin-left:5%;margin-right:25%;border-collapse: collapse;table-layout: fixed"><tr>'
    levels = graded_colours.keys()
    #levels.sort(reverse=True)
    levels.sort()
    for level in levels:
        html_key+='<td style="background-color:'+graded_colours[level][0]+';color:'+graded_colours[level][0]+';text-align:center">|</td>'
    html_key+='</tr>'
    html_key+='<tr>'
    for level in levels:
        if level==levels[0]:
            html_key+='<td style="color:black;text-align:left">'+str(level)+'</td>'
        elif level==levels[-1]:
            html_key+='<td style="color:black;text-align:right">'+str(level)+'</td>'
        elif level==levels[len(levels)/2]:
            html_key+='<td style="color:black;text-align:center">'+str(level)+'</td>'
        else:
            html_key+='<td></td>'
    html_key+='</tr><tr><td style="color:black;text-align:center" colspan="'+str(len(levels))+'">Depth of motif conservation (number of species)</td></tr></table>'

    submenu = ''
    for i,seq in enumerate(sequences):
        submenu+='<a href="#'+str(i+1)+'" style="color:black;text-align:left">&#9654;'+headers_original[i]+'</a>'
    html_output+='<br><div class="dropdown"><span class="dropbtn"><mark style="background-color:LightGrey;color:black;font-size:20px;padding: 15px 15px 15px 15px;text-align: center;border: 3px solid black;">NAVIGATE &#9660</mark></span><div class="dropdown-content">'+submenu+'</div></div><br><br><br></header>\n'



    for i,seq in enumerate(sequences):
        head = headers[i]
        ls = str(len(seq))
        introns = intron_indices[i]
        seq_nodes = kmers_dict[i+1]
        seq_annotations = Annotations[i+1]
        seq = colour_kmers_in_seq_html_2(seq,seq_nodes,kmer_colour_dict,1,0,seq_annotations,introns,project_name,stats_dict,Depth_of_kmers) 
        html_output+='<br><h2 id = "'+str(i+1)+'">'+head+' ('+ls+' bases)</h2><pre style="font-size:22px;line-height:2.0;color:grey">'+seq+'</pre>'
        html_output+=html_key+'<br><br><hr>\n'
                
    html_output+= '</body></html>'


    w= open(outdir+"/"+project_name+"/Html_Files/"+fname,'w',0)
    w.write(html_output)
    w.close()
    ###############################################################





def write_kmers_in_seq_overlap_html(headers,sequences,intron_indices,kmers_dict,outdir,project_name,Annotations,stats_dict,Depth_of_kmers,kmer_colour_dict,fname,select):
    ###############################################################
    #WRITING HTML PAGES TO DISPLAY OVERLAPPING KMERS
    ###############################################################
    
    

    html_output ='<!DOCTYPE html><html><head><style>header{background-color:SkyBlue;text-align:center;} .dropbtn {font-size: 22px;border: none;cursor: pointer;}.dropdown {position: relative;display: inline-block;}.dropdown-content {display: none;position: absolute;background-color: rgb(255, 153, 153);white-space: nowrap; overflow:hidden;box-shadow: 0px 8px 50px 0px rgba(0,0,0,0.2);z-index: 1;}.dropdown-content a {line-height:1.0;font-size: 15px;padding: 12px 16px;text-decoration: none;display: block;}.dropdown:hover .dropdown-content {display: block;}.dropdown-content a:hover {background-color: yellow;}\n'
    html_output+='</style><title>MOTIFS_IN_SEQS_OVERLAP</title></head><body>'
    html_output+='<header><h1 style="font-size:30px">MOTIFS IN SEQUENCES</h1>'
    if fname == 'kmers_in_seqs_level.html':
        html_output+='<h2>Motifs conserved to (and beyond) '+headers[select-1][0:25].strip('>')+' (depth:' +str(select)+')</h2>'

    submenu = ''
    for i,seq in enumerate(sequences):
        submenu+='<a href="#'+str(i+1)+'" style="color:black;text-align:left">&#9654;'+headers[i][0:25].strip('>')+' (depth:'+str(i+1)+')</a>'
    html_output+='<br><div class="dropdown"><span class="dropbtn"><mark style="background-color:LightGrey;color:black;font-size:20px;padding: 15px 15px 15px 15px;text-align: center;border: 3px solid black;">NAVIGATE &#9660</mark></span><div class="dropdown-content">'+submenu+'</div></div><br><br><br></header>\n'

    link = fname.split('.html')[0]+'_graded.html'
    html_output+='<br><a href='+link+' style="color:blue;">Coloured by Conservation >></a><br>'

    for i,seq in enumerate(sequences):
        head = headers[i]
        ls = str(len(seq))
        introns = intron_indices[i]
        seq_nodes = kmers_dict[i+1]
        seq_annotations = Annotations[i+1]
        seq = colour_kmers_in_seq_html_2(seq,seq_nodes,kmer_colour_dict,1,0,seq_annotations,introns,project_name,stats_dict,Depth_of_kmers) 
        html_output+='<h2 id="'+str(i+1)+'">'+head+' ('+ls+' bases)</h2><pre style="font-size:22px;line-height:2.0;color:grey">'+seq+'</pre><hr>\n'
                
    html_output+= '</body></html>'


    w= open(outdir+"/"+project_name+"/Html_Files/"+fname,'w',0)
    w.write(html_output)
    w.close()
    ###############################################################


def write_kmers_in_blocks_overlap_html(headers,sequences,intron_indices,kmers_dict,outdir,project_name,Annotations,stats_dict,Depth_of_kmers,kmer_colour_dict,fname,select):
    #Write kmers_in_blocks_overlap.html
    #html_output ='<!DOCTYPE html><html><head><style>header{background-color:SkyBlue;text-align:center;} .dropbtn {font-size: 22px;border: none;cursor: pointer;}.dropdown {position: relative;display: inline-block;}.dropdown-content {display: none;position: absolute;background-color: rgb(255, 153, 153)  ;white-space: nowrap; overflow:hidden;box-shadow: 0px 8px 50px 0px rgba(0,0,0,0.2);z-index: 1;}.dropdown-content a {line-height:1.0;font-size: 18px;padding: 12px 16px;text-decoration: none;display: block;}.dropdown:hover .dropdown-content {display: block;}.dropdown-content a:hover {background-color: yellow;}</style><title>KMERS_IN_BLOCKS</title></head><body>'
    #html_output+='<header><h1 style="font-size:30px">KMERS IN BLOCKS</h1><h2>All Overlapping KMERS Shown</h2>'


    html_output ='<!DOCTYPE html><html><head><style>header{background-color:SkyBlue;text-align:center;} .dropbtn {font-size: 22px;border: none;cursor: pointer;}.dropdown {position: relative;display: inline-block;}.dropdown-content {display: none;position: absolute;background-color: rgb(255, 153, 153);white-space: nowrap; overflow:hidden;box-shadow: 0px 8px 50px 0px rgba(0,0,0,0.2);z-index: 1;}.dropdown-content a {line-height:1.0;font-size: 15px;padding: 12px 16px;text-decoration: none;display: block;}.dropdown:hover .dropdown-content {display: block;}.dropdown-content a:hover {background-color: yellow;}\n'
    html_output+='</style><title>MOTIFS IN BLOCKS</title></head><body>'
    html_output+='<header><h1 style="font-size:30px">MOTIFS IN BLOCK DIAGRAMS</h1>'
    if fname == 'kmers_in_blocks_level.html':
        html_output+='<h2>Motifs conserved to (and beyond) '+headers[select-1][0:25].strip('>')+' (depth:' +str(select)+')</h2>'    


    submenu = ''
    for i,seq in enumerate(sequences):
        submenu+='<a href="#'+str(i+1)+'" style="color:black;text-align:left">&#9654;'+headers[i][0:25].strip('>')+' (depth:'+str(i+1)+')</a>'
    html_output+='<br><div class="dropdown"><span class="dropbtn"><mark style="background-color:LightGrey;color:black;font-size:20px;padding: 15px 15px 15px 15px;text-align: center;border: 3px solid black;">NAVIGATE &#9660</mark></span><div class="dropdown-content">'+submenu+'</div></div><br><br><br></header>\n'


    for i,seq in enumerate(sequences):
        ls = str(len(seq))
        seq='-'*len(seq)
        head = headers[i]
        introns = intron_indices[i]
        seq_nodes = kmers_dict[i+1]
        seq_annotations = Annotations[i+1]
        seq = colour_kmers_in_seq_html_2(seq,seq_nodes,kmer_colour_dict,1,1,seq_annotations,introns,project_name,stats_dict,Depth_of_kmers)
        html_output+='<h2 id="'+str(i+1)+'">'+head+' ('+ls+' bases)</h2><pre style="font-size:22px;line-height:2.0;color:grey">'+seq+'</pre><hr>\n'
                
    html_output+= '</body></html>'
     
    
    w= open(outdir+"/"+project_name+"/Html_Files/"+fname,'w',0)
    w.write(html_output)
    w.close()





def write_html(headers,sequences,intron_indices,kmers_dict,project_name,Annotations,Modules,Modules_in_sequences,kmer_colour_dict,deep):

    ###############################################################
    if Modules and Modules_in_sequences:
        #Write modules.html
        html_output ='<!DOCTYPE html><html><head><style>header{background-color:SkyBlue;text-align:center;}</style><title>MODULES</title></head><body>'
        html_output+='<header><h1 style="font-size:30px">MODULES</h1><h2>start-KMER-end<span style="color:red">--(number of bases in between)--</span>start-KMER-end</h2><h2>Order is random (to be changed soon): Overlapping KMERS combined<br>As is, only consecutive uniterrupted kmers have been considered</h2></header>'        
    
    
        html_output+=write_modules(Modules,Modules_in_sequences,kmers_dict,headers,kmer_colour_dict,0,Annotations,introns)

        html_output+= '</body></html>'

        w= open("./"+project_name+"/Html_Files/modules.html",'w',0)
        w.write(html_output)
        w.close()
    ###############################################################

    
    ###############################################################
    if Modules and Modules_in_sequences:
        #Write modules_overlap.html
        html_output ='<!DOCTYPE html><html><head><style>header{background-color:SkyBlue;text-align:center;}</style><title>MODULES</title></head><body>'
        html_output+='<header><h1 style="font-size:30px">MODULES</h1><h2>KMER<span style="color:red">--(number of bases in between)--</span>KMER</h2></header>'        
       
        html_output+=write_modules(Modules,Modules_in_sequences,kmers_dict,headers,kmer_colour_dict,1,Annotations,introns)

        html_output+= '</body></html>'

        w= open("./"+project_name+"/Html_Files/modules_overlap.html",'w',0)
        w.write(html_output)
        w.close()
    ###############################################################


def write_lncLOOM_txt(all_kmers,MainGraph,MainGraphLevels,number_of_layers,kmers_depth,LOOM5,LOOM3,LOOM5_Levels,LOOM3_Levels,details5,details3,sequences,headers,outdir,project_name,stats_dict,Annotations):

    output= ["==========================================================================================================\n"]
    output.append("==========================================================================================================\n")
    output.append('LncLOOM Results\n')
    output.append("==========================================================================================================\n")
    output.append("==========================================================================================================\n\n")
    output.append('All Sequences in Dataset (Main Graph Calculation):\n')
    output.append('Name'+' '*51+'Depth'+' '*10+'Length\n')
    output.append("----------------------------------------------------------------------------------------------------------\n")
    for i,head in enumerate(headers):
        output.append(head[0:50]+' '*(55-len(head[0:50]))+str(i+1)+' '*(15-len(str(i+1)))+str(len(sequences[i]))+'\n')
    output.append("==========================================================================================================\n\n")
    output.append("Sequences with extended 5' region (5' Graph Calculation):\n")
    output.append('Name'+' '*51+'Depth'+' '*10+'Length of Extension\n')
    output.append("----------------------------------------------------------------------------------------------------------\n")
    if details5:        
        for seq in details5:
            output.append(headers[seq[0]-1][0:50]+' '*(55-len(headers[seq[0]-1][0:50]))+str(seq[0])+' '*(15-len(str(seq[0])))+str(seq[1])+"\n")
    else:
        output.append('None\n\n')
    output.append("==========================================================================================================\n\n") 
    output.append("Sequences with extended 3' region (3' Graph Calculation):\n")
    output.append('Name'+' '*51+'Depth'+' '*10+'Length of Extension\n')
    output.append("----------------------------------------------------------------------------------------------------------\n")
    if details3:        
        for seq in details3:
            output.append(headers[seq[0]-1][0:50]+' '*(55-len(headers[seq[0]-1][0:50]))+str(seq[0])+' '*(15-len(str(seq[0])))+str(seq[1])+"\n")
    else:
        output.append('None\n\n')

   
    output.append("==========================================================================================================\n\n")    
    output.append('Main Graph Results (Section 1 of 4)\n\n\n')


    neighbourhoods = SM.get_neighbourhoods(MainGraph)


    sort_by_depth = [(nh,neighbourhoods[nh]['Depth'],neighbourhoods[nh][1][0][0]) for nh in neighbourhoods]
    sort_by_depth.sort(key=lambda x:(-x[1],x[2]))

    kmers_dict = SM.get_kmers_per_layer(MainGraphLevels,number_of_layers)
    dict_by_kmer = {}    
    for kmer in kmers_depth:
        dict_by_kmer[kmer]={}
        for depth in range(1,kmers_depth[kmer][0]+1):
            dict_by_kmer[kmer][depth]=[]


    for level in kmers_dict:
        inner_kmers = [x[1] for x in kmers_dict[level]]
        for inner_list in inner_kmers:
            for k in inner_list:
                dict_by_kmer[k[2]][level].append(k)

    top_kmers = {x[0]:x[1] for x in kmers_dict[1]}

    nh_id = 1

    for motifnh in sort_by_depth:

        motif_id_set = motifnh[0]
        depth = motifnh[1]

        kmer_set = top_kmers[(neighbourhoods[motif_id_set][1][0][0][0],neighbourhoods[motif_id_set][1][0][0][1],neighbourhoods[motif_id_set][1][0][0][2])]
        reference = neighbourhoods[motif_id_set][1][0][0][2]
        ref_end = neighbourhoods[motif_id_set][1][0][0][1]
        for ref in neighbourhoods[motif_id_set][1][0][1:]:
            inbetween = '-'*(ref[0]-ref_end-1)
            reference+= inbetween+ref[2]
            ref_end = ref[1]
            kmer_set.extend(top_kmers[(ref[0],ref[1],ref[2])])

        ref_end = len(reference)        

        #output.append("***************************************************************************************************************************************************"+'*'*(ref_end*3)+"\n")
        #output.append('Motif Neighborhood '+str(nh_id)+'   Depth:'+str(depth)+'\n')    
        #output.append("---------------------------------------------------------------------------------------------------------------------------------------------------"+'-'*(ref_end*3)+"\n")
        #output.append("Species"+' '*48+"Start"+' '*10+"End"+' '*12+"Width"+' '*10+"Depth"+' '*10+"Conserved Sites"+' '*(ref_end+10-15)+"Motif Neighborhood\n")
        #output.append("---------------------------------------------------------------------------------------------------------------------------------------------------"+'-'*(ref_end*3)+"\n")

        layers = [x for x in neighbourhoods[motif_id_set].keys() if x!='Depth']
        layers.sort()
        #determine longest region


        longest_region = 0
        layer_references = {}
        for layer in layers:
            limit1 = neighbourhoods[motif_id_set][layer][1]
            limit2 = neighbourhoods[motif_id_set][layer][2]
            extract = sequences[layer-1][limit1-1:limit2]  
            layer_references[layer] = extract.lower()
            if len(extract)>longest_region:
                longest_region = len(extract)

        output.append("***************************************************************************************************************************************************"+'*'*(longest_region*3)+"\n")
        output.append('Motif Neighborhood '+str(nh_id)+'   Depth:'+str(depth)+'\n')    
        output.append("___________________________________________________________________________________________________________________________________________________"+'_'*(longest_region*3)+"\n")
        output.append("\nSpecies"+' '*48+"Start"+' '*10+"End"+' '*12+"Width"+' '*10+"Depth"+' '*10+"Conserved Sites"+' '*(longest_region+10-15)+"Motif Neighborhood\n")
        output.append("___________________________________________________________________________________________________________________________________________________"+'_'*(longest_region*3)+"\n\n")

        for layer in layers:

            sites_for_layer = neighbourhoods[motif_id_set][layer][0]
            #print neighbourhoods[motif_id_set]
            species = headers[layer-1][0:50]
            n1 = sites_for_layer[0]
            start = n1[0]
            #end1 = n1[1]

            n2 = sites_for_layer[-1]
            #start2 = n2[0]
            end = n2[1]
            
            kmer = layer_references[layer]
            ref_index = reference.index(n1[2])
            if len(sites_for_layer)>1:
                last_ref_index = reference.rindex(n2[2])
            else:
                last_ref_index = ref_index

            ref_end_n = ref_end-(last_ref_index+len(n2[2]))

            if start>ref_index:
                ad1 = sequences[layer-1][start-ref_index-1:start-1]
                ad2 = sequences[layer-1][end:end+ref_end_n]
                n_start = start-ref_index
            else:
                ad1 = ' '*(ref_index-start+1)+sequences[layer-1][0:start-1]
                ad2 = sequences[layer-1][end:end+ref_end_n]
                n_start = 1

            if end+ref_end_n>len(sequences[layer-1]):
                n_end = len(sequences[layer-1])
            else:
                n_end=end+ref_end_n

            width = n_end-n_start+1
            #width = len(kmer)
            #kmer_alone = ' '*len(kmer)
            
            #kmer_alone = ' '*ref_index+n[2]+' '*ref_end_n
            kmer_alone = '-'*len(kmer)
            kindex = -1
            kmer_upper = kmer.upper()
            for site in sites_for_layer:
                k = site[2].upper()
                kindex = kmer_upper.index(k,kindex+1)
                kend = kindex+len(k)

                kmer = kmer[0:kindex]+k+kmer[kend:] 
                kmer_alone = kmer_alone[0:kindex]+k+kmer_alone[kend:]
            kmer = ad1.lower()+kmer+ad2.lower()

            output.append(species+' '*(55-len(species)))
            output.append(str(n_start)+' '*(15-len(str(n_start))))
            output.append(str(n_end)+' '*(15-len(str(n_end))))
            output.append(str(width)+' '*(15-len(str(width))))
            output.append(str(layer)+' '*(15-len(str(layer))))
            output.append(kmer_alone+' '*(longest_region+10-len(kmer_alone)))
            output.append(kmer+"\n")
            #output.append(additional)



        output.append("___________________________________________________________________________________________________________________________________________________"+'_'*(longest_region*3)+"\n")
        output.append('\nIndividual Motif Sites in Neighborhood '+str(nh_id)+':\n')
        output.append('----------------------------------------------------------------\n\n')

        motif_id =1
        kmer_set.sort(key = lambda x: (-x[3],x[0]))
        kmer_checked = []
        for item in kmer_set:
            if item[2] in kmer_checked:
                continue
            kmer = item[2]
            kmer_checked.append(kmer)
            level = item[3]
            stats_text = '\nStatistical Significance Not Calculated\n'
            if stats_dict:
                stats = stats_dict[kmer]
                stats_text = '\nE(i)-value='+str(stats[1])+'    P(i)-value='+str(stats[2])+'    E(r)-value='+str(stats[3])+'    E(r)-value='+str(stats[4])+'\n'
            output.append("--------------------------------------------------------------------------------------------------------------------------------------------------"+'-'*(longest_region*3)+"\n")
            output.append('Motif '+str(nh_id)+'.'+str(motif_id)+'   Depth:'+str(level)+'\n')
            output.append(stats_text)
            output.append("--------------------------------------------------------------------------------------------------------------------------------------------------"+'-'*(longest_region*3)+"\n")
            output.append("Species"+' '*48+"Start"+' '*10+"End"+' '*12+"Width"+' '*10+"Depth"+' '*10+"          -Site-          \n")

            output.append("--------------------------------------------------------------------------------------------------------------------------------------------------"+'-'*(longest_region*3)+"\n")

            targetscan = ''
            eclip_blat = ''
            eclip_bed = ''
            site_for_all = dict_by_kmer[kmer]
            layers = [x for x in site_for_all.keys() if site_for_all[x] and x in neighbourhoods[motif_id_set]]
            layers.sort()
            for layer in layers:
                limit1 = neighbourhoods[motif_id_set][layer][1]
                limit2 = neighbourhoods[motif_id_set][layer][2]
                nodes = [x for x in site_for_all[layer] if x[0]>=limit1 and x[1]<=limit2]
                for n in nodes:
                    species = headers[layer-1][0:50]
                    start = n[0]
                    end = n[1]
                    if start>10:
                        ad1 = sequences[layer-1][start-11:start-1]
                        ad2 = sequences[layer-1][end:end+10]
                    else:
                        ad1 = sequences[layer-1][0:start-1]
                        ad2 = sequences[layer-1][end:end+10]
                    output.append(species+' '*(55-len(species)))
                    output.append(str(start)+' '*(15-len(str(start))))
                    output.append(str(end)+' '*(15-len(str(end))))
                    output.append(str(len(kmer))+' '*(15-len(str(len(kmer)))))
                    output.append(str(layer)+' '*(15-len(str(layer))))
                    output.append(ad1.lower()+' '*(10-len(ad1))+'-'+kmer+'-'+ad2.lower()+"\n")
                    if 'TargetScan' in Annotations[layer]:
                        anno=''
                        for match in Annotations[layer]['TargetScan'][(n[0],n[1],n[2])]:
                            anno+=match.split(':')[0]+','
                        if anno!='':
                            anno=species+':    '+anno+'\n'
                            targetscan+=anno

                    if 'eCLIP' in Annotations[layer]:
                        if 'BLAT' in Annotations[layer]['eCLIP']:
                            anno=''
                            for match in Annotations[layer]['eCLIP']['BLAT'][(n[0],n[1],n[2])]:
                                anno+=match[0][3]+','
                            if anno!='':
                                anno=species+':    '+anno+'\n'
                                eclip_blat+=anno

                        if 'BED' in Annotations[layer]['eCLIP']:
                            anno=''
                            for match in Annotations[layer]['eCLIP']['BED'][(n[0],n[1],n[2])]:
                                anno+=match[0][3]+','
                            if anno!='':
                                anno=species+':    '+anno+'\n'
                                eclip_bed+=anno

            if eclip_bed =='':
                eclip_bed = ' None\n'
            if eclip_blat =='':
                eclip_blat = ' None\n'
            if targetscan =='':
                targetscan = ' None\n'

            output.append("--------------------------------------------------------------------------------------------------------------------------------------------------"+'-'*(longest_region*3)+"\n")
            output.append('Motif '+str(nh_id)+'.'+str(motif_id)+' ('+kmer+') Annotations:\n')
            output.append('----------------------------------------------\n')
            output.append('TargetScan Matches:\n')
            output.append(targetscan)
            output.append('\n\n')
            output.append('eCLIP determined binding proteins (based on BLAT alignment):\n')
            output.append(eclip_blat)
            output.append('\n\n')
            output.append('eCLIP determined binding proteins (based on BED File):\n')
            output.append(eclip_bed)
            output.append("--------------------------------------------------------------------------------------------------------------------------------------------------"+'-'*(longest_region*3)+"\n\n\n\n")

            motif_id+=1
        nh_id+=1


    output.append("="*300+"\n")
    output.append("5' Graph Results (Section 2 of 4)\n\n\n")

    neighbourhoods = SM.get_neighbourhoods(LOOM5)


    sort_by_depth = [(nh,neighbourhoods[nh]['Depth'],neighbourhoods[nh][1][0][0]) for nh in neighbourhoods]
    sort_by_depth.sort(key=lambda x:(-x[1],x[2]))

    kmers_dict = SM.get_kmers_per_layer(LOOM5_Levels,number_of_layers)
    dict_by_kmer = {}    
    for kmer in kmers_depth:
        dict_by_kmer[kmer]={}
        for depth in range(1,kmers_depth[kmer][0]+1):
            dict_by_kmer[kmer][depth]=[]


    for level in kmers_dict:
        inner_kmers = [x[1] for x in kmers_dict[level]]
        for inner_list in inner_kmers:
            for k in inner_list:
                dict_by_kmer[k[2]][level].append(k)

    top_kmers = {x[0]:x[1] for x in kmers_dict[1]}
    #for nh in neighbourhoods:
        #print neighbourhoods[nh][1]

     
    nh_id = 1

    for motifnh in sort_by_depth:

        motif_id_set = motifnh[0]
        depth = motifnh[1]

        kmer_set = top_kmers[(neighbourhoods[motif_id_set][1][0][0][0],neighbourhoods[motif_id_set][1][0][0][1],neighbourhoods[motif_id_set][1][0][0][2])]
        reference = neighbourhoods[motif_id_set][1][0][0][2]
        ref_end = neighbourhoods[motif_id_set][1][0][0][1]
        for ref in neighbourhoods[motif_id_set][1][0][1:]:
            inbetween = '-'*(ref[0]-ref_end-1)
            reference+= inbetween+ref[2]
            ref_end = ref[1]
            kmer_set.extend(top_kmers[(ref[0],ref[1],ref[2])])

        ref_end = len(reference) 

        #output.append("***************************************************************************************************************************************************"+'*'*(ref_end*3)+"\n")
        #output.append('Motif Neighborhood '+str(nh_id)+'   Depth:'+str(depth)+'\n')    
        #output.append("---------------------------------------------------------------------------------------------------------------------------------------------------"+'-'*(ref_end*3)+"\n")
        #output.append("Species"+' '*48+"Start"+' '*10+"End"+' '*12+"Width"+' '*10+"Depth"+' '*10+"Conserved Sites"+' '*(ref_end+10-15)+"Motif Neighborhood\n")
        #output.append("---------------------------------------------------------------------------------------------------------------------------------------------------"+'-'*(ref_end*3)+"\n")

        layers = [x for x in neighbourhoods[motif_id_set].keys() if x!='Depth']
        layers.sort()
        #determine longest region


        longest_region = 0
        layer_references = {}
        for layer in layers:
            limit1 = neighbourhoods[motif_id_set][layer][1]
            limit2 = neighbourhoods[motif_id_set][layer][2]
            extract = sequences[layer-1][limit1-1:limit2]  
            layer_references[layer] = extract.lower()
            if len(extract)>longest_region:
                longest_region = len(extract)

        output.append("***************************************************************************************************************************************************"+'*'*(longest_region*3)+"\n")
        output.append('Motif Neighborhood '+str(nh_id)+'   Depth:'+str(depth)+'\n')    
        output.append("___________________________________________________________________________________________________________________________________________________"+'_'*(longest_region*3)+"\n")
        output.append("\nSpecies"+' '*48+"Start"+' '*10+"End"+' '*12+"Width"+' '*10+"Depth"+' '*10+"Conserved Sites"+' '*(longest_region+10-15)+"Motif Neighborhood\n")
        output.append("___________________________________________________________________________________________________________________________________________________"+'_'*(longest_region*3)+"\n\n")
        for layer in layers:

            sites_for_layer = neighbourhoods[motif_id_set][layer][0]
            species = headers[layer-1][0:50]
            n1 = sites_for_layer[0]
            start = n1[0]
            #end1 = n1[1]

            n2 = sites_for_layer[-1]
            #start2 = n2[0]
            end = n2[1]
            
            kmer = layer_references[layer]

            ref_index = reference.index(n1[2])
            last_ref_index = reference.rindex(n2[2])
            ref_end_n = ref_end-(last_ref_index+len(n2[2]))

            if start>ref_index:
                ad1 = sequences[layer-1][start-ref_index-1:start-1]
                ad2 = sequences[layer-1][end:end+ref_end_n]
                n_start = start-ref_index
            else:
                ad1 = ' '*(ref_index-start+1)+sequences[layer-1][0:start-1]
                ad2 = sequences[layer-1][end:end+ref_end_n]
                n_start = 1

            if end+ref_end_n>len(sequences[layer-1]):
                n_end = len(sequences[layer-1])
            else:
                n_end=end+ref_end_n

            width = n_end-n_start+1
            #width = len(kmer)
            kmer_alone = ' '*len(kmer)
            
            #kmer_alone = ' '*ref_index+n[2]+' '*ref_end_n
            kmer_alone = '-'*len(kmer)
            kindex = -1
            kmer_upper = kmer.upper()
            for site in sites_for_layer:
                k = site[2].upper()
                kindex = kmer_upper.index(k,kindex+1)
                kend = kindex+len(k)

                kmer = kmer[0:kindex]+k+kmer[kend:] 
                kmer_alone = kmer_alone[0:kindex]+k+kmer_alone[kend:]
            kmer = ad1.lower()+kmer+ad2.lower()

            output.append(species+' '*(55-len(species)))
            output.append(str(n_start)+' '*(15-len(str(n_start))))
            output.append(str(n_end)+' '*(15-len(str(n_end))))
            output.append(str(width)+' '*(15-len(str(width))))
            output.append(str(layer)+' '*(15-len(str(layer))))
            output.append(kmer_alone+' '*(longest_region+10-len(kmer_alone)))
            output.append(kmer+"\n")
            #output.append(additional)



        output.append("___________________________________________________________________________________________________________________________________________________"+'_'*(longest_region*3)+"\n")
        output.append('\nIndividual Motif Sites in Neighborhood '+str(nh_id)+':\n')
        output.append('----------------------------------------------------------------\n\n')
        motif_id =1
        kmer_set.sort(key = lambda x: (-x[3],x[0]))
        kmer_checked = []
        for item in kmer_set:
            if item[2] in kmer_checked:
                continue
            kmer = item[2]
            kmer_checked.append(kmer)
            level = item[3]
            stats_text = '\nStatistical Significance Not Calculated\n'
            if stats_dict:
                stats = stats_dict[kmer]
                stats_text = '\nE(i)-value='+str(stats[1])+'    P(i)-value='+str(stats[2])+'    E(r)-value='+str(stats[3])+'    E(r)-value='+str(stats[4])+'\n'
            output.append("--------------------------------------------------------------------------------------------------------------------------------------------------"+'-'*(ref_end*3)+"\n")
            output.append('Motif '+str(nh_id)+'.'+str(motif_id)+'   Depth:'+str(level)+'\n')
            output.append(stats_text)
            output.append("--------------------------------------------------------------------------------------------------------------------------------------------------"+'-'*(ref_end*3)+"\n")
            output.append("Species"+' '*48+"Start"+' '*10+"End"+' '*12+"Width"+' '*10+"Depth"+' '*10+"          -Site-          \n")

            output.append("--------------------------------------------------------------------------------------------------------------------------------------------------"+'-'*(ref_end*3)+"\n")

            targetscan = ''
            eclip_blat = ''
            eclip_bed = ''
            site_for_all = dict_by_kmer[kmer]
            layers = [x for x in site_for_all.keys() if site_for_all[x] and x in neighbourhoods[motif_id_set]]
            layers.sort()
            for layer in layers:
                limit1 = neighbourhoods[motif_id_set][layer][1]
                limit2 = neighbourhoods[motif_id_set][layer][2]
                nodes = [x for x in site_for_all[layer] if x[0]>=limit1 and x[1]<=limit2]
                for n in nodes:
                    species = headers[layer-1][0:50]
                    start = n[0]
                    end = n[1]
                    if start>10:
                        ad1 = sequences[layer-1][start-11:start-1]
                        ad2 = sequences[layer-1][end:end+10]
                    else:
                        ad1 = sequences[layer-1][0:start-1]
                        ad2 = sequences[layer-1][end:end+10]
                    output.append(species+' '*(55-len(species)))
                    output.append(str(start)+' '*(15-len(str(start))))
                    output.append(str(end)+' '*(15-len(str(end))))
                    output.append(str(len(kmer))+' '*(15-len(str(len(kmer)))))
                    output.append(str(layer)+' '*(15-len(str(layer))))
                    output.append(ad1.lower()+' '*(10-len(ad1))+'-'+kmer+'-'+ad2.lower()+"\n")
                    if 'TargetScan' in Annotations[layer]:
                        anno=''
                        for match in Annotations[layer]['TargetScan'][(n[0],n[1],n[2])]:
                            anno+=match.split(':')[0]+','
                        if anno!='':
                            anno=species+':    '+anno+'\n'
                            targetscan+=anno

                    if 'eCLIP' in Annotations[layer]:
                        if 'BLAT' in Annotations[layer]['eCLIP']:
                            anno=''
                            for match in Annotations[layer]['eCLIP']['BLAT'][(n[0],n[1],n[2])]:
                                anno+=match[0][3]+','
                            if anno!='':
                                anno=species+':    '+anno+'\n'
                                eclip_blat+=anno

                        if 'BED' in Annotations[layer]['eCLIP']:
                            anno=''
                            for match in Annotations[layer]['eCLIP']['BED'][(n[0],n[1],n[2])]:
                                anno+=match[0][3]+','
                            if anno!='':
                                anno=species+':    '+anno+'\n'
                                eclip_bed+=anno

            if eclip_bed =='':
                eclip_bed = ' None\n'
            if eclip_blat =='':
                eclip_blat = ' None\n'
            if targetscan =='':
                targetscan = ' None\n'

            output.append("--------------------------------------------------------------------------------------------------------------------------------------------------"+'-'*(ref_end*3)+"\n")
            output.append('Motif '+str(nh_id)+'.'+str(motif_id)+' ('+kmer+') Annotations:\n')
            output.append('----------------------------------------------\n')
            output.append('TargetScan Matches:\n')
            output.append(targetscan)
            output.append('\n\n')
            output.append('eCLIP determined binding proteins (based on BLAT alignment):\n')
            output.append(eclip_blat)
            output.append('\n\n')
            output.append('eCLIP determined binding proteins (based on BED File):\n')
            output.append(eclip_bed)
            output.append("--------------------------------------------------------------------------------------------------------------------------------------------------"+'-'*(ref_end*3)+"\n\n\n\n")

            motif_id+=1
        nh_id+=1

    output.append("="*300+"\n")
    output.append("3' Graph Results (Section 3 of 4)\n\n\n")

    neighbourhoods = SM.get_neighbourhoods(LOOM3)

    sort_by_depth = [(nh,neighbourhoods[nh]['Depth'],neighbourhoods[nh][1][0][0]) for nh in neighbourhoods]
    sort_by_depth.sort(key=lambda x:(-x[1],x[2]))

    kmers_dict = SM.get_kmers_per_layer(LOOM3_Levels,number_of_layers)
    dict_by_kmer = {}    
    for kmer in kmers_depth:
        dict_by_kmer[kmer]={}
        for depth in range(1,kmers_depth[kmer][0]+1):
            dict_by_kmer[kmer][depth]=[]


    for level in kmers_dict:
        inner_kmers = [x[1] for x in kmers_dict[level]]
        for inner_list in inner_kmers:
            for k in inner_list:
                dict_by_kmer[k[2]][level].append(k)

    top_kmers = {x[0]:x[1] for x in kmers_dict[1]}
    #print top_kmers
    #for nh in neighbourhoods:
        #print neighbourhoods[nh][1]
     
    nh_id = 1

    for motifnh in sort_by_depth:

        motif_id_set = motifnh[0]
        depth = motifnh[1]

        kmer_set = top_kmers[(neighbourhoods[motif_id_set][1][0][0][0],neighbourhoods[motif_id_set][1][0][0][1],neighbourhoods[motif_id_set][1][0][0][2])]
        reference = neighbourhoods[motif_id_set][1][0][0][2]
        ref_end = neighbourhoods[motif_id_set][1][0][0][1]
        for ref in neighbourhoods[motif_id_set][1][0][1:]:
            inbetween = '-'*(ref[0]-ref_end-1)
            reference+= inbetween+ref[2]
            ref_end = ref[1]
            kmer_set.extend(top_kmers[(ref[0],ref[1],ref[2])])

        ref_end = len(reference) 

        #output.append("***************************************************************************************************************************************************"+'*'*(ref_end*3)+"\n")
        #output.append('Motif Neighborhood '+str(nh_id)+'   Depth:'+str(depth)+'\n')    
        #output.append("---------------------------------------------------------------------------------------------------------------------------------------------------"+'-'*(ref_end*3)+"\n")
        #output.append("Species"+' '*48+"Start"+' '*10+"End"+' '*12+"Width"+' '*10+"Depth"+' '*10+"Conserved Sites"+' '*(ref_end+10-15)+"Motif Neighborhood\n")
        #output.append("---------------------------------------------------------------------------------------------------------------------------------------------------"+'-'*(ref_end*3)+"\n")

        layers = [x for x in neighbourhoods[motif_id_set].keys() if x!='Depth']
        layers.sort()
        #determine longest region


        longest_region = 0
        layer_references = {}
        for layer in layers:
            limit1 = neighbourhoods[motif_id_set][layer][1]
            limit2 = neighbourhoods[motif_id_set][layer][2]
            extract = sequences[layer-1][limit1-1:limit2]  
            layer_references[layer] = extract.lower()
            if len(extract)>longest_region:
                longest_region = len(extract)

        output.append("***************************************************************************************************************************************************"+'*'*(longest_region*3)+"\n")
        output.append('Motif Neighborhood '+str(nh_id)+'   Depth:'+str(depth)+'\n')    
        output.append("___________________________________________________________________________________________________________________________________________________"+'_'*(longest_region*3)+"\n")
        output.append("\nSpecies"+' '*48+"Start"+' '*10+"End"+' '*12+"Width"+' '*10+"Depth"+' '*10+"Conserved Sites"+' '*(longest_region+10-15)+"Motif Neighborhood\n")
        output.append("___________________________________________________________________________________________________________________________________________________"+'_'*(longest_region*3)+"\n\n")
        for layer in layers:

            sites_for_layer = neighbourhoods[motif_id_set][layer][0]
            species = headers[layer-1][0:50]
            n1 = sites_for_layer[0]
            start = n1[0]
            #end1 = n1[1]

            n2 = sites_for_layer[-1]
            #start2 = n2[0]
            end = n2[1]
            
            kmer = layer_references[layer]

            ref_index = reference.index(n1[2])
            last_ref_index = reference.rindex(n2[2])
            ref_end_n = ref_end-(last_ref_index+len(n2[2]))

            if start>ref_index:
                ad1 = sequences[layer-1][start-ref_index-1:start-1]
                ad2 = sequences[layer-1][end:end+ref_end_n]
                n_start = start-ref_index
            else:
                ad1 = ' '*(ref_index-start+1)+sequences[layer-1][0:start-1]
                ad2 = sequences[layer-1][end:end+ref_end_n]
                n_start = 1

            if end+ref_end_n>len(sequences[layer-1]):
                n_end = len(sequences[layer-1])
            else:
                n_end=end+ref_end_n

            width = n_end-n_start+1
            #width = len(kmer)
            kmer_alone = ' '*len(kmer)
            
            #kmer_alone = ' '*ref_index+n[2]+' '*ref_end_n
            kmer_alone = '-'*len(kmer)
            kindex = -1
            kmer_upper = kmer.upper()
            for site in sites_for_layer:
                k = site[2].upper()
                kindex = kmer_upper.index(k,kindex+1)
                kend = kindex+len(k)

                kmer = kmer[0:kindex]+k+kmer[kend:] 
                kmer_alone = kmer_alone[0:kindex]+k+kmer_alone[kend:]
            kmer = ad1.lower()+kmer+ad2.lower()

            output.append(species+' '*(55-len(species)))
            output.append(str(n_start)+' '*(15-len(str(n_start))))
            output.append(str(n_end)+' '*(15-len(str(n_end))))
            output.append(str(width)+' '*(15-len(str(width))))
            output.append(str(layer)+' '*(15-len(str(layer))))
            output.append(kmer_alone+' '*(longest_region+10-len(kmer_alone)))
            output.append(kmer+"\n")
            #output.append(additional)



        output.append("___________________________________________________________________________________________________________________________________________________"+'_'*(longest_region*3)+"\n")
        output.append('\nIndividual Motif Sites in Neighborhood '+str(nh_id)+':\n')
        output.append('----------------------------------------------------------------\n\n')
        motif_id =1
        kmer_set.sort(key = lambda x: (-x[3],x[0]))
        kmer_checked = []
        for item in kmer_set:
            if item[2] in kmer_checked:
                continue
            kmer = item[2]
            kmer_checked.append(kmer)
            level = item[3]
            stats_text = '\nStatistical Significance Not Calculated\n'
            if stats_dict:
                stats = stats_dict[kmer]
                stats_text = '\nE(i)-value='+str(stats[1])+'    P(i)-value='+str(stats[2])+'    E(r)-value='+str(stats[3])+'    E(r)-value='+str(stats[4])+'\n'
            output.append("--------------------------------------------------------------------------------------------------------------------------------------------------"+'-'*(ref_end*3)+"\n")
            output.append('Motif '+str(nh_id)+'.'+str(motif_id)+'   Depth:'+str(level)+'\n')
            output.append(stats_text)
            output.append("--------------------------------------------------------------------------------------------------------------------------------------------------"+'-'*(ref_end*3)+"\n")
            output.append("Species"+' '*48+"Start"+' '*10+"End"+' '*12+"Width"+' '*10+"Depth"+' '*10+"          -Site-          \n")

            output.append("--------------------------------------------------------------------------------------------------------------------------------------------------"+'-'*(ref_end*3)+"\n")

            targetscan = ''
            eclip_blat = ''
            eclip_bed = ''
            site_for_all = dict_by_kmer[kmer]
            layers = [x for x in site_for_all.keys() if site_for_all[x] and x in neighbourhoods[motif_id_set]]
            layers.sort()
            for layer in layers:
                limit1 = neighbourhoods[motif_id_set][layer][1]
                limit2 = neighbourhoods[motif_id_set][layer][2]
                nodes = [x for x in site_for_all[layer] if x[0]>=limit1 and x[1]<=limit2]
                for n in nodes:
                    species = headers[layer-1][0:50]
                    start = n[0]
                    end = n[1]
                    if start>10:
                        ad1 = sequences[layer-1][start-11:start-1]
                        ad2 = sequences[layer-1][end:end+10]
                    else:
                        ad1 = sequences[layer-1][0:start-1]
                        ad2 = sequences[layer-1][end:end+10]
                    output.append(species+' '*(55-len(species)))
                    output.append(str(start)+' '*(15-len(str(start))))
                    output.append(str(end)+' '*(15-len(str(end))))
                    output.append(str(len(kmer))+' '*(15-len(str(len(kmer)))))
                    output.append(str(layer)+' '*(15-len(str(layer))))
                    output.append(ad1.lower()+' '*(10-len(ad1))+'-'+kmer+'-'+ad2.lower()+"\n")
                    if 'TargetScan' in Annotations[layer]:
                        anno=''
                        for match in Annotations[layer]['TargetScan'][(n[0],n[1],n[2])]:
                            anno+=match.split(':')[0]+','
                        if anno!='':
                            anno=species+':    '+anno+'\n'
                            targetscan+=anno

                    if 'eCLIP' in Annotations[layer]:
                        if 'BLAT' in Annotations[layer]['eCLIP']:
                            anno=''
                            for match in Annotations[layer]['eCLIP']['BLAT'][(n[0],n[1],n[2])]:
                                anno+=match[0][3]+','
                            if anno!='':
                                anno=species+':    '+anno+'\n'
                                eclip_blat+=anno

                        if 'BED' in Annotations[layer]['eCLIP']:
                            anno=''
                            for match in Annotations[layer]['eCLIP']['BED'][(n[0],n[1],n[2])]:
                                anno+=match[0][3]+','
                            if anno!='':
                                anno=species+':    '+anno+'\n'
                                eclip_bed+=anno

            if eclip_bed =='':
                eclip_bed = ' None\n'
            if eclip_blat =='':
                eclip_blat = ' None\n'
            if targetscan =='':
                targetscan = ' None\n'

            output.append("--------------------------------------------------------------------------------------------------------------------------------------------------"+'-'*(ref_end*3)+"\n")
            output.append('Motif '+str(nh_id)+'.'+str(motif_id)+' ('+kmer+') Annotations:\n')
            output.append('----------------------------------------------\n')
            output.append('TargetScan Matches:\n')
            output.append(targetscan)
            output.append('\n\n')
            output.append('eCLIP determined binding proteins (based on BLAT alignment):\n')
            output.append(eclip_blat)
            output.append('\n\n')
            output.append('eCLIP determined binding proteins (based on BED File):\n')
            output.append(eclip_bed)
            output.append("--------------------------------------------------------------------------------------------------------------------------------------------------"+'-'*(ref_end*3)+"\n\n\n\n")

            motif_id+=1
        nh_id+=1


    output.append("="*300+"\n")
    output.append("Motifs per sequence (Section 4 of 4)\n\n\n")
    layers = all_kmers.keys()
    layers.sort()
    output.append('-'*150+'\n')
    output.append('Start'+' '*10+'End'+' '*12+'Significance'+' '*20+'Motif'+' '*10+'\t\t\tTargetScan\t\t\teCLIP\n'+"-"*150+'\n')
    for layer in layers:
        output.append('\n'+headers[layer-1]+'\n\n')
        for kmer in all_kmers[layer]:
            stats_text = 'Not Calculated'
            if stats_dict:
                stats = stats_dict[kmer[0][2]]
                stats_text = str(stats[1])+','+str(stats[2])+','+str(stats[3])+','+str(stats[4])
            output.append(str(kmer[0][0])+' '*(15-len(str(kmer[0][0]))))
            output.append(str(kmer[0][1])+' '*(15-len(str(kmer[0][1]))))
            output.append(stats_text+' '*(32-len(stats_text)))
            output.append(str(kmer[0][2]))

            targetscan = ''
            if 'TargetScan' in Annotations[layer]:    
                for match in Annotations[layer]['TargetScan'][kmer[0]]:
                    targetscan+=match.split(':')[0]+','

            eCLIP = ''
            if 'eCLIP' in Annotations[layer]:
                if 'BLAT' in Annotations[layer]['eCLIP']:
                    matches = [] 
                    for match in Annotations[layer]['eCLIP']['BLAT'][kmer[0]]:
                        if not match[0][3] in matches: 
                            eCLIP+=match[0][3]+','
                            matches.append(match[0][3])


                if 'BED' in Annotations[layer]['eCLIP']:
                    matches = []
                    for match in Annotations[layer]['eCLIP']['BED'][kmer[0]]:
                        if not match[0][3] in matches: 
                            eCLIP+=match[0][3]+','
                            matches.append(match[0][3])

            output.append('\t\t\t'+targetscan+'\t\t\t'+eCLIP+'\n')
        output.append('_'*150+'\n')

    w= open(outdir+"/"+project_name+"/lncLOOM_Results.txt",'w',0)
    for l in output:
        w.write(l)
    w.close()


def print_trackRGB(kmers_in_query,kmer_graded_colours,chrm,strand,num_of_exons,size_of_exons,seq_exons,chromosome_coords,outdir,project_name,header,sequence):

    kmer_base_colours = {}
    #determine color of bases in each kmer
    nodes = []
    for kmer_set in kmers_in_query:
        kmer_id = kmer_set[0]
        nodes.append(kmer_id)
        full_start = kmer_set[0][0]
        overlaps = kmer_set[1]
        kmer_colours = {}
        kmer_names = {}
        for kmer in overlaps:
            relative_start = kmer[0]-full_start
            colour = kmer_graded_colours[kmer[2]][0]
            colour = colour.lstrip('rgb(')
            colour = colour.rstrip(')')
            for index in range(len(kmer[2])):
                relative_index = index+relative_start
                kmer_colours[relative_index]=colour
                kmer_names[relative_index]=kmer[2]

        start_colour = kmer_colours[0]
        #name = kmer_names[0]
        count_cons = 1
        kmer_in_steps = []
        reference = kmer_set[0][2]
        cursor = 0
        for i in range(1,len(reference)):
            next_colour = kmer_colours[i]
            if next_colour==start_colour:
                count_cons+=1
            else:
                kmer_in_steps.append((start_colour,count_cons,reference[cursor:i]))
                start_colour = next_colour
                cursor=i
                count_cons = 1
                #name = kmer_names[i]

        kmer_in_steps.append((start_colour,count_cons,reference[cursor:i+1]))
        kmer_base_colours[kmer_id] = kmer_in_steps

    #get chromosome coordinates for each kmer (reference to long combined kmer)



    #num_of_exons: holds the total number of exons in seq
    #size_of exons: holds the length of each exon
    #seq_exons: holds the start postion of each exon relative to the query sequence - index 0 based
    #chromosome_coords: holds the start position of each exon relative to the chromosome - index 0 based
    #nodes: holds an ordered list of tuples for each node (start,kmer) - index here is 1 based

    #1) Obtain positions of each node relative to chromosome co-ordinates (remember nodes are index 1 based)

    #Determine which nodes are located in which exon and if any nodes span two exons

    end_of_seq = chromosome_coords[-1]+size_of_exons[-1]
  
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
        node_start = nodes_starts[n] #node_start is 0 indexed: same as blat
        node_end = node_start+len(node[2]) #node_end is 1 indexed: same as blat
        
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


    #make a list of chrom start and end for each colour step in each kmer
    features = []
    for node in nodes:
        steps = kmer_base_colours[node]
        chromosome_pos = nodes_in_chromosome[node]
        carry_over = 0
        s = 0
        for exon in chromosome_pos:
            start = exon[1]
            end = exon[2]
            feat_end = start+steps[s][1]-carry_over
            while feat_end<=end and s<len(steps)-1:
                features.append([start,feat_end,steps[s][0],steps[s][2]])
                s+=1
                start=feat_end
                feat_end = start+steps[s][1]

            if feat_end<=end:
                features.append([start,feat_end,steps[s][0],steps[s][2]])
            else:
                features.append([start,end,steps[s][0],steps[s][2]])
                carry_over = feat_end-end

    output = 'browser position '+chrm+':'+str(features[0][0])+'-'+str(end_of_seq)+'\n'
    output+='track name="'+project_name+'" description="LncLOOM Conserved Motifs" itemRgb="On" visibility="1" maxItems="'+str(len(features))+'"\n'
    for feat in features:
        output+=chrm+'\t'+str(feat[0])+'\t'+str(feat[1])+'\t'+feat[3]+'\t0\t'+strand+'\t'+str(feat[0])+'\t'+str(feat[1])+'\t'+feat[2]+'\n'

    w = open(outdir+'/'+project_name+'/LOOM_'+header.strip('>')+'_RGB.bed','w')
    w.write(output)
    w.close()


def print_trackScore(kmers_in_query,kmer_graded_colours,chrm,strand,num_of_exons,size_of_exons,seq_exons,chromosome_coords,outdir,project_name,query_head,headers,depth_of_kmers,sequence):

    total = len(headers)
    scores = range(0,1000,1000/total)

 
    kmer_base_colours = {}
    #determine color of bases in each kmer
    nodes = []
    for kmer_set in kmers_in_query:
        overlaps = kmer_set[1]
        for kmer in overlaps:
            depth = depth_of_kmers[kmer[2]]
            score = scores[depth[0]-1]
            #name = headers[depth[0]-1]
            name = kmer[2]
            kmer_base_colours[kmer] = (score,name)
            nodes.append(kmer)

    #get chromosome coordinates for each kmer (reference to long combined kmer)



    #num_of_exons: holds the total number of exons in seq
    #size_of exons: holds the length of each exon
    #seq_exons: holds the start postion of each exon relative to the query sequence - index 0 based
    #chromosome_coords: holds the start position of each exon relative to the chromosome - index 0 based
    #nodes: holds an ordered list of tuples for each node (start,kmer) - index here is 1 based

    #1) Obtain positions of each node relative to chromosome co-ordinates (remember nodes are index 1 based)

    #Determine which nodes are located in which exon and if any nodes span two exons

    end_of_seq = chromosome_coords[-1]+size_of_exons[-1]
  
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
        node_start = nodes_starts[n] #node_start is 0 indexed: same as blat
        node_end = node_start+len(node[2]) #node_end is 1 indexed: same as blat
        
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


    #make a list of chrom start and end for each colour step in each kmer
    features = []
    for node in nodes:
        details = kmer_base_colours[node]
        chromosome_pos = nodes_in_chromosome[node]
        carry_over = 0
        s = 0
        for exon in chromosome_pos:
            start = exon[1]
            end = exon[2]
            features.append([exon[1],exon[2],details[0],details[1]])


    output = 'browser position '+chrm+':'+str(features[0][0])+'-'+str(end_of_seq)+'\n'
    output+='track name="'+project_name+'" description="LncLOOM Conserved Motifs" useScore="1" color="100,50,0" visibility="1" maxItems="'+str(len(features))+'"\n'


    for feat in features:
        output+=chrm+'\t'+str(feat[0])+'\t'+str(feat[1])+'\t'+feat[3]+'\t'+str(feat[2])+'\t'+strand+'\t'+str(feat[0])+'\t'+str(feat[1])+'\n'

    w = open(outdir+'/'+project_name+'/LOOM_'+query_head.strip('>')+'.bed','w')
    w.write(output)
    w.close()



















    
            
            





