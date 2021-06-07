# -*- coding: utf-8 -*-
"""
Created on Wed Dec 30 22:36:44 2020

@author: Mostafa
"""


import re
from Bio import SeqIO
import csv
import os


#patterns to search for 

seq1 = "CCCTAACCCTAA"
seq2 = "CCTAACCCTAAC"
seq3 = "CTAACCCTAACC"
seq4 = "TAACCCTAACCC"
seq5 = "AACCCTAACCCT"
seq6 = "ACCCTAACCCTA"

seq7 = "TTAGGGTTAGGG"
seq8 = "TAGGGTTAGGGT"
seq9 = "AGGGTTAGGGTT"
seq10 = "GGGTTAGGGTTA"
seq11 = "GGTTAGGGTTAG"
seq12 = "GTTAGGGTTAGG"

###  import the final location of each Terminal MoTeRs + 50 location 
ends = []
with open('~/Terminal_MoTeR_Positions_V4.csv', encoding='utf-8-sig') as csv_file:
    csv_reader = csv.reader(csv_file, delimiter='\t')
    for row in csv_reader:
        ends.append(row)
#print(ends)

### genomes folder
os.chdir("/Users/mostafa/Google Drive/Final_MinION_Genomes/Final_genomes/Final_genomes/")


# find the relevant files in this directory
files = [f for f in os.listdir('.') if f.endswith('.fasta')]



counts = []
found_dim=[]

for file in files: 
    print(file)
    name = re.sub('\.fasta$', '', file)
    #print(name)
    #file = "/Users/mostafa/Google Drive/Final_MinION_Genomes/CD156_Final.fasta"

    
    pattern = [seq1,seq2,seq3,seq4,seq5,seq6,seq7,seq8,seq9,seq10,seq11,seq12] #
    regex = re.compile(r'({})+'.format('|'.join(pattern)))
    #print(regex)
    parse=[i for i in SeqIO.parse(file, "fasta")]
    #print(parse)
    
    

    for seq_record in parse:
        
        Tetramer = 0
        Trimer =0
        Dimer = 0
        lenghts = 0
        print(seq_record.id)
        
        results = [(match.start(), match.end()) for match in regex.finditer(str(seq_record.seq))]

        
        #check if the patterns are bigger than 2 motifs (12bp)
        ls = []        
        for ind in results:       
            motif_ln = ind[1]-ind[0]
            if motif_ln > 11:
                #print(motif_ln)
                ls.append(ind)
                
        #print((ls))
        
        # checking if the position are farter than 10bp from each other        
        for index,tupl in enumerate(ls):
            if len(ls)<=1: #check length>1 tuple
                pass
            elif (index!=0):
                distance=tupl[0]-ls[index-1][1]# distance of consecutive elements
                if distance<50:
                    ls[index]=(ls[index-1][0],ls[index][1]) #change current tuple
                    ls[index-1]=0 #make previous tuple zero
        lss=[el for el in ls if el!=0] # remove zeros  
                    
                
        #print((lss))
                
            
        for ind in lss:
            #print(ind)

            index = ind[0]        
            motif_ln = ind[1]-index
            #print(motif_ln)
               
            if motif_ln > 23:
                found_dim.append((file, "Tetramer", index+1, ind[1],seq_record.id, seq_record.seq[index:ind[1]], seq_record.seq[index-50:index], seq_record.seq[ind[1]:ind[1]+50])) #, seq_record.seq[index-50:index], seq_record.seq[ind[1]:ind[1]+50]
                Tetramer +=1
                #print(index+1, ind[1],seq_record.id, seq_record.seq[index-50:index])
            elif motif_ln > 17 and motif_ln <24:
                found_dim.append((file,"Trimer", index+1, ind[1],seq_record.id, seq_record.seq[index:ind[1]], seq_record.seq[index-50:index], seq_record.seq[ind[1]:ind[1]+50]))
                Trimer +=1
                #print(Trimer)
            elif motif_ln > 11 and motif_ln <18:
                found_dim.append((file,"Dimer", index+1, ind[1] ,seq_record.id, seq_record.seq[index:ind[1]], seq_record.seq[index-50:index], seq_record.seq[ind[1]:ind[1]+50]))
                Dimer +=1  
                              

#### removing dims that sitting in the tel repeats

final = []
for i in ends:
    for rec in found_dim:
        if i[0] == rec[0] and i[1] == rec[4]:
            #print(i, "######################################")
            #print(rec)
            if int(rec[2]) > int(i[2]) and int(rec[3]) < int(i[3]):
                final.append(rec)
                
                #print(i, "######################################")
                #print(rec)            
        
#print(len(final))


##########summarizing the data      
#### counting the di-, tri- and tetra- mers
import pandas as pd 

counting = []
for i in final:
    counting.append((i[0], i[1]))
    
#print(counting)
df = pd.DataFrame(counting, columns=["Genome","Mer"])
#print(df)

count_series = df.groupby(["Genome","Mer"]).size()
new_df = count_series.to_frame(name = 'size').reset_index()
#print(new_df)   

new_df.to_csv('~/Summary_Internal_Tels_counts5_50bpGap_4.csv', index=False)
        
