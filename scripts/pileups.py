#!/usr/bin/env python

import sys
import pandas as pd
import pysam 


BED       = pd.read_csv("input/GOI.bed", sep="\t", names=["chrom", "start", "stop", "anno"])
RES_PATH  = sys.argv[1]
BAM       = pysam.AlignmentFile(sys.argv[2], "rb")
SAMPLE_ID = sys.argv[3] 



BED["anno"] = BED["anno"].str.split("|").str[1]

# numbering the different targets of a gene
checked=[]
list_of_lengths=[]

for gene in BED["anno"].tolist():
    if gene in checked:
        continue
    checked.append(gene)
    for i in range(1, BED["anno"].tolist().count(gene)+1):
        list_of_lengths.append(i)
        
BED["gene"] = BED["anno"]+"_"+pd.Series(list_of_lengths).astype(str)

df = pd.DataFrame()
for _, (chrom,start,stop,_, gene) in BED.iterrows(): # for each target, pileup of each position
    curr_df = pd.DataFrame()
    
    curr_df["pos"] = [i for i in range(start, stop+1)]
    curr_df["chrom"]=chrom
    curr_df["anno"]= gene
    curr_df.set_index(["pos"], inplace=True)
    curr_df[["A", "C", "G", "T", "N", "DEL", "COV", "INS", "REF"]]= 0, 0, 0, 0, 0, 0, 0, 0, "N"
    
    print(chrom, start, stop, gene)
    
    for read in BAM.fetch(str(chrom), start-100, stop+100):

        genomic_pos = read.reference_start
        read_pos=0
        
        for operation, length in read.cigartuples: # iterate through cigartuples

            if operation == 0: # match
                if len(set([i for i in range(genomic_pos, genomic_pos+length+1)]).intersection(set(curr_df.index))):
                    for indv_read, indv_genome in enumerate(range(genomic_pos, genomic_pos+length+1)):
                        base = read.query_sequence
                        if indv_genome in curr_df.index:
                            try:
                                curr_df.at[indv_genome, base[read_pos+indv_read]] += 1
                                curr_df.at[indv_genome, "REF"] = base[read_pos+indv_read]
                            except IndexError:
                                pass
                genomic_pos += length
                read_pos += length

            elif operation == 1: # insertion
                if genomic_pos in curr_df.index:
                    for num_bases, _ in enumerate(range(read_pos, read_pos+length+1)):
                        if genomic_pos+num_bases in curr_df.index:
                            if genomic_pos+num_bases in curr_df.index:
                                curr_df.at[genomic_pos+num_bases, "INS"] +=1
                read_pos+=length
            elif operation== 2: # deletion
                if len(set([i for i in range(genomic_pos, genomic_pos+length+1)]).intersection(set(curr_df.index))):
                    for _, indv_genome in enumerate(range(genomic_pos, genomic_pos+length+1)):
                        if indv_genome in curr_df.index:
                            curr_df.at[indv_genome, "DEL"] += 1
                genomic_pos+=length 
            elif operation==3: # refskip
                genomic_pos += length
            elif operation == 4: # softclip
                read_pos += length
            else: # hardclips and trash
                pass
            

    curr_df["COV"] = curr_df["A"]+curr_df["C"]+curr_df["G"]+curr_df["T"]+curr_df["N"]
    
    df = pd.concat([curr_df, df])

df.reset_index(inplace=True)

df_overview = pd.DataFrame()
for anno in set(df["anno"]):
      
    helper = df[df["anno"]==anno]
    mean_cov = helper["COV"].mean() 
    median_cov = helper["COV"].median()

    bases = len(helper)
    covered = len(helper[helper["COV"]>0])
    try:
        x= bases/covered
    except ZeroDivisionError:
        x=0
        
    curr = {"anno":anno, "vertical_mean": helper["COV"].mean(), "vertical_median": helper["COV"].median(), "bases":bases, "covered":covered, "horizontal_coverage":x}
    df_overview = pd.concat([df_overview, pd.DataFrame.from_dict(curr.items()).set_index(0).T])
    
writer = pd.ExcelWriter(f"results/{SAMPLE_ID}.xlsx", engine="xlsxwriter")

# used bed file
BED.to_excel(writer, index=False, sheet_name="bed_file")

# coverage overview
df_overview.sort_values("anno", inplace=True)
df_overview.set_index("anno", inplace=True)
df_overview.to_excel(writer,  sheet_name="coverage")

# noticeable positions
noticed=pd.DataFrame()

# make makeshift customized vcf file
vcf=pd.DataFrame()

# allele frequency 0.1
for index, (pos, chrom, anno, a, c, g, t, n, dele, cov, ins, ref) in df.iterrows():
    
    dict_bases={"A":a, "C":c, "G":g, "T":t, "DEL":dele, "INS":ins}

    if dict_bases[ref] == cov or ref == "N":
        continue
    
    dict_data={}
    
    if a>0:
        dict_data["A"] = round(a/cov,2)
    
    if c>0:
        dict_data["C"] = round(c/cov,2)
   
    if g>0:
        dict_data["G"] = round(g/cov,2)
    
    if t>0:
        dict_data["T"] = round(t/cov,2)        
   
    # remove reference 
    if dict_data[ref] >= 0.9:
        continue
    
    del dict_data[ref]

    to_remove=[]
    for key, val in dict_data.items():
        if val < 0.1:
            to_remove.append(key)
    
    for key in to_remove:
        del dict_data[key]

    dict_data = {k: v for k, v in reversed(sorted(dict_data.items(), key=lambda item: item[1]))}
    
    # vcf information
    ALT=",".join(dict_data.keys())
    DP=cov
    AF=",".join(list(map(str, dict_data.values())))
    
    GT="0/1" # to keep everything

    vcf_dict = {"#CHROM":chrom, "POS":pos, "REF":ref, "ALT":ALT, "INFO":f"DP={DP};AF={AF}", "FORMAT":f"DP:AF:GT", "SAMPLE":f"{DP}:{AF}:{GT}", "ID":".", "QUAL":100, "FILTER":""}
    
    vcf = pd.concat([vcf, pd.DataFrame.from_dict(vcf_dict.items()).set_index(0).T])

vcf = vcf[["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "SAMPLE"]]
vcf.sort_values(["#CHROM", "POS"], inplace=True)
vcf.to_csv(f"results/{SAMPLE_ID}.vcf", sep="\t", index=False)
vcf.to_excel(writer, index=False, sheet_name="snv")


# deletions and insertions
df_indel = df[(df["INS"] >= 1) | (df["DEL"] >= 1)]
df_indel["freq_del"] = df_indel["DEL"]/df_indel["COV"]
df_indel["freq_ins"] = df_indel["INS"]/df_indel["COV"]
df_indel = df_indel[(df_indel["freq_del"] >=0.1) | (df_indel["freq_ins"] >= 0.1)]

df_indel = df_indel[["anno", "pos", "chrom", "REF", "COV", "DEL","INS"]]
df_indel.to_excel(writer, index=False, sheet_name="indel")

# detalis with each position
df = df[["anno", "chrom", "pos", "COV", "A", "C", "G", "T", "N", "DEL", "INS"]]
df.to_excel(writer, index=False, sheet_name="details")

writer.close()