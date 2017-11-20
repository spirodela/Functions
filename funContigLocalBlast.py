#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 25 15:02:59 2017

@author: yi.yan
"""

# -*- coding: utf-8 -*-
"""
Created on Tue Sep 12 9:14:36 2017
This is a function that automatically runs blast nt on a input sequence 
input: sFileName, Name of the input fasta sequence file 
output: Score, Identity % and Coverage Percent of the top 5 blast hits
@author: yihen
"""

def funLocalBlast(sFastaFileName,sGBKFileName,dbName):
    """Import packages used """
    from Bio.Blast.Applications import NcbiblastnCommandline
    from Bio import SeqIO
    from Bio.SeqUtils import GC
    import subprocess
    import xlsxwriter
    from funReadBlast import funReadBlast
    from funANICalc import funANICalc 
    from funBlastANI2XLS import funBlastANI2XLS
    
    #sFastaFileName = '/Users/yi.yan/Documents/FDA-ARGOS/Batch6/PFDA1_Batch6_Misidentified_Contaminated/CNH_804.fasta'
    #sGBKFileName = "TestFolderGenBank/AMERTCC_44.annotation.20161209.gbk"
    #sGBKFileName = 'N/A'
    #dbName = "ref_prok_rep_genomes"
    
    columnTitleRow = ["FDAARGOS_ID",#0
                 "Num_Contig",#1 
                 "Assembly_Size",#2
                 "N_50", #3
                 "Largest_Contig_Size", #4
                 "Contig_ID", #5
                 "Contig_Length", #6
                 "Contig_GC", #7
                 "Proposed Organism", #8
                 "Blast_Hit", #9
                 "ACCESSION", #10
                 "Score", #11
                 "Percent_Query_Identity",#12
                 "Percent_Query_Coverage", #13
                 "Scientific_Name", #14
                 "Query_ANI_Coverage", #15
                 "Subject_ANI_Coverage", #16
                 "Query_ANI_Length", #17
                 "Subject_ANI_Length", #18
                 "Query_ANI_HD", #19
                 "Subject_ANI_HD", #20
                 "Query_ANI_Identity", #21
                 "Subject_ANI_Identity", #22
                 "Query_ANI_SE", #23
                 "Subject_ANI_SE"]
            
    sFileName = sFastaFileName[0:-6]+'.xlsx' 
    
    lARGOSID = sFastaFileName.split("/")
    sARGOSID = lARGOSID[-1][0:-6]
    
    """Import Fasta sequence from assembly file"""
    
    lSeqRecord = []
    for seq_record in SeqIO.parse(sFastaFileName, "fasta"):
        lSeqRecord.append(seq_record)
    
    """Import Annotation"""
    all_species = []

    if sGBKFileName == "N/A":
        for seq_record in lSeqRecord :
            all_species.append('N/A-N/A')
    else:
        f = open(sGBKFileName,'r',errors = 'ignore')
        for line in f :
            if "ORGANISM" in line:
                print(line)
                sSpecie = line
                all_species.append(sSpecie)
        f.close
    """Calculate Contig Statistics"""
    lSize = []
    lGC = []
    
    for seq_record in lSeqRecord: 
        lSize.append(len(seq_record.seq))
        lGC.append(GC(seq_record.seq))
        
    nTotalAssemblySize = sum(lSize)
    nNumContig = len(lSize)
    nLargestContig = max(lSize)
    
    #nTotalGC = np.multiply(lGC,lSize)
    #nTotalPercentGC = sum(nTotalGC)/nTotalAssemblySize
    
    nThreshold = 0.5*nTotalAssemblySize 
    lTempSize = sorted(lSize,reverse=True)
    
    nSize = 0
    count = 0
    
    while nSize <= nThreshold:
        nSize = nSize + lTempSize[count]
        out = count
        count = count + 1
        
    nN50 = lTempSize[out]    
  
#Run Blast 

    sOutFileName = sARGOSID + ".txt"
    
    
    blastn_cline = NcbiblastnCommandline(task = "megablast", \
                                         query = sFastaFileName, \
                                         db = dbName,\
                                         evalue = 0.001, \
                                         max_target_seqs = 5, \
                                         outfmt = "\"6 " +\
                                         "qseqid "+\
                                         "qlen "+\
                                         "sscinames "+\
                                         "sacc "+\
                                         "stitle "+\
                                         "length "+\
                                         "score "+\
                                         "pident "+\
                                         "qcovs\"",
                                         out = sOutFileName)
    
    print('run Refseq Blast')
    process = subprocess.Popen("export BLASTDB=/Users/yi.yan/Documents/db/:$BLASTDB"\
                               +"&&/usr/local/ncbi/blast/bin/"\
                               +str(blastn_cline),\
                               shell=True,\
                               stdout = subprocess.PIPE,\
                               stderr = subprocess.PIPE)
    proc_out, proc_err = process.communicate()
    

    tblComplete = funReadBlast(sOutFileName,all_species,sARGOSID,nNumContig,nTotalAssemblySize,\
                 nN50,nLargestContig,lGC,lSize)        
#Run ANI 
    
    print('Calculating ANI')
    FinalTbl = funANICalc(tblComplete,lSeqRecord,dbName)        
    s = sorted(FinalTbl, key=lambda x:(x[6],x[11]),reverse=True)
    workbook = xlsxwriter.Workbook(sFileName)
    funBlastANI2XLS(workbook,s,dbName,columnTitleRow)
    #BLAST NT 
    
    sOutFileName = sARGOSID + ".txt"
    
    
    blastn_cline = NcbiblastnCommandline(task = "megablast", \
                                         query = sFastaFileName, \
                                         db = "nt",\
                                         evalue = 0.001, \
                                         max_target_seqs = 5, \
                                         outfmt = "\"6 " +\
                                         "qseqid "+\
                                         "qlen "+\
                                         "sscinames "+\
                                         "sacc "+\
                                         "stitle "+\
                                         "length "+\
                                         "score "+\
                                         "pident "+\
                                         "qcovs\"",
                                         out = sOutFileName)
    
    print('Run BLAST NT')
    process = subprocess.Popen("export BLASTDB=/Users/yi.yan/Documents/db/:$BLASTDB"\
                               +"&&/usr/local/ncbi/blast/bin/"\
                               +str(blastn_cline),\
                               shell=True,\
                               stdout = subprocess.PIPE,\
                               stderr = subprocess.PIPE)
    proc_out, proc_err = process.communicate()
    
    
    tblComplete = funReadBlast(sOutFileName,all_species,sARGOSID,nNumContig,nTotalAssemblySize,\
                 nN50,nLargestContig,lGC,lSize)        
#Run ANI 
    
    print('Calculating ANI')
    FinalTbl = funANICalc(tblComplete,lSeqRecord,"nt")        
    s = sorted(FinalTbl, key=lambda x:(x[6],x[11]),reverse=True)
    funBlastANI2XLS(workbook,s,'NT',columnTitleRow)

        
    workbook.close()