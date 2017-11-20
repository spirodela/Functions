#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  3 20:41:44 2017

@author: yi.yan
"""

def funGenusLocalBlast(sFastaFileName,sGBKFileName,dbName):
    """Import packages used """
    from Bio.Blast.Applications import NcbiblastnCommandline
    from Bio import SeqIO
    from Bio.SeqUtils import GC
    import subprocess
    import xlsxwriter
    from funReadBlast import funReadBlast
    from funANICalc import funANICalc 
    from collections import Counter
    from ete3 import NCBITaxa

    #sFastaFileName = "TestFolderFasta/AMERTCC_31.fasta"
    #sGBKFileName = "TestFolderGenBank/AMERTCC_31.annotation.20161209.gbk"
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
            
    sFileName = sFastaFileName[0:-6]+'_Genus_Blast.xlsx' 
    workbook = xlsxwriter.Workbook(sFileName)
    
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
    ncbi = NCBITaxa()

    Genus = all_species[0].split()
    Genus = Genus[1]
    name2taxid = ncbi.get_name_translator([Genus])
    sOrganism = "\"txid" + str(name2taxid[Genus][0]) + " [ORGN]\""

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
                                         entrez_query = sOrganism,
                                         remote = 1,
                                         out = sOutFileName)
    
    
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
    

    FinalTbl = funANICalc(tblComplete,lSeqRecord,'nt')  
    
    s = sorted(FinalTbl, key=lambda x:(x[6],x[11]),reverse=True)
    
    ContigNames = [i[5] for i in s]
    lContigName = list(set(ContigNames))
    SummaryTbl = []
    for Name in lContigName:
        SummaryTbl.append(s[ContigNames.index(Name)])
    
    SummaryTbl = sorted(SummaryTbl, key=lambda x:(x[6],x[11]),reverse=True)

    worksheet = workbook.add_worksheet('NT_Genus')
    
    for i in range(0,len(s[0])):
        worksheet.write(0,i,columnTitleRow[i])
        
    for i in range(0,len(s)):
        for j in range(0,len(s[0])):
            worksheet.write(i+1,j,s[i][j])
    #Specie Distribution         
    nRowStart = i + 5 
    lSciNames = [i[14] for i in s]
    lGenus = []
    lSpecies = []
    for i in lSciNames :
        temp = i.replace('[','')
        temp = temp.replace(']','')
        lGenus.append(temp.split()[0])
        lSpecies.append(temp.split()[0] + ' ' + temp.split()[1])
        
    c= Counter(lGenus)
    
    nItem = len(lGenus)
    Genus = list(c.keys())
    
    tempHeader = ['Genus','Count','Percentage']
    
    for i in range(0,3):    
        worksheet.write(nRowStart,i,tempHeader[i])
    
    for i in range(0,len(Genus)):
        worksheet.write(nRowStart+i,0,Genus[i])
        worksheet.write(nRowStart+i,1,c.get(Genus[i]))
        worksheet.write(nRowStart+i,2,c.get(Genus[i])/nItem)
    
    nRowStart = nRowStart + i + 5    
    
    c= Counter(lSpecies)

    nItem = len(lSpecies)
    Species = list(c.keys())
    
    tempHeader = ['Genus','Count','Percentage']
    
    for i in range(0,3):    
        worksheet.write(nRowStart,i,tempHeader[i])
    
    for i in range(0,len(Species)):
        worksheet.write(nRowStart+i,0,Species[i])
        worksheet.write(nRowStart+i,1,c.get(Species[i]))
        worksheet.write(nRowStart+i,2,c.get(Species[i])/nItem)
  #Summary       
    worksheet = workbook.add_worksheet('NT_Summary_Genus')
    
    for i in range(0,len(s[0])):
        worksheet.write(0,i,columnTitleRow[i])
        
    for i in range(0,len(SummaryTbl)):
        for j in range(0,len(s[0])):
            worksheet.write(i+1,j,SummaryTbl[i][j])
    
#Specie Distribution         
    nRowStart = i + 5 
    lSciNames = [i[14] for i in SummaryTbl]
    lGenus = []
    lSpecies = []
    for i in lSciNames :
        temp = i.replace('[','')
        temp = temp.replace(']','')
        lGenus.append(temp.split()[0])
        lSpecies.append(temp.split()[0] + ' ' + temp.split()[1])
        
    c= Counter(lGenus)
    
    nItem = len(lGenus)
    Genus = list(c.keys())
    
    tempHeader = ['Genus','Count','Percentage']
    
    for i in range(0,3):    
        worksheet.write(nRowStart,i,tempHeader[i])
    
    for i in range(0,len(Genus)):
        worksheet.write(nRowStart+i,0,Genus[i])
        worksheet.write(nRowStart+i,1,c.get(Genus[i]))
        worksheet.write(nRowStart+i,2,c.get(Genus[i])/nItem)
    
    nRowStart = nRowStart + i + 5    
    
    c= Counter(lSpecies)

    nItem = len(lSpecies)
    Species = list(c.keys())
    
    
    for i in range(0,len(Species)):
        worksheet.write(nRowStart+i,0,Species[i])
        worksheet.write(nRowStart+i,1,c.get(Species[i]))
        worksheet.write(nRowStart+i,2,c.get(Species[i])/nItem)   
        
    workbook.close()