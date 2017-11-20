# -*- coding: utf-8 -*-
"""
Created on Tue Sep 12 9:14:36 2017
This is a function that automatically runs blast nt on a input sequence 
input: sFileName, Name of the input fasta sequence file 
output: Score, Identity % and Coverage Percent of the top 5 blast hits
@author: yihen
"""

def funContigBlast(sFastaFileName,sGBKFileName,csv1,csv2):
    """Import packages used """
    from Bio.Blast import NCBIWWW
    from Bio.Blast import NCBIXML 
    from Bio import SeqIO
    from Bio.SeqUtils import GC
    import numpy as np
    from ete3 import NCBITaxa
    sFastaFileName = "TestFolderFasta/AMERTCC_31.fasta"
    sGBKFileName = "TestFolderGenBank/AMERTCC_31.annotation.20161209.gbk"
    
    lARGOSID = sFastaFileName.split("//")
    sARGOSID = lARGOSID[-1][0:-6]
    
    """Import Fasta sequence from assembly file"""
    
    lSeqRecord = []
    for seq_record in SeqIO.parse(sFastaFileName, "fasta"):
        lSeqRecord.append(seq_record)
    
    """Import Annotation"""
    all_species = []
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
    
    nTotalGC = np.multiply(lGC,lSize)
    nTotalPercentGC = sum(nTotalGC)/nTotalAssemblySize
    
    nThreshold = 0.5*nTotalAssemblySize 
    lTempSize = sorted(lSize,reverse=True)
    
    nSize = 0
    count = 0
    
    while nSize <= nThreshold:
        nSize = nSize + lTempSize[count]
        out = count
        count = count + 1
        
    nN50 = lTempSize[out]
    
    count = 0;
    
    """Obtain NCBI TaxID of the Genus """
    ncbi = NCBITaxa()
    
    for seq_record in lSeqRecord:
        
        """Check Organism Annotation """
    
        temp = all_species[count].split()
        query = temp[1]
        name2taxid = ncbi.get_name_translator([query])
        
        """Modify Genus ID to be compatible with blast input"""
    
        sOrganism = "txid" + str(name2taxid[query][0]) + " [ORGN]"
        
    
        result_handle = NCBIWWW.qblast(program= "blastn", \
                                       megablast = 'TRUE', \
                                       database= "ref_prok_rep_genomes", \
                                       sequence= seq_record.seq, \
                                       format_type="XML", \
                                       entrez_query= sOrganism)
    
        count = count + 1
        
        blast_record = NCBIXML.read(result_handle)
        
        count2 = 0
        
        for blast_description in blast_record.descriptions:
            title = blast_description.title
            
            nQueryCoverage = 0
            nIdentity = 0
            
            for hsp in blast_record.alignments[count2].hsps:
                if hsp.expect <= 0.05: 
                    nIdentity = nIdentity + hsp.identities 
                    nQueryCoverage = nQueryCoverage + hsp.align_length - hsp.gaps
            
            nIdentity = nIdentity/nQueryCoverage 
            nQueryCoverage = nQueryCoverage/lSize[count-1]
            
            row = sARGOSID \
            + "," + str(nNumContig) \
            + "," + str(nTotalAssemblySize) \
            + "," + str(nN50) \
            + "," + str(nLargestContig) \
            + "," + str(nTotalPercentGC) \
            + "," + seq_record.id\
            + "," + str(lSize[count-1]) \
            + "," + title.replace(","," ") \
            + "," + str(blast_description.score) \
            + "," + str(blast_description.e) \
            + "," + str(blast_description.num_alignments) \
            + "," + str(nIdentity) \
            + "," + str(nQueryCoverage) + '\n'
            
            csv1.write(row)
            
            if count2 == 0:
                csv2.write(row)
                
            count2 = count2 + 1
