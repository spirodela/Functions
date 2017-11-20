#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  3 18:42:25 2017

@author: yi.yan
"""
def funANICalc(tblComplete,lSeqRecord,dbName):
    import subprocess 
    from Bio import SeqIO

    ANIOut = []
    for temp in tblComplete :
        seqName = temp[5]
        ANIRow = []
        for tempSeq in lSeqRecord :
            if seqName == tempSeq.id :
                print('Calculate Specific ANI')
                SeqIO.write(tempSeq,"ANIInput/Input_1.fasta","fasta")
        process = subprocess.Popen("export BLASTDB=/Users/yi.yan/Documents/db/:$BLASTDB"\
                                   +"&&/usr/local/ncbi/blast/bin/"\
                                   +"blastdbcmd "\
                                   +" -db "\
                                   +dbName\
                                   +" -entry "\
                                   +temp[-5]\
                                   +" -out "\
                                   +" ANIInput/Input_2.fasta",\
                                   shell=True,\
                                   stdout = subprocess.PIPE,\
                                   stderr = subprocess.PIPE)
        proc_out, proc_err = process.communicate()
                
        process = subprocess.Popen("export PATH=/usr/local/ncbi/blast/bin/:$PATH"+\
                                   "&&average_nucleotide_identity.py" +\
                                   " -i ANIInput"+\
                                   " -o ANIOutput"+\
                                   " -f "+\
                                   " -m ANIb",
                                   shell=True,\
                                   stdout = subprocess.PIPE,\
                                   stderr = subprocess.PIPE)
        proc_out, proc_err = process.communicate()
        
        file = open("ANIOutput/ANIb_alignment_coverage.tab","r")
        Lines = file.readlines()
        ANIRow.append(Lines[1].split()[2])
        ANIRow.append(Lines[2].split()[1])
        file = open("ANIOutput/ANIb_alignment_lengths.tab","r")
        Lines = file.readlines()
        ANIRow.append(Lines[1].split()[2])
        ANIRow.append(Lines[2].split()[1])
        file = open("ANIOutput/ANIb_hadamard.tab","r")
        Lines = file.readlines()
        ANIRow.append(Lines[1].split()[2])
        ANIRow.append(Lines[2].split()[1])
        file = open("ANIOutput/ANIb_percentage_identity.tab","r")
        Lines = file.readlines()
        ANIRow.append(Lines[1].split()[2])
        ANIRow.append(Lines[2].split()[1])
        file = open("ANIOutput/ANIb_similarity_errors.tab","r")
        Lines = file.readlines()
        ANIRow.append(Lines[1].split()[2])
        ANIRow.append(Lines[2].split()[1])
        ANIOut.append(ANIRow.copy())
                 
    
    
    count = 0 
    for temp in tblComplete:
        tblComplete[count][6] = float(tblComplete[count][6])
        count = count+1 
   
    FinalTbl = []
    count = 0 
    for temp in tblComplete:
        tempRow = tuple(temp + ANIOut[count])
        FinalTbl.append(tempRow)
        count = count + 1
    return FinalTbl 