#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  3 18:42:22 2017

@author: yi.yan
"""

def funReadBlast(sOutFileName,all_species,sARGOSID,nNumContig,nTotalAssemblySize,\
                 nN50,nLargestContig,lGC,lSize):
    import asciitable 
    tblBlast = asciitable.read(sOutFileName, Reader = asciitable.NoHeader, delimiter ='\t')
    
    """Get Name of Contigs""" 
    lContigName = list(set(tblBlast['col1']))
    
    tblComplete = [];
    for j in range(0,len(lContigName)):
        print('Processing BLAST Result for Contig')
        iContigIndices = [i for i,x in enumerate(list(tblBlast['col1'])) if x == lContigName[j]]
        lAccesions = list(set(tblBlast[iContigIndices]['col4']))
        
        temp = all_species[j]
    
        lTempRow = [sARGOSID] #0
        lTempRow.append(str(nNumContig)) #1
        lTempRow.append(str(nTotalAssemblySize)) #2
        lTempRow.append(str(nN50)) #3
        lTempRow.append(str(nLargestContig)) #4
        lTempRow.append(lContigName[j]) #5
        lTempRow.append(str(tblBlast[iContigIndices[0]][1])) #6
        lTempRow.append(str(lGC[lSize.index(tblBlast[iContigIndices[0]][1])])) #7
        lTempRow.append(temp[0:-3]) #8
        lTempRow.append('Hit Name') #9
        lTempRow.append('Hit ACC') #10
        lTempRow.append('Hit Score') #11
        lTempRow.append('Hit PIdent') #12
        lTempRow.append('Hit PCOV') #13
        lTempRow.append('Scientific_Names') #14

        for k in range(0,len(lAccesions)):
            iAccesionIndices = [i for i,x in enumerate(list(tblBlast['col4'])) if x == lAccesions[k]]
            iAccesionIndices = list(set(iAccesionIndices).intersection(iContigIndices))
            name = tblBlast[iAccesionIndices[0]]['col5']
            name = name.replace(",","_")
            lTempRow[9] = name
            lTempRow[10] = lAccesions[k]
            lTempRow[-2] = tblBlast[iAccesionIndices[0]]['col9']
            lTempRow[-1] = tblBlast[iAccesionIndices[0]]['col3']
            nPident = 0
            nAsize = 0
            nScore = 0 
            for l in range(0,len(iAccesionIndices)):
                nPident = nPident + tblBlast[iAccesionIndices[l]]['col8']*tblBlast[iAccesionIndices[l]]['col6']
                nAsize = nAsize + tblBlast[iAccesionIndices[l]]['col6']
                nScore = nScore + tblBlast[iAccesionIndices[l]]['col7']
                
            nPident = nPident/nAsize
            lTempRow[11] = nScore
            lTempRow[12] = nPident
            tblComplete.append(lTempRow.copy())
    return tblComplete 