#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 14 21:26:27 2017

@author: yi.yan
"""
def funProcessCountXL(sDir,lFiles,suffix1,suffix2,suffix3,numSeq):
    import xlsxwriter 
    import xlrd
    import math
    from scipy.stats import norm
    
    lSamples = [] 
    for file in lFiles :
        if file[-3:] == 'txt':
            tpIndex = file.index('BLAST') + 5 
            lSamples.append(file[0:tpIndex])
    
    lSamples = list(set(lSamples))
    for sample in lSamples :
        wFileName = sDir + sample + '_Summary.xlsx'
        workbook = xlsxwriter.Workbook(wFileName)
        worksheet = workbook.add_worksheet()
        
        sFile1 = sDir + sample + suffix1 
        sFile2 = sDir + sample + suffix2
        sFile3 = sDir + sample + suffix3
        bookNT = xlrd.open_workbook(sFile1)
        NTsheet = bookNT.sheet_by_index(0)
        bookARGOS = xlrd.open_workbook(sFile2)
        ARGOSsheet = bookARGOS.sheet_by_index(0)
        bookFDAARGOS = xlrd.open_workbook(sFile3)
        FDAARGOSsheet = bookFDAARGOS.sheet_by_index(0)
        
        NTSpecies = NTsheet.col_values(0) 
        ARGOSSpecies = ARGOSsheet.col_values(0)
        lSharedSpecies = list(set(NTSpecies[1:]).intersection(ARGOSSpecies[1:]))
        
        lTestStat = [] 
        for species in lSharedSpecies :
            numNT = NTsheet.cell(NTSpecies.index(species),2).value
            numARGOS = ARGOSsheet.cell(ARGOSSpecies.index(species),2).value
            p = (numNT + numARGOS)/(numSeq + numSeq)
            p1 = numNT/numSeq 
            p2 = numARGOS/numSeq 
            SE = math.sqrt(p * ( 1 - p) * ((1/numSeq) + (1/numSeq)) )
            z = (p1-p2)/SE 
            p_values = norm.sf(abs(z))*2
            lTestStat.append([species, numNT,numARGOS,p_values,numSeq])
            
        worksheet.write(0,0,'NT-Result')
        for i in range(0,NTsheet.nrows):
            for j in range(0,NTsheet.ncols):
                worksheet.write(i,j+1,NTsheet.cell(i,j).value)
        worksheet.write(0,7,'NT/ARGOS-Result')
        for i in range(0,ARGOSsheet.nrows):
            for j in range(0,ARGOSsheet.ncols):
                worksheet.write(i,j+8,ARGOSsheet.cell(i,j).value)
        worksheet.write(0,14,'Shared Results')
        
        newTestStat = sorted(lTestStat,key=lambda x: x[3])
        TestHeader = ['Shared Species Name','num in NT','num in ARGOS','p-val','num-Seq']
        for i in range (0,5):
            worksheet.write(0,i+15,TestHeader[i])
        for i in range (0,5):
            for j in range(0,len(newTestStat)):
                worksheet.write(j+1,i+15,newTestStat[j][i])
        worksheet.write(0,20,'ARGOS Only')
        for i in range(0,FDAARGOSsheet.nrows):
            for j in range(0,FDAARGOSsheet.ncols):
                worksheet.write(i,j+21,FDAARGOSsheet.cell(i,j).value)       
        workbook.close()