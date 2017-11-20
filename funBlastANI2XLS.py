#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 18 09:18:19 2017

@author: yi.yan
"""
def funBlastANI2XLS(workbook,s,dbName,columnTitleRow):
    import xlsxwriter
    from collections import Counter

    ContigNames = [i[5] for i in s]
    lContigName = list(set(ContigNames))
    SummaryTbl = []
    for Name in lContigName:
        SummaryTbl.append(s[ContigNames.index(Name)])
    
    SummaryTbl = sorted(SummaryTbl, key=lambda x:(x[6],x[11]),reverse=True)

    worksheet = workbook.add_worksheet(dbName)
    
    for i in range(0,len(s[0])):
        worksheet.write(0,i,columnTitleRow[i])
        
    for i in range(0,len(s)):
        for j in range(0,len(s[0])):
            worksheet.write(i+1,j,s[i][j])
    #Specie Distribution         
    nRowStart = i + 5 
    lSciNames = [i[14] for i in s]
    lContigLength = [i[6] for i in s]

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
    
    print('Calculating Specie Distribution')
    tempHeader = ['Genus','Count','Percentage','Total Contig Length']
    
    for i in range(0,4):    
        worksheet.write(nRowStart,i,tempHeader[i])
        

    for i in range(0,len(Genus)):
        worksheet.write(nRowStart+i+1,0,Genus[i])
        lGenusIndices = [k for k,x in enumerate(lGenus) if x == Genus[i]]

        worksheet.write(nRowStart+i+1,1,c.get(Genus[i]))
        worksheet.write(nRowStart+i+1,2,c.get(Genus[i])/nItem)
        worksheet.write(nRowStart+i+1,3,sum([lContigLength[j] for j in lGenusIndices]))
    
    nRowStart = nRowStart + i + 5    
    
    c= Counter(lSpecies)

    nItem = len(lSpecies)
    Species = list(c.keys())
        
    for i in range(0,3):    
        worksheet.write(nRowStart,i,tempHeader[i])
    
    for i in range(0,len(Species)):
        worksheet.write(nRowStart+i+1,0,Species[i])
        lSpeciesIndices = [k for k,x in enumerate(lSpecies) if x == Species[i]]

        worksheet.write(nRowStart+i+1,1,c.get(Species[i]))
        worksheet.write(nRowStart+i+1,2,c.get(Species[i])/nItem)
        worksheet.write(nRowStart+i+1,3,sum([lContigLength[j] for j in lSpeciesIndices]))

  #Summary       
    worksheet = workbook.add_worksheet(dbName +'_Summary')
    
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
    lContigLength = [i[6] for i in SummaryTbl]

    for i in lSciNames :
        temp = i.replace('[','')
        temp = temp.replace(']','')
        lGenus.append(temp.split()[0])
        lSpecies.append(temp.split()[0] + ' ' + temp.split()[1])
        
    c= Counter(lGenus)
    
    nItem = len(lGenus)
    Genus = list(c.keys())
        
    for i in range(0,4):    
        worksheet.write(nRowStart,i,tempHeader[i])
    
    for i in range(0,len(Genus)):
        worksheet.write(nRowStart+i,0,Genus[i])
        lGenusIndices = [k for k,x in enumerate(lGenus) if x == Genus[i]]

        worksheet.write(nRowStart+i+1,1,c.get(Genus[i]))
        worksheet.write(nRowStart+i+1,2,c.get(Genus[i])/nItem)
        worksheet.write(nRowStart+i+1,3,sum([lContigLength[j] for j in lGenusIndices]))

    nRowStart = nRowStart + i + 5    
    
    c= Counter(lSpecies)

    nItem = len(lSpecies)
    Species = list(c.keys())
    
    
    for i in range(0,len(Species)):
        worksheet.write(nRowStart+i+1,0,Species[i])
        lSpeciesIndices = [k for k,x in enumerate(lSpecies) if x == Species[i]]

        worksheet.write(nRowStart+i+1,1,c.get(Species[i]))
        worksheet.write(nRowStart+i+1,2,c.get(Species[i])/nItem)  
        worksheet.write(nRowStart+i+1,3,sum([lContigLength[j] for j in lSpeciesIndices]))
