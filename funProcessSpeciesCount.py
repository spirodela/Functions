#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct  6 12:43:13 2017

@author: yi.yan
"""
def funProcessSpeciesCount(sFileName):
    import xlsxwriter
    import math 
    from scipy.stats import t 
    
    alpha = 1 - 0.99
    CriticalProbability = 1 - alpha/2 
     
    
    #sFileName = 'PFDA1_20170908_M00708_IL100092259_TAATGCG_L001_R1_trimmed_10000_Count.txt'
    lTable = []
    Header = ('Specie Name','Tax ID','Count','Percentage','Total Hit','99% Confidence Interval')
    file_object  = open(sFileName, 'r')   
    
    TotalHit = 0
    for line in file_object.readlines(): 
        temp = line.split() 
        TotalHit = TotalHit + int(temp[0]) 
    
    file_object.close()
    
    DegreesFreedom = TotalHit - 1 
    tStat = t.ppf(CriticalProbability,DegreesFreedom)
    file_object  = open(sFileName, 'r')   
    for line in file_object.readlines(): 
        temp = line.split() 
        SpecieName = temp[1] 
        for j in range(2,len(temp)-1):
            SpecieName = SpecieName + ' ' + temp[j]
        TaxID = temp[-1] 
        Count = int(temp[0])
        Percentage = Count/TotalHit 
        SE = math.sqrt(Percentage*(1-Percentage)/TotalHit)
        CI = SE*tStat*100 
        Percentage = Percentage*100 
        tempTup = (SpecieName,TaxID,Count,Percentage,TotalHit,CI)
        lTable.append(tempTup)
    
    s = sorted(lTable,key=lambda x:(x[2]),reverse=True)
    
    workbook = xlsxwriter.Workbook(sFileName[0:-4]+'.xlsx')
    worksheet = workbook.add_worksheet()
    
    for j in range(0,len(Header)):
        worksheet.write(0,j,Header[j])
        
    for i in range(0,len(s)):
        for j in range(0,len(s[i])):
            worksheet.write(i+1,j,s[i][j])
    
    workbook.close()