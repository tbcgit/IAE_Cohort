#-*- coding: utf-8 -*-
#!/usr/bin/env python

# script split in unique files by sample a table of ASV counts

def DepurateLineBreak(ListFromFile):
    NewList=[]
    
    for element in ListFromFile:
        if '\n' in element:
            NewElement=element.replace('\n','')
            NewList.append(NewElement)
        else:
            NewList.append(element)
            
    return(NewList)
            


infile = open('ASVs_counts_DEF.tsv','r')
AllLines = infile.readlines()
SamplesList = DepurateLineBreak(AllLines[0].split('\t'))


for sample in SamplesList[1:]:
    sampleMatrix = []
    index = SamplesList.index(sample)
    Filename = str(sample)+'.tsv'
    output = open(Filename,'w')
    output.write('TAXONOMY'+'\t'+str(sample)+'\n')
    for line in AllLines[1:]:
        count = line.split('\t')[index]
        sampleMatrix.append([line.split('\t')[0], count])
    for pair in sampleMatrix:
        output.write(str(pair[0])+'\t'+str(pair[1])+'\n')
    output.close()
        
    
infile.close()