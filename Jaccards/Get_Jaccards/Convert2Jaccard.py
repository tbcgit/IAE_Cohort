#-*- coding: utf-8 -*-
#!/usr/bin/env python
#AQUEST CODI NO VA BÉ DEL TOT, La segona linea de l'output sobra i s'ha de llevar sencera MANUALMENT

def Reduce(line): 
    
    lineLIST = line.split('\t')
    lineOut=[]
    for i in range(len(lineLIST)):
        if i == 0: #per evitar que procese el primer element, que és el nom del taxon s/SILVA
            lineOut.append(lineLIST[0])
        else:
            if lineLIST[i] in [0,'0']:
                lineOut.append('0')
            else:
                lineOut.append('1')
	#print lineOut

    return lineOut


inputName ='TaxSilva_GEN.csv'
infile=open(inputName ,'r')

contador=0
lines2write=[]
for linea in infile:
    if contador == 0:
        lines2write.append(linea)
        contador=contador + 1
    
    if contador > 0:
        newlineREDUCED = Reduce(linea)
        lines2write.append(newlineREDUCED)
        contador = contador + 1

#print lines2write
outputName=inputName.split('.')[0]+'_transformed.tsv'
output=open(outputName,'w')

contador=0
for linea in lines2write:
    if contador == 0:
        output.write(linea)
        contador = contador + 1
    else:
        line2write = '\t'.join(linea) + '\n'
        output.write(line2write)
        contador = contador + 1

print 'Transformació OK'

output.close()


        
