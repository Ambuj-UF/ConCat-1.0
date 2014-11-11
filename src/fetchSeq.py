# -*- coding: utf-8 -*-

################################################################################################################
# Codon sequence extraction program                                                                            #
#                                                                                                              #
# Copyright (C) {2014}  {Ambuj Kumar, Kimball-Brain lab group, Biology Department, University of Florida}      #
#                                                                                                              #
# This program is free software: you can redistribute it and/or modify                                         #
# it under the terms of the GNU General Public License as published by                                         #
# the Free Software Foundation, either version 3 of the License, or                                            #
# (at your option) any later version.                                                                          #
#                                                                                                              #
# This program is distributed in the hope that it will be useful,                                              #
# but WITHOUT ANY WARRANTY; without even the implied warranty of                                               #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                                                #
# GNU General Public License for more details.                                                                 #
#                                                                                                              #
# This program comes with ABSOLUTELY NO WARRANTY;                                                              #
# This is free software, and you are welcome to redistribute it                                                #
# under certain conditions;                                                                                    #
#                                                                                                              #
################################################################################################################

import os
import sys
import urllib2

try:
    from Bio import Entrez
    from Bio import SeqIO
except ImportError, e:
    sys.exit("Biopython not found")


def xmlcreate(ID):
    url = "http://www.ncbi.nlm.nih.gov/gene/" + ID
    xmlData = urllib2.urlopen(url)
    with open("export.xml",'wb+') as fp:
        for line in xmlData:
            fp.write(line)


def xmlparser():
    IDs = list()
    handle = open('export.xml', 'r').readlines()
    for lines in handle:
        if '→' in lines and "/nuccore/" in lines and 'Location' not in lines:
            IDs.append(lines.split('→')[0].split('/')[2].split('"')[0])
    
    return IDs


def cdsExt(ID, geneName):
    retdata = Entrez.efetch(db="nucleotide", id=ID, rettype='gb', retmode='text').read()
    with open("Align/" + geneName.split('.')[0] + ".log", "a") as fp:
        if 'LOW QUALITY PROTEIN' in retdata:
            fp.write('%s CDS is of low quality\n' %ID)
    
    data = retdata.split('\n')
    for obj in data:
        if '     CDS             ' in obj:
            try:
                cdsRange = [int(obj.lstrip('     CDS             ').split('..')[0]), int(obj.lstrip('     CDS             ').split('..')[1])]
            except ValueError:
                try:
                    cdsRange = [int(obj.lstrip('     CDS             ').split('..')[0].lstrip('<')), int(obj.lstrip('     CDS             ').split('..')[1])]
                except ValueError:
                    try:
                        cdsRange = [int(obj.lstrip('     CDS             ').split('..')[0]), int(obj.lstrip('     CDS             ').split('..')[1].lstrip('>'))]
                    except ValueError:
                        print("Problem found while extracting cds from %s. Please report this issue to ambuj (at) ufl (dot) edu" %obj)
                        continue

    recData = Entrez.efetch(db="nucleotide", id=ID, rettype="gb", warning=False)
    record = SeqIO.read(recData, 'genbank')
    record.seq = record.seq[cdsRange[0]-1:cdsRange[1]]

    return record


def mrnaExt(ID):
    recData = Entrez.efetch(db="nucleotide", id=ID, rettype="gb", warning=False)
    record = SeqIO.read(recData, 'genbank')
    return record


def cdsImport(geneName, group):
    
    inpTerm = geneName + "[sym] AND " + group + "[orgn]"
    Entrez.email = 'sendambuj@gmail.com'
    
    print("Importing CDS sequences for %s gene" %geneName)
    handle = Entrez.esearch(db="gene", term=inpTerm, rettype='xml', RetMax=300, warning=False)
    records = Entrez.read(handle)
    idList = records["IdList"]
    
    outRecord = list()
    
    for ids in idList:
        xmlcreate(ids)
        refIds = xmlparser()
        os.remove('export.xml')
        recordList = list()
        for inIDs in refIds:
            recordList.append(cdsExt(inIDs, geneName))

        try:
            longestRec = recordList[0]
        except:
            continue
        for rec in recordList:
            longestRec = rec if len(rec.seq) > len(longestRec.seq) else longestRec
        print(longestRec.description)
        outRecord.append(longestRec)

    with open("Align/" + geneName + '.fas', 'w') as fp:
        SeqIO.write(outRecord, fp, 'fasta')

    fdata = open("Align/" + geneName + '.fas', 'r').readlines()
    with open("Align/" + geneName + '.fas', 'w') as fp:
        for lines in fdata:
            if '>' in lines and 'PREDICTED' in lines:
                newLine = '>' + lines.split(' ')[2] + '_' + lines.split(' ')[3] + '|' + lines.split(' ')[0].lstrip('>')
                fp.write('%s\n'%newLine)
            elif '>' in lines and 'PREDICTED' not in lines:
                newLine = '>' + lines.split(' ')[1] + '_' + lines.split(' ')[2] + '|' + lines.split(' ')[0].lstrip('>')
                fp.write('%s\n'%newLine)
            else:
                fp.write('%s'%lines)



def mrnaImport(geneName, group):
    inpTerm = geneName + "[sym] AND " + group + "[orgn]"
    Entrez.email = 'sendambuj@gmail.com'
    
    handle = Entrez.esearch(db="gene", term=inpTerm, rettype='xml', RetMax=300, silent=True)
    records = Entrez.read(handle)
    idList = records["IdList"]
    
    outRecord = list()
    for ids in idList:
        xmlcreate(ids)
        refIds = xmlparser()
        os.remove('export.xml')
        recordList = list()
        for inIDs in refIds:
            recordList.append(mrnaExt(inIDs))
        
        try:
            longestRec = recordList[0]
        except:
            continue
        for rec in recordList:
            longestRec = rec if len(rec.seq) > len(longestRec.seq) else longestRec
        print(longestRec.description)
        outRecord.append(longestRec)
    
    with open("Align/" + geneName + '.fas', 'w') as fp:
        SeqIO.write(outRecord, fp, 'fasta')
    
    fdata = open("Align/" + geneName + '.fas', 'r').readlines()
    with open("Align/" + geneName + '.fas', 'w') as fp:
        for lines in fdata:
            if '>' in lines and 'PREDICTED' in lines:
                newLine = '>' + lines.split(' ')[2] + '_' + lines.split(' ')[3] + '|' + lines.split(' ')[0].lstrip('>')
                fp.write('%s\n'%newLine)
            elif '>' in lines and 'PREDICTED' not in lines:
                newLine = '>' + lines.split(' ')[1] + '_' + lines.split(' ')[2] + '|' + lines.split(' ')[0].lstrip('>')
                fp.write('%s\n'%newLine)
            else:
                fp.write('%s'%lines)




def oneGeneCdsImport(geneName, group):

    inpTerm = geneName + "[sym] AND " + group + "[orgn]"
    Entrez.email = 'sendambuj@gmail.com'
        
    print("Importing %s %s gene CDS sequence" %(group, geneName))
    handle = Entrez.esearch(db="gene", term=inpTerm, rettype='xml', RetMax=300, warning=False)
    records = Entrez.read(handle)
    ids = records["IdList"][0]
    
    xmlcreate(ids)
    refIds = xmlparser()
    os.remove('export.xml')
    recordList = list()
    for inIDs in refIds:
        recordList.append(cdsExt(inIDs, geneName))
        
    longestRec = recordList[0]
    for rec in recordList:
        longestRec = rec if len(rec.seq) > len(longestRec.seq) else longestRec

    longestRec.id = geneName

    return longestRec


def oneGeneMrnaImport(geneName, group):
    
    inpTerm = geneName + "[sym] AND " + group + "[orgn]"
    Entrez.email = 'sendambuj@gmail.com'
    
    print("Importing %S %S gene sequences" %(group, geneName))
    handle = Entrez.esearch(db="gene", term=inpTerm, rettype='xml', RetMax=300, warning=False)
    records = Entrez.read(handle)
    ids = records["IdList"][0]
    
    xmlcreate(ids)
    refIds = xmlparser()
    os.remove('export.xml')
    recordList = list()
    for inIDs in refIds:
        recordList.append(mrnaExt(inIDs))
    
    longestRec = recordList[0]
    for rec in recordList:
        longestRec = rec if len(rec.seq) > len(longestRec.seq) else longestRec
    
    return longestRec




