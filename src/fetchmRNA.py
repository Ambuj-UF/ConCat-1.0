# -*- coding: utf-8 -*-

################################################################################################################
# MRNA sequence extraction program                                                                             #
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



def xmlcreateMRNA(ID):
    url = "http://www.ncbi.nlm.nih.gov/gene/" + ID
    xmlData = urllib2.urlopen(url)
    with open("export.xml",'wb+') as fp:
        for line in xmlData:
            fp.write(line)

def xmlparserMRNA():
    IDs = list()
    handle = open('export.xml', 'r').readlines()
    for lines in handle:
        if '→' in lines and 'Location' not in lines:
            IDs.append(lines.split('→')[0].split('/')[2].split('"')[0])
    
    return IDs

def mrnaExt(ID, geneName):
    retdata = Entrez.efetch(db="nucleotide", id=ID, rettype='gb', retmode='text').read()
    data = retdata.split('\n')
    recData = Entrez.efetch(db="nucleotide", id=ID, rettype="gb")
    record = rec = SeqIO.read(recData, 'genbank')
    return record

def fetchMRNA(geneName, group):
    inpTerm = geneName + "[sym] AND " + group + "[orgn]"
    Entrez.email = 'sendambuj@gmail.com'

    handle = Entrez.esearch(db="gene", term=inpTerm, rettype='xml', RetMax=300)
    records = Entrez.read(handle)
    idList = records["IdList"]

    outRecord = list()
    for ids in idList:
        xmlcreateMRNA(ids)
        refIds = xmlparserMRNA()
        os.remove('export.xml')
        recordList = list()
        for inIDs in refIds:
            recordList.append(mrnaExt(inIDs, geneName))

        try:
            longestRec = recordList[0]
        except:
            continue
        for rec in recordList:
            longestRec = rec if len(rec.seq) > len(longestRec.seq) else longestRec
        print longestRec.description
        outRecord.append(longestRec)

    with open("../Align/" + geneName + '.fas', 'w') as fp:
        SeqIO.write(outRecord, fp, 'fasta')

    fdata = open("../Align/" + geneName + '.fas', 'r').readlines()
    with open("../Align/" + geneName + '.fas', 'w') as fp:
        for lines in fdata:
            if '>' in lines and 'PREDICTED' in lines:
                newLine = '>' + lines.split(' ')[3] + '_' + lines.split(' ')[2] + '|' + lines.split(' ')[0].lstrip('>')
                fp.write('%s\n'%newLine)
            elif '>' in lines and 'PREDICTED' not in lines:
                newLine = '>' + lines.split(' ')[2] + '_' + lines.split(' ')[1] + '|' + lines.split(' ')[0].lstrip('>')
                fp.write('%s\n'%newLine)
            else:
                fp.write('%s'%lines)














