################################################################################################################
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

from Bio.Alphabet import IUPAC, Gapped
from Bio import SeqIO
import glob
import os


def RCVcal(combine): # combined is a multiple sequence alignment object
    
    # This section removes all the conserved sites from multiple sequence alignment object
    similarityCount = {}
    posMatrix = []
    n=0
    while n < len(combine[0]):
        similarityCount[n] = 0
        n = n + 1
    i = 1
    while i <= len(combine[0]):
        j = 1
        while j <= len(combine):
            character = combine[j-1][i-1]
            if character == combine[0][i-1]:
                similarityCount[i-1] = similarityCount[i-1] + 1
            j = j + 1
        if similarityCount[i-1] == len(combine):
            posMatrix.append(i)
        i = i + 1

    cycles = 0

    while cycles < len(posMatrix)-1:
        if cycles == 0:
            Position = posMatrix[1]
            termPos =  Position - len(combine[1])
            varFirst = Position-1
            combine = combine[:, :varFirst] + combine[:, termPos:]
        else:
            Position = posMatrix[cycles+1]
            termPos = (Position-cycles) - len(combine[1])
            varFirst = Position-cycles-1
            combine = combine[:, :varFirst] + combine[:, termPos:]
        cycles = cycles + 1

    ########################################################################################
    # RCV calculation begins here

    numA = 0
    numC = 0
    numG = 0
    numT = 0

    for i, val in enumerate(combine):
        numA = numA + val.seq.count('A')
        numC = numC + val.seq.count('C')
        numG = numG + val.seq.count('G')
        numT = numT + val.seq.count('T')

    countDict = dict()
    for i, val in enumerate(combine):
        countDict[combine[i].id] = ([val.seq.count('A'), val.seq.count('C'), val.seq.count('G'), val.seq.count('T')])

    def abs(number):
        if number > 0 or number == 0:
            pass
        else:
            number = - number
        return number

    rcvCalc = 0

    nTaxa = len(combine)

    for key, val in countDict.items():
        rcvCalc = rcvCalc + abs(val[0] - (numA/nTaxa)) + abs(val[1] - (numC/nTaxa)) + abs(val[2] - (numG/nTaxa)) + abs(val[3] - (numT/nTaxa))

    totalRCV = float(rcvCalc)/float(len(combine)*len(combine[0]))

    return totalRCV


def Convert(input, output, filename):
    formDict = {
        'fasta': '*.fas',
        'nexus': '*.nex',
        'phylip': '*.phy',
        'phylip-sequential': '*.phy',
        'phylip-relaxed': '*.phy'
    }
    
    if input == 'fasta' and output == 'nexus':
        alignment = AlignIO.read(open(filename), "fasta", alphabet=Gapped(IUPAC.protein))
        g = open(filename.split(".")[0] + '.nex', 'w')
        g.write (alignment.format("nexus"))
        g.close()
    
    else:
        try:
            handle = open(filename, 'rU')
            record = list(SeqIO.parse(handle, input))
            fp = open(filename.split('.')[0] + '.' + formDict[output].split('.')[1], 'w')
            SeqIO.write(record, fp, output)
            fp.close()
            handle.close()
        except:
            print "Bad Alignment\n"


def ConvertAll(input):
    formDict = {
        'fasta': '*.fas',
        'nexus': '*.nex',
        'phylip': '*.phy',
        'phylip-sequential': '*.phy',
        'phylip-relaxed': '*.phy'
    }
    os.chdir('Input')
    files = glob.glob(formDict[input])
    
    if input == 'fasta' and output == 'nexus':
        for filename in files:
            alignment = AlignIO.read(open(filename), "fasta", alphabet=Gapped(IUPAC.protein))
            g = open(filename.split(".")[0] + '.nex', 'w')
            g.write (alignment.format("nexus"))
            g.close()
    
    else:
        for filename in files:
            try:
                handle = open(filename, 'rU')
                record = list(SeqIO.parse(handle, input))
                fp = open(filename.split('.')[0] + '.nex', 'w')
                SeqIO.write(record, fp, "nexus")
                fp.close()
                handle.close()
            except:
                print "Bad Alignment %s\n" %filename

    os.chdir('..')



