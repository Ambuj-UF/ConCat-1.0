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
from Bio import AlignIO
from Handler import *
import operator


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
        numA = numA + val.seq.count('A') + val.seq.count('a')
        numC = numC + val.seq.count('C') + val.seq.count('c')
        numG = numG + val.seq.count('G') + val.seq.count('g')
        numT = numT + val.seq.count('T') + val.seq.count('t')

    countDict = dict()
    for i, val in enumerate(combine):
        countDict[combine[i].id] = ([val.seq.count('A') + val.seq.count('a'), val.seq.count('C') + val.seq.count('c'), val.seq.count('G') + val.seq.count('g'), val.seq.count('T') + val.seq.count('t')])

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


def RCVprotCal(combine): # combined is a multiple sequence alignment object
    
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
    numB = 0
    numI = 0
    numL = 0
    numF = 0
    numN = 0
    numS = 0
    
    for i, val in enumerate(combine):
        numA = numA + val.seq.count('D') + val.seq.count('d') + val.seq.count('E') + val.seq.count('e')
        numB = numB + val.seq.count('R') + val.seq.count('r') + val.seq.count('K') + val.seq.count('k')
        numI = numI + val.seq.count('I') + val.seq.count('i') + val.seq.count('V') + val.seq.count('v')
        numL = numL + val.seq.count('L') + val.seq.count('l') + val.seq.count('M') + val.seq.count('m')
        numF = numF + val.seq.count('F') + val.seq.count('f') + val.seq.count('W') + val.seq.count('w') + val.seq.count('Y') + val.seq.count('y')
        numN = numN + val.seq.count('N') + val.seq.count('n') + val.seq.count('Q') + val.seq.count('q')
        numS = numS + val.seq.count('S') + val.seq.count('s') + val.seq.count('T') + val.seq.count('t')
    
    countDict = dict()
    for i, val in enumerate(combine):
        countDict[combine[i].id] = ([val.seq.count('D') + val.seq.count('d') + val.seq.count('E') + val.seq.count('e'), \
                                     val.seq.count('R') + val.seq.count('r') + val.seq.count('K') + val.seq.count('k'), \
                                     val.seq.count('I') + val.seq.count('i') + val.seq.count('V') + val.seq.count('v'), \
                                     val.seq.count('L') + val.seq.count('l') + val.seq.count('M') + val.seq.count('m'), \
                                     val.seq.count('F') + val.seq.count('f') + val.seq.count('W') + val.seq.count('w') + val.seq.count('Y') + val.seq.count('y'), \
                                     val.seq.count('N') + val.seq.count('n') + val.seq.count('Q') + val.seq.count('q'), \
                                     val.seq.count('S') + val.seq.count('s') + val.seq.count('T') + val.seq.count('t')])
    
    def abs(number):
        if number > 0 or number == 0:
            pass
        else:
            number = - number
        return number
    
    rcvCalc = 0
    
    nTaxa = len(combine)
    
    for key, val in countDict.items():
        rcvCalc = rcvCalc + abs(val[0] - (numA/nTaxa)) + \
                            abs(val[1] - (numB/nTaxa)) + \
                            abs(val[2] - (numI/nTaxa)) + \
                            abs(val[3] - (numL/nTaxa)) + \
                            abs(val[4] - (numF/nTaxa)) + \
                            abs(val[5] - (numN/nTaxa)) + \
                            abs(val[6] - (numS/nTaxa))
    
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

    print "Final output saved in %s" %filename.split('.')[0] + '.' + formDict[output].split('.')[1]

def ConvertAll(inp_format):
    os.chdir('Input')
    files = glob.glob("*.*")
    
    if inp_format == 'fasta':
        for filename in files:
            try:
                alignment = AlignIO.read(open(filename), "fasta", alphabet=Gapped(IUPAC.protein))
                g = open(filename.split(".")[0] + '.nex', 'w')
                g.write(alignment.format("nexus"))
                g.close()
            
            except ValueError:
                continue
    
    else:
        for filename in files:
            try:
                handle = open(filename, 'rU')
                record = list(SeqIO.parse(handle, inp_format))
                fp = open(filename.split('.')[0] + '.nex', 'w')
                SeqIO.write(record, fp, "nexus")
                fp.close()
                handle.close()
            except:
                print "Bad Alignment %s\n" %filename

    os.chdir('..')

def fastEvol(combined, cutOff):
    if cutOff == None:
        listPos = []
    else:
        msa = MultipleSeqAlignment(NexusHandler(1).combineToRecord(combined))
        charList = []
        i = 1
        while i <= len(msa[0]):
            j = 1
            tempList = []
            while j <= len(msa):
                tempList.append(msa[j-1][i-1])
                j = j + 1
            charList.append(tempList)
            i = i + 1

        OVdict = dict()
        for i, val in enumerate(charList):
            val = list(''.join(val).replace('?', ''))
            posVal = []

            if len(set(val)) == 1:
                posVal.append(0)

            else:
                outCounter = 0
                inCounter = 1
                for inval in val:
                    inCounter = 1 + outCounter
                    while inCounter < len(val):
                        if inval == val[inCounter]:
                            posVal.append(0)
                        else:
                            posVal.append(1)
                        inCounter = inCounter + 1
            
                    outCounter = outCounter + 1


            k = (math.pow(len(val), 2) - len(val))/2

            OVdict['Position_%s' %i] = (sum(posVal)/k)
            
            listVal = sorted(OVdict.values())
            for j, val in enumerate(listVal):
                if j < len(listVal) - 1:
                    if listVal[j+1] - listVal[j] > cutOff:
                        cutVal = listVal[i]                    
                    else:
                        cutVal = listVal[-1]
                                         
        listPos = [[x, val] for x, val in OVdict.items() if val > cutVal]
                                         
    return listPos


















