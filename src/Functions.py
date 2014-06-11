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
import math
import functools


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
            Position = posMatrix[0]
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
        cutVal = listVal[-1]
        for j, inval in enumerate(listVal):
            if j < len(listVal) - 1:
                if listVal[j+1] - listVal[j] > cutOff:
                    cutVal = listVal[j]
                    break
                                         
        listPos = [[x, val] for x, val in OVdict.items() if val > cutVal]
                                         
    return listPos


def binAll(rcvRange,
           entropyRange,
           combined,
           RCVdict,
           entropyDict,
           gcDict,
           gcRange
           ):
    lineListRcv = []
    lineListEntropy = []
    lineListGC = []
    if rcvRange != None:
        rStart = float(rcvRange.split('-')[0])
        rEnd = float(rcvRange.split('-')[1])
        
        for key, val in RNAdict.items():
            try:
                sink = entropyDict[key]
            except KeyError:
                entropyDict[key] = (['NA'])
            try:
                sink = gcDict[key]
            except KeyError:
                gcDict[key] = (['NA'])
            if val >= rStart and val <= rEnd:
                try:
                    lineListRcv.append("BIN_RCV %s = %s-%s [RCV Score = %s] [Entropy = %s] [GC Content (in percentage) = %s]" %(key, combined.charsets[key][0], combined.charsets[key][-1], val, entropyDict[key], gcDict[key]))
                except KeyError:
                    continue

    if entropyRange != None:
        rStart = float(entropyRange.split('-')[0])
        rEnd = float(entropyRange.split('-')[1])

        for key, val in entropyDict.items():
            try:
                sink = RCVdict[key]
            except KeyError:
                RCVdict[key] = (['NA'])
            try:
                sink = gcDict[key]
            except KeyError:
                gcDict[key] = (['NA'])

            if val >= rStart and val <= rEnd:
                try:
                    lineListEntropy.append("BIN_Entropy %s = %s-%s [RCV Score = %s] [Entropy = %s] [GC Content (in percentage) = %s]" %(key, combined.charsets[key][0]+1, combined.charsets[key][-1]+1, RCVdict[key], val, gcDict[key]))
                except KeyError:
                    continue

    if gcRange != None:
        rStart = float(gcRange.split('-')[0])
        rEnd = float(gcRange.split('-')[1])
            try:
                sink = entropyDict[key]
            except KeyError:
                entropyDict[key] = (['NA'])
            try:
                sink = RCVdict[key]
            except KeyError:
                RCVdict[key] = (['NA'])

        for key, val in gcDict.items():
            if val >= rStart and val <= rEnd:
                try:
                    lineListGC.append("BIN_GC %s = %s-%s [RCV Score = %s] [Entropy = %s] [GC Content (in percentage) = %s]" %(key, combined.charsets[key][0]+1, combined.charsets[key][-1]+1, RCVdict[key], entropyDict[key], val))
                except KeyError:
                    continue

    return [lineListRcv, lineListEntropy, lineListGC]



def GCcontent(combined):
    GCdict = dict()
    msa = MultipleSeqAlignment(NexusHandler(1).combineToRecord(combined))
    for key, val in combined.charsets.items():
            gcCount = 0
            try:
                msaGene = msa[:, combined.charsets[key][0]:combined.charsets[key][-1]]
                for inval in msaGene:
                    gcCount = gcCount + inval.seq.count('G') + inval.seq.count('C')
                GCdict[key] = (float(gcCount)/(len(msaGene)*len(msaGene[1]))*100)
            except KeyError:
                continue

    return GCdict





def percentile(N, percent, key=lambda x:x):
    """
        Find the percentile of a list of values.
        
        @parameter N - is a list of values. Note N MUST BE already sorted.
        @parameter percent - a float value from 0.0 to 1.0.
        @parameter key - optional key function to compute value from each element of N.
        
        @return - the percentile of the values
        """
    if not N:
        return None
    k = (len(N)-1) * percent
    f = math.floor(k)
    c = math.ceil(k)
    if f == c:
        return key(N[int(k)])
    d0 = key(N[int(f)]) * (c-k)
    d1 = key(N[int(c)]) * (k-f)
    return d0+d1



def gcEntropyPer(RCVdict,
                 entropyDict,
                 gcDict,
                 combined,
                 which
                 ):
    """
        Finds the values grouped in percentile from input dictionary data.
        
        @parameter RCVdict - is RCV dictionary element.
        @parameter entropyDict - is Entropy dictionary element.
        @parameter gcDict - is GC content dictionary element.
        @parameter combined - is a 3D nexus data matrix
        @parameter which - takes dictionary type as input (It could be GC, RCV or Entropy dictionary)
        
        """

    if which == 'rcv':
        supDict = RCVdict
    
    elif which == 'ent':
        supDict = entropyDict

    else:
        supDict = gcDict

    for key, val in supDict.items():
        supList.append(val)

    sup25 = percentile(supList.sort(), 0.25)
    sup75 = percentile(supList.sort(), 0.75)

    v0to25 = [x for x in supList if x >= 0 and x < sup25]
    v25to75 = [x for x in supList if x >= ent25 and x < sup75]
    v75to100 = [x for x in supList if x >= sup75]

    for key, val in supDict.items():
        try:
            sink = RCVdict[key]
        except KeyError:
            RCVdict[key] = (['NA'])
        try:
            sink = entropyDict[key]
        except KeyError:
            entropyDict[key] = (['NA'])
        try:
            sink = gcDict[key]
        except KeyError:
            gcDict[key] = (['NA'])

        if val in v0to25:
            entPlist0to25.append("BIN_RCV %s = %s-%s [RCV Score = %s] [Entropy = %s] [GC Content (in percentage) = %s]"\
                                 %(key, combined.charsets[key][0]+1, combined.charsets[key][-1]+1, RCVdict[key], entropyDict[key], gcDict[key]))
        elif val in v25to75:
            entPlist25to75.append("BIN_RCV %s = %s-%s [RCV Score = %s] [Entropy = %s] [GC Content (in percentage) = %s]"\
                                  %(key, combined.charsets[key][0]+1, combined.charsets[key][-1]+1, RCVdict[key], entropyDict[key], gcDict[key]))
        else:
            entPlist75to100.append("BIN_RCV %s = %s-%s [RCV Score = %s] [Entropy = %s] [GC Content (in percentage) = %s]"\
                                   %(key, combined.charsets[key][0]+1, combined.charsets[key][-1]+1, RCVdict[key], entropyDict[key], gcDict[key]))

    supPdict['0_to_25'] = (entPlist0to25)
    supPdict['25_to_75'] = (entPlist25to75)
    supPdict['75_to_100'] = (entPlist75to100)

    return supPdict



def binPercent(RCVdict, entropyDict, gcDict, combined):
    """
        returns the list of dictionaries that has charset data grouped under percentile bins.
        
        @parameter RCVdict - is RCV dictionary element.
        @parameter entropyDict - is Entropy dictionary element.
        @parameter gcDict - is GC content dictionary element.
        @parameter combined - is a 3D nexus data matrix
        
        """
    
    for key, val in RCVdict:
        RCVdict[key] = (val[0])

    rcvPdict = gcEntropyPer(RCVdict,
                            entropyDict,
                            gcDict,
                            combined,
                            which='rcv'
                            )
    
    entPdict = gcEntropyPer(RCVdict,
                            entropyDict,
                            gcDict,
                            combined,
                            which='ent'
                            )
                            
    gcPdict = gcEntropyPer(RCVdict,
                           entropyDict,
                           gcDict,
                           combined,
                           which='gc'
                           )

    return [rcvPdict, entPdict, gcPdict]







