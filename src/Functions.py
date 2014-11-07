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

import sys
import glob
import os
import operator
import math
import functools

try:
    import matplotlib.mlab as mlab
    import matplotlib.pyplot as plt
except:
    sys.exit("ConCat requires 1.8 or higher version of Numpy")

from array import array
from Bio import AlignIO
from Handler import *
from Bio.Alphabet import IUPAC, Gapped
from Bio import SeqIO



def RCVcal(combine):
    """
        Calculates RCV values for DNA alignment
        @parameter combine - combine is a multiple sequence alignment object
        Returns - Dictionary element with RCV values and their corresponding alignment IDs
        """
    # This section removes all the conserved sites from multiple sequence alignment object
    similarityCount = {}; posMatrix = []; n=0
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

    numA = 0; numC = 0; numG = 0; numT = 0

    for i, val in enumerate(combine):
        numA = numA + val.seq.count('A') + val.seq.count('a')
        numC = numC + val.seq.count('C') + val.seq.count('c')
        numG = numG + val.seq.count('G') + val.seq.count('g')
        numT = numT + val.seq.count('T') + val.seq.count('t')

    countDict = dict()
    for i, val in enumerate(combine):
        countDict[combine[i].id] = ([val.seq.count('A') + val.seq.count('a'), val.seq.count('C') + val.seq.count('c'),\
                                     val.seq.count('G') + val.seq.count('g'), val.seq.count('T') + val.seq.count('t')])

    def abs(number):
        if number > 0 or number == 0:
            pass
        else:
            number = - number
        return number

    rcvCalc = 0; nTaxa = len(combine)
    for key, val in countDict.items():
        rcvCalc = rcvCalc + abs(val[0] - (numA/nTaxa)) + abs(val[1] - (numC/nTaxa)) + abs(val[2] - (numG/nTaxa)) + abs(val[3] - (numT/nTaxa))

    totalRCV = float(rcvCalc)/float(len(combine)*len(combine[0]))
    return totalRCV


def RCVprotCal(combine):
    """ 
        Calculates RCV values for protein alignment
        @parameter combine - combine is a multiple sequence alignment object
        Returns - Dictionary element with RCV values and their corresponding alignment IDs
        """
    # This section removes all the conserved sites from multiple sequence alignment object
    similarityCount = {}; posMatrix = []; n=0
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
    
    numA = 0; numB = 0; numI = 0; numL = 0; numF = 0; numN = 0; numS = 0; numC = 0; numH = 0; numU = 0; numG = 0; numP = 0
    for i, val in enumerate(combine):
        numA = numA + val.seq.count('D') + val.seq.count('d') + val.seq.count('E') + val.seq.count('e')
        numB = numB + val.seq.count('R') + val.seq.count('r') + val.seq.count('K') + val.seq.count('k')
        numI = numI + val.seq.count('I') + val.seq.count('i') + val.seq.count('V') + val.seq.count('v')
        numL = numL + val.seq.count('L') + val.seq.count('l') + val.seq.count('M') + val.seq.count('m')
        numF = numF + val.seq.count('F') + val.seq.count('f') + val.seq.count('W') + val.seq.count('w') + val.seq.count('Y') + val.seq.count('y')
        numN = numN + val.seq.count('N') + val.seq.count('n') + val.seq.count('Q') + val.seq.count('q')
        numS = numS + val.seq.count('S') + val.seq.count('s') + val.seq.count('T') + val.seq.count('t')
        numC = numC + val.seq.count('C') + val.seq.count('c')
        numC = numH + val.seq.count('H') + val.seq.count('h')
        numU = numU + val.seq.count('A') + val.seq.count('A')
        numG = numG + val.seq.count('G') + val.seq.count('g')
        numP = numP + val.seq.count('P') + val.seq.count('p')
    
    countDict = dict()
    for i, val in enumerate(combine):
        countDict[combine[i].id] = ([val.seq.count('D') + val.seq.count('d') + val.seq.count('E') + val.seq.count('e'), \
                                     val.seq.count('R') + val.seq.count('r') + val.seq.count('K') + val.seq.count('k'), \
                                     val.seq.count('I') + val.seq.count('i') + val.seq.count('V') + val.seq.count('v'), \
                                     val.seq.count('L') + val.seq.count('l') + val.seq.count('M') + val.seq.count('m'), \
                                     val.seq.count('F') + val.seq.count('f') + val.seq.count('W') + val.seq.count('w') + val.seq.count('Y') + val.seq.count('y'), \
                                     val.seq.count('N') + val.seq.count('n') + val.seq.count('Q') + val.seq.count('q'), \
                                     val.seq.count('S') + val.seq.count('s') + val.seq.count('T') + val.seq.count('t'),\
                                     val.seq.count('C') + val.seq.count('c'),\
                                     val.seq.count('H') + val.seq.count('h'),\
                                     val.seq.count('A') + val.seq.count('A'),\
                                     val.seq.count('G') + val.seq.count('g'),\
                                     val.seq.count('P') + val.seq.count('p')])
    
    def abs(number):
        if number > 0 or number == 0:
            pass
        else:
            number = - number
        return number
    
    rcvCalc = 0; nTaxa = len(combine)
    for key, val in countDict.items():
        rcvCalc = rcvCalc + abs(val[0] - (numA/nTaxa)) + abs(val[1] - (numB/nTaxa)) + abs(val[2] - (numI/nTaxa)) + abs(val[3] - (numL/nTaxa)) + \
                            abs(val[4] - (numF/nTaxa)) + abs(val[5] - (numN/nTaxa)) + abs(val[6] - (numS/nTaxa)) + abs(val[7] - (numC/nTaxa)) + \
                            abs(val[8] - (numH/nTaxa)) + abs(val[9] - (numU/nTaxa)) + abs(val[10] - (numG/nTaxa)) + abs(val[11]) - (numG/nTaxa) + \
                            abs(val[12] - (numP/nTaxa))

    totalRCV = float(rcvCalc)/float(len(combine)*len(combine[0]))
    return totalRCV


def Convert(input, output, filename):
    """
        File format conversion program (fasta, strict-phylip, sequential-phylip, relaxed-phylip and nexus).
        @parameter input - Input file format.
        @parameter output - Output file format.
        @parameter filename - Input filename.
        """
    formDict = {
        'fasta': '*.fas',
        'nexus': '*.nex',
        'phylip': '*.phy',
        'phylip-sequential': '*.phy',
        'phylip-relaxed': '*.phy'
    }

    os.chdir('..')
    
    if input == 'fasta' and output == 'nexus':
        alignment = AlignIO.read(open(filename), "fasta", alphabet=Gapped(IUPAC.protein))
        g = open(filename.split(".")[0] + '.nex', 'w')
        g.write(alignment.format("nexus")); g.close()

    else:
        try:
            handle = open(filename, 'rU'); record = list(SeqIO.parse(handle, input))
            fp = open(filename.split('.')[0] + '.' + formDict[output].split('.')[1], 'w')
            SeqIO.write(record, fp, output); fp.close(); handle.close()
        
        except:
            print "Bad Alignment\n"

    print "Final output saved in %s" %filename.split('.')[0] + '.' + formDict[output].split('.')[1]

def ConvertAll(inp_format):
    """
        File format conversion program (fasta, strict-phylip, sequential-phylip and relaxed-phylip to nexus).
        @parameter inp_format - Input file format.
        """
    
    os.chdir('Input'); files = glob.glob("*.*")
    for filename in files:
        if '.nex' not in filename:
            try:
                alignment = AlignIO.read(open(filename), inp_format, alphabet=Gapped(IUPAC.protein))
                g = open(filename.split(".")[0] + '.nex', 'w')
                g.write(alignment.format("nexus")); g.close()
            except ValueError:
                print("Error raised in importing %s file" %filename)
                continue

    os.chdir('..')


def fastEvol(combined, cutOff):
    """
        Detects fast evolving sites in final concatenated alignment
        @parameter combined - is a 3D nexus data matrix.
        @cutOff - OV cutoff values supplied by user.
        Returns - Fast evolving sites and their corresponding OV values.
        """
    if cutOff == None:
        listPos = []
    else:
        msa = MultipleSeqAlignment(NexusHandler(1).combineToRecord(combined))
        charList = []; i = 1
        while i <= len(msa[0]):
            j = 1; tempList = []
            while j <= len(msa):
                tempList.append(msa[j-1][i-1])
                j = j + 1
            charList.append(tempList)
            i = i + 1

        OVdict = dict()
        for i, val in enumerate(charList):
            val = list(''.join(val).replace('?', '')); posVal = []
            if len(set(val)) == 1:
                posVal.append(0)
            else:
                outCounter = 0; inCounter = 1
                for inval in val:
                    inCounter = 1 + outCounter
                    while inCounter < len(val):
                        posVal.append(0) if inval == val[inCounter] else posVal.append(1)
                        inCounter = inCounter + 1
            
                    outCounter = outCounter + 1

            k = (math.pow(len(val), 2) - len(val))/2
            OVdict['Position_%s' %i] = (sum(posVal)/k)
        
        listPos = [[x, val] for x, val in OVdict.items() if val > cutOff]
    return listPos


def binAll(rcvRange, entropyRange, combined, RCVdict, entropyDict, gcDict, gcRange):
    """
        Creat RCV, GC and Entropy dataset for user defines bin range.
        
        @parameter rcvRange - Range of RCV values passed by user to create RCV bin.
        @parameter entropyRange - Range of Entropy values passed by user to create Entropy bin.
        @parameter gcRange - Range of percetange GC content passed by user to create GC bin.
        @parameter RCVdict - Dictionary element that conctains RCV values for all the input alignment files.
        @parameter entropyDict - Dictionary element that conatins entropy values for all the input alignment files.
        @parameter gcDict - Dictionary element that contains percentage GC content for all the input alignment files.
        @parameter combined - is a 3D nexus data matrix.
        
        Returns - list of values with alignment IDs and their corresponding RCV, GC content and Entropy values in user defined bin range.
        """
    lineListRcv = []; lineListEntropy = []; lineListGC = []
    if rcvRange != None:
        rStart = float(rcvRange.split('-')[0]); rEnd = float(rcvRange.split('-')[1])
        for key, val in RCVdict.items():
            try:
                sink = entropyDict[key]
            except KeyError:
                entropyDict[key] = (['NA'])
            try:
                sink = gcDict[key]
            except KeyError:
                gcDict[key] = (['NA'])
            if float(val.lstrip('[ ').rstrip(' ]')) >= rStart and float(val.lstrip('[ ').rstrip(' ]')) <= rEnd:
                try:
                    lineListRcv.append("BIN_RCV %s = %s-%s [RCV Score = %s] [Entropy = %s] [GC Content (in percentage) = %s]" %(key, combined.charsets[key][0], combined.charsets[key][-1], val, entropyDict[key], gcDict[key]))
                except KeyError:
                    continue

    if entropyRange != None:
        rStart = float(entropyRange.split('-')[0]); rEnd = float(entropyRange.split('-')[1])
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
        rStart = float(gcRange.split('-')[0]); rEnd = float(gcRange.split('-')[1])
        for key, val in gcDict.items():
            try:
                sink = entropyDict[key]
            except KeyError:
                entropyDict[key] = (['NA'])
            try:
                sink = RCVdict[key]
            except KeyError:
                RCVdict[key] = (['NA'])
            if val >= rStart and val <= rEnd:
                try:
                    lineListGC.append("BIN_GC %s = %s-%s [RCV Score = %s] [Entropy = %s] [GC Content (in percentage) = %s]" %(key, combined.charsets[key][0]+1, combined.charsets[key][-1]+1, RCVdict[key], entropyDict[key], val))
                except KeyError:
                    continue

    return [lineListRcv, lineListEntropy, lineListGC]


def GCcontent(combined):
    """
        Calculates GC content in alignment
        @parameter combined - is a 3D nexus data matrix.
        Returns - Dictionary with GC values
        """
    GCdict = dict()
    msa = MultipleSeqAlignment(NexusHandler(1).combineToRecord(combined))
    for key, val in combined.charsets.items():
        if key != 'RNA_Stem' and key != 'RNA_Loop' and key != "'RNA_Stem'" and key != "'RNA_Loop'":
            gcCount = 0
            try:
                msaGene = msa[:, combined.charsets[key][0]:combined.charsets[key][-1]]
                for inval in msaGene:
                    gcCount = gcCount + inval.seq.count('G') + inval.seq.count('g') + inval.seq.count('C') + inval.seq.count('c')
                GCdict[key] = (float(gcCount)/(len(msaGene)*len(msaGene[1]))*100)
            except KeyError:
                continue

    gcHist(GCdict)
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
    k = (len(N)-1) * percent; f = math.floor(k); c = math.ceil(k)
    if f == c:
        return key(N[int(k)])
    d0 = key(N[int(f)]) * (c-k); d1 = key(N[int(c)]) * (k-f)
    return d0+d1


def allPerentile(dictRCV, dictEntropy, dictGC, combined, which):
    """
        Finds the values grouped in percentile from input dictionary data.
        @parameter RCVdict - is RCV dictionary element.
        @parameter entropyDict - is Entropy dictionary element.
        @parameter gcDict - is GC content dictionary element.
        @parameter combined - is a 3D nexus data matrix
        @parameter which - takes dictionary type as input (It could be GC, RCV or Entropy dictionary)
        """
    RCVdict = dictRCV; entropyDict = dictEntropy; gcDict = dictGC
    supDict = RCVdict if which == 'rcv' else entropyDict if which == 'ent' else gcDict
    supList = []
    for key, val in supDict.items():
        supList.append(val) if isinstance(val, str) == False and hasattr(val, '__iter__') == False else None

    supPlist0to25 = []; supPlist25to75 = []; supPlist75to100 = []
    supList.sort()
    sup25 = percentile(supList, 0.25); sup75 = percentile(supList, 0.75)
    v0to25 = [x for x in supList if x >= 0 and x < sup25]
    v25to75 = [x for x in supList if x >= sup25 and x < sup75]
    v75to100 = [x for x in supList if x >= sup75]

    for key, val in supDict.items():
        try: sink = RCVdict[key]
        except KeyError: RCVdict.update({key: ['NA']})
        try: sink = entropyDict[key]
        except KeyError: entropyDict.update({key: ['NA']})
        try: sink = gcDict[key]
        except KeyError: gcDict.update({key: ['NA']})

        if val in v0to25:
            try:
                supPlist0to25.append("BIN %s = %s-%s [RCV Score = %s] [Entropy = %s] [GC Content (in percentage) = %s]"\
                                 %(key, combined.charsets[key][0]+1, combined.charsets[key][-1]+1, RCVdict[key], entropyDict[key], gcDict[key]))
            except KeyError:
                continue
        elif val in v25to75:
            try:
                supPlist25to75.append("BIN %s = %s-%s [RCV Score = %s] [Entropy = %s] [GC Content (in percentage) = %s]"\
                                  %(key, combined.charsets[key][0]+1, combined.charsets[key][-1]+1, RCVdict[key], entropyDict[key], gcDict[key]))
            except KeyError:
                continue
        else:
            try:
                supPlist75to100.append("BIN %s = %s-%s [RCV Score = %s] [Entropy = %s] [GC Content (in percentage) = %s]"\
                                   %(key, combined.charsets[key][0]+1, combined.charsets[key][-1]+1, RCVdict[key], entropyDict[key], gcDict[key]))
            except KeyError:
                continue
    
    supPdict = dict()
    supPdict['0_to_25'] = (supPlist0to25); supPdict['25_to_75'] = (supPlist25to75); supPdict['75_to_100'] = (supPlist75to100)
    return supPdict


def binPercent(RCVdict, entropyDict, gcDict, combined, calRCVvalue, runShannon, runGC):
    """
        returns the list of dictionaries that has charset data grouped under percentile bins.
        @parameter RCVdict - is RCV dictionary element.
        @parameter entropyDict - is Entropy dictionary element.
        @parameter gcDict - is GC content dictionary element.
        @parameter combined - is a 3D nexus data matrix
        """
    for key, val in RCVdict.items():
        RCVdict[key] = (float(val.lstrip('[ ').rstrip(' ]')))

    rcvPdict = allPerentile(RCVdict, entropyDict, gcDict, combined, which='rcv') if calRCVvalue == True else dict()
    entPdict = allPerentile(RCVdict, entropyDict, gcDict, combined, which='ent') if runShannon == True else dict()
    gcPdict = allPerentile(RCVdict, entropyDict, gcDict, combined, which='gc') if runGC == True else dict()
    return [rcvPdict, entPdict, gcPdict]


def removePerBin(filename):
    """
        Creates a dictionary element that contains gene name and their corresponding percentile bins.
        """
    RCVbinList = []; ENTbinList = []; GCbinList = []
    fr = open(filename[0], 'r'); data = fr.readlines()
    
    Flag = False
    for lines in data:
        Flag = False if '[Entropy Bin]' in lines else Flag
        RCVbinList.append(lines.rstrip('\n').lstrip('\tBIN ')) if Flag == True else None
        Flag = True if '[RCV Bin]' in lines else Flag
    
    for lines in data:
        Flag = False if '[GC Bin]' in lines else Flag
        ENTbinList.append(lines.rstrip('\n').lstrip('\tBIN ')) if Flag == True else None
        Flag = True if '[Entropy Bin]' in lines else Flag

    for lines in data:
        Flag = False if '[User percentile GC Bin]' or 'end' in lines else Flag
        GCbinList.append(lines.rstrip('\n').lstrip('\tBIN ')) if Flag == True else None
        Flag = True if '[GC Bin]' in lines else Flag

    binList = [RCVbinList, ENTbinList, GCbinList]
    for listval in binList:
        listval.remove('[25th to 75th percentile] '); listval.remove('[75th to 100 percentile] '); listval.remove('[0 to 25th percentile] ')
    
    for i, listval in binList:
        for j, inval in enumerate(listval):
            if inval != '':
                listval[j] = inval.split(' ')[0]
        binList[i] = listval
    
    binDict = [{}, {}, {}]; binKey = ['RCV', 'ENT', 'GC']
    for j, outVal in enumerate(binList):
        nameList = [[], [], []]; pos = [i for i, val in enumerate(outVal) if val == '']; counter = -1
        for i, val in enumerate(outVal):
            if i in pos:
                counter = counter + 1
                continue
            else:
                nameList[counter].append(val)
        
        binDict[j][binKey[j]] = (nameList)
    
    newBinDict = dict()
    for i, outVal in enumerate(binDict):
        for key, val in outVal.items():
            for j, inval in enumerate(val):
                val[j] = {'25-75': inval} if j == 0 else {'75-100': inval} if j == 1 else {'0-25': inval} if j == 2 else None
            outVal[key] = (val)
            for inkey, item in outVal.items():
                newBinDict[inkey] = {}
                for i, dictVal in enumerate(item):
                    newBinDict[inkey].update(dictVal)
    return newBinDict



def gcUserBin(part, gcDict):
    """
        Create Bin data for user defined GC partition range.
        """
    gcList = []
    for key, val in gcDict.items():
        gcList.append(val)
    
    gcList.sort(); npart = 100/part; counter = 1; myDict = dict()
    
    rangeList = range(0, 100, part)
    rangeList.append(100)
    
    for i, ranges in enumerate(rangeList):
        try:
            rangeList[i] = str(ranges) + '-' + str(rangeList[i+1])
        except IndexError:
            rangeList[i] = str(ranges) + '-' + '100'

    if rangeList[-1] == '100-100':
        rangeList = rangeList[:-1]

    for ranges in rangeList:
        myDict["Percentile[" + ranges + "]"] = list()
        for val in gcList:
            if val >= percentile(gcList, float(ranges.split('-')[0])/100) and val < percentile(gcList, float(ranges.split('-')[1])/100):
                myDict["Percentile[" + ranges + "]"].append(val)
    myDict["Percentile[" + rangeList[-1] + "]"].append(max(gcList))


    retDict = dict()

    for key, val in myDict.items():
        for i, inval in enumerate(val):
            for inkey, gcVal in gcDict.items():
                if inval == gcVal:
                    val[i] = inkey
        myDict[key] = (val)
    return myDict


def remUsrBin(combined, part):
    """Creates a list of gene keys for which the GC value lies in the user defined range"""
    gcDict = GCcontent(combined); remList = []
    for key, val in gcDict.items():
        remList.append(key) if val >= int(part.split('-')[0]) and val <= int(part.split('-')[1]) else None
    return remList


######################################################################################################

"""Sets of functions required for calculating average, variance and standard deviation from list of data."""

def average(s): return sum(s) * 1.0 / len(s)

def variance(avg, s): return map(lambda x: (x - avg)**2, s)

def stDev(variance): return math.sqrt(average(variance))

######################################################################################################

def gcHist(gcDict):
    """
        Histogram plot of GC values
        @parameter gcDict -  Dictionary of GC values
        """
    l = []
    for key, val in gcDict.items():
        l.append(int(val))
    gcArray = array('l', l)
    avg = average(l); var = variance(avg, l); sigma = stDev(var); num_bins = 20
    n, bins, patches = plt.hist(gcArray, num_bins, normed=1, facecolor='green', alpha=0.5)
    y = mlab.normpdf(bins, float(sum(l))/len(l), sigma)
    plt.plot(bins, y, 'r--')
    plt.xlabel('GC Values'); plt.ylabel('(Percentage Gene Count)/100')
    plt.title(r'Histogram of GC Content: #mean=%s,  #sigma=%s \n\n' %(float(sum(l))/len(l), sigma))
    plt.subplots_adjust(left=0.15)
    plt.savefig('GCplot.png')


def setUpdate(sets, positions):
    for key in sets:
        for i, val in enumerate(positions):
            x1=[x-1 for x in sets[key] if x > val-i]
            x2=[x for x in sets[key] if x < val-i]
            sets[key] = (x2+x1)
    
    return sets











