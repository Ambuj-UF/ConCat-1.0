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





import os
import re
import glob
import time
import shutil
import subprocess
import math
import collections

from sets import Set
from collections import Counter

from Utils.Bio import SeqIO
from Utils.Bio import AlignIO
from Utils.Bio.Seq import Seq
from Utils.Bio.Nexus import Nexus
from Utils.Bio.SeqRecord import SeqRecord
from Utils.Bio.Alphabet import generic_dna
from Utils.Bio.Align import MultipleSeqAlignment
from Utils.Bio.Alphabet.IUPAC import IUPACAmbiguousDNA



class myDict(dict):
    
    def __init__(self):
        self = dict()
    
    def add(self, key, value):
        self[key] = value


class BaseHandle:
    
    """
       BaseHandle class deals with handeling and storing alignment IDs and detecting 
       spelling errors in the taxon name. Additionaly it contains functions for RY
       coding and generating NexML file output.
        
    """
    
    def __init__(self, file_format):
        self.file_format = file_format

    ################################################################################
    #                               ID Import                                      #
    ################################################################################


    def fileOpenConc(self):
        
        """
            This function is used in opening input alignment files for concatenation.
            It can handle file formats other than Nexus.
            
        """
        
        file_format = self.file_format
        
        extList = ["*.fas", "*.phy", "*.phy", "*.phy"]
        typeList = ["fasta", "phylip", "phylip-sequential", "phylip-relaxed"]
        fileList = glob.glob(extList[file_format-1])
    
        recordList = []
        for filename in fileList:
            handle = open(filename, "rU")
            record = list(SeqIO.parse(handle, typeList[file_format - 1]))
            recordList.append(record)
    
        return recordList


    
    def fileOpenID(self):
        
        """
           Quick Record Import. This program creates a list of records for all files.
           It can handle Nexus record objects.
           
        """
        
        file_format = self.file_format
        
        extList = ["*.fas", "*.nex", "*.phy", "*.phy", "*.phy"]
        typeList = ["fasta", "nexus", "phylip", "phylip-sequential", "phylip-relaxed"]
        
        fileList = glob.glob(extList[file_format - 1])
    
        dict = {}
    
        for filename in fileList:
            handle = open(filename, "rU")
            idList = []
            for record in SeqIO.parse(handle, typeList[file_format - 1]):
                idList.append(record.id)
            gene = filename.split(".")[0]
            dict[gene] = idList
            handle.close()
        
        return dict


    
    def storeId(self):
        
        """
           This function imports all the accession IDs linked with the species name in alignment matrix.
           
        """
        
        idDict = self.fileOpenID()
        print("%s" %idDict)
        with open("AccessionID.txt", 'w') as fp:
            fp.write("Gene  Species  Database    Accession ID \n")
    
            for key in idDict:
                idList = idDict[key]
            
            for data in idList:
                data = data.split('|')
                n = 1
                while n < len(data) - 2:
                    fp.write(key + "\t" + data[0] + "\t" + data[n] + "\t" + data[n + 1] + "\n")
                    n = n + 2
        
                fp.write("\n")
            
    
        message = "Accession numbers are stored in file [AccessionID.txt] \n"
    
        return message



    def fuzyName(self):
        
        """
           This function is meant to check spelling errors in the taxon name in alignment files
           
        """
        
        idDict = self.fileOpenID()
        idList = []
        keyList = []
    
        for key, val in idDict.items():
            try:
                idList.append([x.split('|')[0] for x in val])
            except KeyError:
                idList.append(val)
        
            keyList.append(key)
    
        counter = 1
    
        for value in idList:
            n = counter
            while n < len(idList):
                for i, j in zip(value, idList[n]):
                    if 1.0 > float([x == y for (x, y) in zip(i, j)].count(True))/len(j) > 0.8:
                        print("WARNING             |    %s|    Found %s in gene %s but gene %s (has %s)" %(time.strftime("%c"), i, keyList[counter-1], keyList[n], j))
                        #print("Found " + i + " in gene " + keyList[counter-1] + " but gene " + keyList[n] + " (has " + j+")\n")
            
                n = n + 1
            counter = counter + 1



    def alignOutput(self, combine):
        
        """
           alignOutput creates an output file in user defined file format
           
           @parameter combine - concatenated alignment matrix
           
        """
        
        output_format = self.file_format
        if output_format == 1:
            filecompname = "Result1.fasta"
            file_Write = open(filecompname, 'w')
            SeqIO.write(combine, file_Write, "fasta")
            file_Write.close()
            #This section is for cleaning any unknown description tag from the final fasta file
        
            fin = open("Result1.fasta", "r")
            fout = open("Result.fasta", "w+")
        
            input_data = fin.readlines()
        
            for line in input_data:
                if "<unknown description>" in line:
                    line = line.replace("<unknown description>", "")
                fout.write(line)
            fin.close()
            fout.close()
        
        elif output_format == 2:
            file_Write = open("Result.phy", 'w')
            SeqIO.write(combine, file_Write, "phylip")
            file_Write.close()
    
        elif output_format == 3:
            file_Write = open("Result.phy", 'w')
            SeqIO.write(combine, file_Write, "phylip-sequential")
            file_Write.close()
    
        elif output_format == 4:
            file_Write = open("Result.phy", 'w')
            SeqIO.write(combine, file_Write, "phylip-relaxed")
            file_Write.close()
    
        else:
            sys.exit("You have enetered wrong value \n Program Terminated...")



    def taxaRemove(self, combine):
        
        """
           This function is used to remove taxa from the alignment matrix.
           
           @parameter combine - final concatenated alignment matrix.
           
           Returns - alignment matrix without specific taxa.
           
        """
        
        rem_taxa = self.taxonList(combine)
    
        n=1
    
        while n <= len(combine):
            if combine[n-1].id == rem_taxa:
                align1 = combine[ :n-1, :]
                align2 = combine[n:len(combine), :]
                align1.extend(align2)
            n = n + 1
    
        combine = align1
        return combine




    def autoGapRemove(self, combine):
        
        """
           Remove positions that has gap as well as missing data for all taxa
           
           @parameter combine - final concatenated alignment matrix
           
        """
        
        similarityCount = {}
        n=0
    
        while n < len(combine[0]):
            similarityCount[n] = 0
            n = n + 1
    
        posMatrix = [""]
    
        i = 1
    
        while i <= len(combine[0]):
            j = 1
            while j <= len(combine):
                character = combine[j-1][i-1]
                if character == "_" or character == "?" or character == "-":
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

        return combine



    def nexML(self, filename):
        
        """
            Produces concatenated alignment file in NexML format.
        """
        
        fp = open('Results.xml', 'w')
        handleXML = open(filename, 'rU')
        recordsXML = list(SeqIO.parse(handleXML, "nexus"))
        SeqIO.write(recordsXML, fp, "seqxml")
        fp.close()
        handleXML.close()


    def RYcoding(self, file, position, msaObject):
        
        """
           RY-coding program: It replaces A & G to R and C & T to Y either user defined positions
           or at all the positions. It depends on user selection.
           
           @parameter file -
           @parameter position - user defined position to perform RY coding in alignment matrix
           @parameter msaObject - Input multiple sequence alignment matri data
           
           Return - Multiple sequence alignment object with RY coding
        """
        
        def ReplaceThird(self, string, position):
            for i in range(position, len(string), 3):
                if string[i] == 'A' or string[i] == 'G':
                    string = string[:i-1] + "R" + string[i:]
                elif string[i] == 'C' or string[i] == 'T':
                    string = string[:i-1] + "Y" + string[i:]
            return string
                    
        def ReplaceAll(self, string):
            for i in range(1, len(string), 1):
                if string[i] == 'A' or string[i] == 'G':
                    string = string[:i-1] + "R" + string[i:]
                elif string[i] == 'C' or string[i] == 'T':
                    string = string[:i-1] + "Y" + string[i:]
            return string

        handle = open("Results.nex", "rU")
        records = list(SeqIO.parse(handle, "nexus"))
        handle.close()
        
        msa = msaObject

        seqlist = []
        idlist = []
        data = []
        x = 0
        while x < len(msa):
            sequence = ""
            y=0
            idlist.append(msa[x].id)
            while y < len(msa[1]):
                sequence = sequence + msa[x][y]
                y = y + 1
            seqlist.append(sequence)
            x = x + 1
            
        
        newSeqList = []
        
        if position == 'all':
            for seqData in seqlist:
                newSeqData = ReplaceAll(self, seqData)
                newSeqList.append(newSeqData)
        else:
            for seqData in seqlist:
                newSeqData = ReplaceThird(self, seqData, int(position))
                newSeqList.append(newSeqData)


        counter = 0
        while counter < len(newSeqList):
            data.append(SeqRecord(Seq(newSeqList[counter], generic_dna),\
                                  id = records[counter].id, name = records[counter].name,\
                                  description = records[counter].description))
                                  
            counter = counter + 1

        newmsa = MultipleSeqAlignment(data)

        return newmsa




#############################################################################################################################
#                   taxanomyClass is designed to add or remove taxanomic classes from the sequence ID's                     #
#############################################################################################################################


class taxanomyClass:
    
    """
       taxanomyClass is designed to add or remove taxanomic names from the sequence ID's
       It takes taxonomy data from Taxonomy.csv file available in the ConCat package 
       home directory.
        
    """
    
    def __init__(self, taxDict, combined):
        self.taxDict = taxDict
        self.combined = combined
    
    def addTaxanomy(self):
        
        """
           Add taxonomic data to the current sequence Ids
        """
        
        combined = self.combined
        taxDict = self.taxDict
        
        for key, data in combined.matrix.items():
            if key.split('_')[0] + '_' + key.split('_')[1] in list(taxDict.keys()):
                taxID = taxDict[key.split('_')[0] + '_' + key.split('_')[1]]
                combined.matrix[key + '_' + taxID[0]] = combined.matrix.pop(key)
        for i, labels in enumerate(combined.taxlabels):
            if labels.split('_')[0] + '_' + labels.split('_')[1] in list(taxDict.keys()):
                taxID = taxDict[labels.split('_')[0] + '_' + labels.split('_')[1]]
                combined.taxlabels[i] = labels + '_' + taxID[0]

        return combined

    def remTaxanomy(self):
        
        """
            Removes taxonomic data from the current sequence Ids
        """

        combined = self.combined
        taxDict = self.taxDict
        for key, data in combined.matrix.items():
            if key.count('_') >= 1:
                combined.matrix['_'.join(key.split('_')[:-1])] = combined.matrix.pop(key)
        for i, labels in enumerate(combined.taxlabels):
            if labels.count('_') >= 1:
                combined.taxlabels[i] = '_'.join(labels.split('_')[:-1])
        
        combined.taxsets.clear()

        return combined



class NexusHandler:
    
    """
        Nexus file handling operations.
        
        This Class handles overall architecture of the nexus alignment concatenation process.
        It contains functions for creating, mapping and updating RNA structure data to the
        alignment matrix. Furthermore, it handles all the sequence annotation and update,
        character set handling and update, gap removal, alignment entropy calculation, 
        and the programs for calculating missing data, storing and handling database Ids.
        
    """
    
    def __init__(self, filename):
        self.filename = filename

    


    def RNAfoldConsensus(self):
        
        """
           Creates RNA structure data from consensus alignment using RNAfold program.
           Output is stored in RNAConsensus.txt file
           
        """
        
        os.chdir("RNAdata")
        fileList = glob.glob('*.nex')

        newFileList = []
    
        for name in fileList:
            file_name = name.split('.')[0]
            newName = file_name + '.aln'
            newFileList.append(newName)
        
        recordList = self.fileOpenConcNex()
        os.chdir("..")
        n = 0
        while n < len(fileList):
            record = recordList[n]
            file_Write = open(newFileList[n], 'w')
            SeqIO.write(record, file_Write, "clustal")
            file_Write.close()
            n = n + 1
                              
        fp = open("RNAConsensus.txt", 'w')

        for name in newFileList:
            print("RNA structure     |    %s|    %s" %(time.strftime("%c"), name))
            fp.write("[ %s ]\n" % name.split('.')[0])
            fp.write(subprocess.check_output("RNAalifold < %s" % name, shell = True))
            fp.write('\n\n')
            
        fp.close()


    def fileOpenConcNex(self):
        
        """
           This functions creates a list of alignment records from the files stored in RNAdata
           directory.
           
           Returns - List of alignment records.
        """
        
        fileList = glob.glob("*.nex")
        recordList = []
        for filename in fileList:
            handle = open(filename, "rU")
            record = list(SeqIO.parse(handle, "nexus"))
            recordList.append(record)
            
        return recordList



    def RNAtoNexus(self, combined, RNAstrucData):
        
        """
           This functions creates a dictionary element of RNA structure coordinates for
           the alignment files marked for RNA structure prediction.
           
           @parameter combined - 3D nexus data matrix
           @parameter RNAstrucData - RNA data obtained from ConCat block (Supplied by user).
           
           Returns - Dictionary element with RNA structure coordinate and their corresponding
           alignment file names.
           
        """

        rnaDict = {}
        rnaList = list()

        f = open("RNAConsensus.txt", 'r')
        fdata = f.readlines()
        n = 0
        while n < len(fdata):
            rnaList = []
            
            if '|' in fdata[n]:
                geneName = fdata[n].split(' ')[1].split('|')[0]  + '.nex'
            
            else:
                geneName = fdata[n].split(' ')[1] + '.nex'

            rnaData = fdata[n+2].split(' ')[0]

            posLoop = [x.start() for x in re.finditer('\.', rnaData)]
            rnaList.append(posLoop)
            
            posStem = [x.start() for x in re.finditer('\(', rnaData)] + [x.start() for x in re.finditer('\)', rnaData)]
            rnaList.append(posStem)

            rnaDict[geneName] = (rnaList)

            n = n + 5

        f.close()
        
        def insert_loop(string, index):
            return string[:index] + '.' + string[index:]
        
        if (RNAstrucData and True) or False == True:
            
            for key, val in RNAstrucData.items():
                if ',' in val:
                    taxName = val.split(',')[0]
                    rnaDataRaw = val.split(',')[1]
                    startPos = int(val.split(',')[2])
                    if taxName in combined.taxlabels:
                        seq = combined.matrix[taxName][combined.charsets[key][0]: combined.charsets[key][-1]]
                        try:
                            gapPos = [x.start() for x in re.finditer('\-', str(seq))]
                        
                            for i, pos in enumerate(gapPos):
                                rnaDataRaw = insert_loop(rnaDataRaw, pos + i)
                        except TypeError:
                            raise ValueError("%s sequence of %s is empty or not readbale \n" % (taxName, key))
                        
                        rnaData = rnaDataRaw
                    
                    else:
                        print("Error in %s ConCat Block Taxon name defined in RNA_Struc variable. Taxa not found in alignment. Please check the spelling \n" % key)
            
                else:
                    taxName = val.split(',')[0]
                    rnaDataRaw = val.split(',')[1]
                    startPos = 0
                    if taxName in combined.taxlabels:
                        seq = combined.matrix[taxName][combined.charsets[key][0]: combined.charsets[key][-1]]
                        try:
                            gapPos = [x.start() for x in re.finditer('\-', str(seq))]
                        
                            for i, pos in enumerate(gapPos):
                                rnaDataRaw = insert_loop(rnaDataRaw, pos + i)
                    
                        except TypeError:
                            raise ValueError("%s sequence of %s is empty or not readbale \n" % (taxName, key))
                            
                        rnaData = rnaDataRaw

            
                    else:
                        print("Error in %s ConCat Block Taxon name defined in RNA_Struc variable. Taxa not found in alignment. Please check the spelling \n" % key)


                posLoop = [x.start() for x in re.finditer('\.', rnaData)]
                rnaList.append([x + startPos for x in posLoop])
            
                posStem = [x.start() for x in re.finditer('\(', rnaData)] + [x.start() for x in re.finditer('\)', rnaData)]
                rnaList.append([x + startPos for x in posStem])
            
                rnaDict[key] = (rnaList)
        
        setsDict = combined.charsets
        
        for key in rnaDict:
            for inkey in setsDict:
                if key == inkey:
                    for i, val in enumerate(rnaDict[key]):
                        for j, inval in enumerate(val):
                            rnaDict[key][i][j] = val[j] + setsDict[inkey][0]

        listLoop = []
        listStem = []

        for keys in rnaDict:
            listLoop.append(rnaDict[keys][0])
            listStem.append(rnaDict[keys][1])

        def unlistList(item):
            inList = []
            for val in item:
                for inval in val:
                    inList.append(inval)

            return inList

        listLoop = unlistList(listLoop)
        listStem = unlistList(listStem)

        procRNAdict = {'RNA_Loop': listLoop, 'RNA_Stem': listStem}
        
        for files in glob.glob("*.aln"):
            os.remove(files)

        return procRNAdict



    def rnaSetUpdate(self, positions, procRNAdict):
        
        """
           Updates RNA structure coordinates after gap/missing site removal.
        """
        
        for key in procRNAdict:
            for i, val in enumerate(positions):
                x1=[x-1 for x in procRNAdict[key] if x > val-i]
                x2=[x for x in procRNAdict[key] if x < val-i]
                procRNAdict[key] = (x2+x1)

        return procRNAdict
    


    def ProtEntropy(self, msa):
        
        """
           Calculates entropy of protein alignment. It helps in identifying the alignment quality.
            
        """
        
        totEntropy = []
        n = 1
        while n <= len(msa[0]):
            seq = msa[:, n-1:n]
            charList = [[d for d in c] for c in seq]
            for i, val in enumerate(charList):
                for d in val:
                    charList[i] = d
                                             
            alphabet = [e for e in Set(charList)]
            if '?' in alphabet:
                alphabet.remove('?')
                                             
            aminoDict = {
                        'u': ['D', 'E'],
                        'b': ['R', 'K'],
                        'i': ['I', 'V'],
                        'l': ['L', 'M'],
                        'f': ['F', 'W', 'Y'],
                        'n': ['N', 'Q'],
                        's': ['S', 'T'],
                        'c': ['C'],
                        'h': ['H'],
                        'a': ['A'],
                        'g': ['G'],
                        'p': ['P']
                                             
            }
                                             
            for i, val in enumerate(alphabet):
                for key, val in aminoDict.items():
                    if val in aminoDict[key]:
                        alphabet[i] = key
                                             
            freqList = []
            for symbol in alphabet:
                ctr = 0
            for sym in charList:
                if sym == symbol:
                    ctr = ctr + 1
            freqList.append(float(ctr) / len(charList))
                                             
            colEntropy =[]
            for values in freqList:
                colEntropy.append(values*math.log(values))
                                             
            colEnt = -sum(colEntropy)
                                             
            totEntropy.append(colEnt)
            n = n + 1
                                             
        Entropy = (sum(totEntropy)/(len(msa[0]) + 1))
        return Entropy
                                             

                                             

    def ShannonEntropy(self, msa):
        
        """
            Calculates entropy of DNA alignment. It helps in identifying the alignment quality.
            
        """
        
        totEntropy = []
        n = 1
        while n <= len(msa[0]):
            seq = msa[:, n-1:n]
            charList = [[d for d in c] for c in seq]
            for i, val in enumerate(charList):
                for d in val:
                    charList[i] = d
        
            alphabet = [e for e in Set(charList)]
            if '?' in alphabet:
                alphabet.remove('?')
                                             
            
            freqList = []
            for symbol in alphabet:
                ctr = 0
                for sym in charList:
                    if sym == symbol:
                        ctr = ctr + 1
                freqList.append(float(ctr) / len(charList))

            colEntropy =[]
            for values in freqList:
                colEntropy.append(values*math.log(values))

            colEnt = -sum(colEntropy)

            totEntropy.append(colEnt)
            n = n + 1
        
        Entropy = (sum(totEntropy)/(len(msa[0]) + 1))
        return Entropy


    def entropyCal(self, combined):
        
        """
            Calls for entropy calculation is regulated by this function.
            
            @parameter combined - 3D matrix of nexus data
            
            Returns - List of two dictionary element.
            1. entropyDict - Contains entropy data for alignments having entropy value > 0.7.
            2. entropyValDict - Contains entropy data for all alignment files.
            
        """
        
        print("Entropy Calculation |    %s|    -------Initiating Entropy Calculation-------" %(time.strftime("%c")))
                                             
        dna = set("ATGC-N?")
        
        def checkDNA(seq, alphabet=dna):
            """Checks if string is a DNA sequence"""
            
            leftover = set(seq.upper()) - alphabet
            return not leftover
        
        entropyValDict = dict()
        entropyDict = {}
        msaObject = MultipleSeqAlignment(self.combineToRecord(combined))
        totalLength = len(list(combined.charsets))
    
        for i, mkeys in enumerate(combined.charsets):
            if mkeys.count(".nex") == 1:
                if float(i)/totalLength*100 < 10:
                    print("Entropy Calculation |    %s|    %.2f percent completed             |    %s" %(time.strftime("%c"), float(i)/totalLength*100, mkeys))
                elif 10 <= float(i)/totalLength*100 < 100:
                    print("Entropy Calculation |    %s|   %.2f percent completed             |    %s" %(time.strftime("%c"), float(i)/totalLength*100, mkeys))
                else:
                    print("Entropy Calculation |    %s|  %.2f percent completed             |    %s" %(time.strftime("%c"), float(i)/totalLength*100, mkeys))
                
                start = combined.charsets[mkeys][0]
                end = combined.charsets[mkeys][-1]
                if checkDNA(msaObject[:, start:end][1].seq) == True:
                    entropyValue = self.ShannonEntropy(msaObject[:, start:end])

                else:
                    entropyValue = self.ProtEntropy(msaObject[:, start:end])
                                             
                if entropyValue > 1:
                    entropyDict[mkeys] = (combined.charsets[mkeys])
                else:
                    pass
    
                entropyValDict[mkeys] = (entropyValue)
                
        return [entropyDict, entropyValDict]



    def combineToRecord(self, combined):
        
        """
           Magic function of this program. It changes 3D nexus data matrix to SeqRecord object.
           Used when handling and editing multiple sequence alignemnt data stored in the nexus 
           matrix.
           
        """
        
        try:
            matrixNew = combined.matrix
        except AttributeError:
            raise IOError("Files not found in Input directory")
        
        listID = combined.taxlabels
        records = []
        
        for ids in listID:
            records.append(matrixNew[ids])
            
        for i, sequence in enumerate(records):
            newrecord = SeqRecord(records[i], id = listID[i])
            records[i] = newrecord
    
    
        return SeqRecord(records)



    def charEdit(self, combined, positions):
        
        """
           This function edits the character sets after gap removal.
           This function is a life saver for ConCat.
            
        """
        
        listGap=[]
        n = 0
        while n < len(combined.charsets):
            list1 = combined.charsets[list(combined.charsets.keys())[n]]
            for val in positions:
                if val in list1:
                    list1.remove(val)
    
            listGap.append(list1)
            n = n + 1


        for m,values in enumerate(positions):
            for j, inval in enumerate(listGap):
                for i, pos in enumerate(inval):
                    if pos > values - m:
                        listGap[j][i] = listGap[j][i] - 1

        p = 0
        while p < len(combined.charsets):
            combined.charsets[list(combined.charsets.keys())[p]] = (listGap[p])
            p = p + 1

        return combined

    def msaToMatrix(self, msa, combined):
        
        """
           Another magic function of ConCat. It is used for pushing edited multiple
           sequence alignment matrix into 3D nexus matrix.
           
           @parameter msa - edited multiple sequence alignment object
           @parameter combined - Nexus 3D matrix data
           
           Returns - Nexus matrix with updated sequence alignment data.
        """
        
        for i, values in enumerate(msa):
            combined.matrix[values.id] = msa[i].seq
        
        return combined
    


    def entropyDictUpdate(self, entropyDict, positions):
        
        """
           Updates entropy dictionary data after gap removal.
           
           @parameter entropyDict - Dictionary of entropy data.
           @parameter positions - coloumns with gap + missing data at all taxa positions
           
           Returns - Dictionary with updated entropy value
        """
        
        listGap = []
        n = 0
        while n < len(entropyDict):
            list1 = entropyDict[list(entropyDict.keys())[n]]
            for val in positions:
                if val in list1:
                    list1.remove(val)
        
            listGap.append(list1)
            n = n + 1
    
        for m,values in enumerate(positions):
            for j, inval in enumerate(listGap):
                for i, pos in enumerate(inval):
                    if pos > values - m:
                        listGap[j][i] = listGap[j][i] - 1
    
        p = 0
        while p < len(entropyDict):
            entropyDict[list(entropyDict.keys())[p]] = (listGap[p])
            p = p + 1

        for key, val in entropyDict.items():
            entropyDict[key+"_Bad_Alignment"] = entropyDict.pop(key)

        return entropyDict


    def validate(self, seq, alphabet = set("?")):
        """ Checks if the sequence is missing data """
        
        leftover = set(seq.upper()) - alphabet
        return not leftover



    def missingScan(self, combined):
        """Cretes a list of missing taxon data to store in Taxsets"""
        
        MSA = MultipleSeqAlignment(self.combineToRecord(combined))

        missingList = {}
        totalLength = len(list(combined.charsets))
        counter = 0
        
        for i,key in enumerate(combined.charsets):
            if key.count(".nex") == 1:
                if float(i)/totalLength*100 < 10:
                    print("Missing Scan        |    %s|    %.2f percent completed             |    Scanning missing taxa in %s" %(time.strftime("%c"), float(i)/totalLength*100, key))
                elif 10 <= float(i)/totalLength*100 < 100:
                    print("Missing Scan        |    %s|   %.2f percent completed             |    Scanning missing taxa in %s" %(time.strftime("%c"), float(i)/totalLength*100, key))
                else:
                    print("Missing Scan        |    %s|  %.2f percent completed             |    Scanning missing taxa in %s" %(time.strftime("%c"), float(i)/totalLength*100, key))
            
                try:
                    start = combined.charsets[key][0]
                    end = combined.charsets[key][-1]
        
                    n = 0
                    while n < len(MSA):
                        dataSeq = MSA[n][start:end].seq
                        if self.validate(dataSeq) == True:
                            try:
                                missingList[key].append(MSA[n].id)
                            except KeyError:
                                missingList[key] = [(MSA[n].id)]
                        n = n + 1
    
                except IndexError:
                    pass

        if missingList == {}:
            print("Running ConCat-build|    %s|    No missing taxa found" %(time.strftime("%c")))
            #print("No missing taxa found")
        
        else:
            
            print("Missing Scan        |    %s|    --------------- Processing Output ---------------" %(time.strftime("%c")))

            for key, val in missingList.items():
                missingList["Missing_"+key] = missingList.pop(key)
        
            if combined.taxsets == {}:
                combined.taxsets = dict()
                for key, val in missingList.items():
                    combined.taxsets[key] = (val)

        
            else:
                combined.taxsets.update(missingList)
            
            print("Missing Scan        |    %s|    ---------- Missing scan step completed ----------" %(time.strftime("%c")))
        
        return combined



    def managePipes(self):
        
        """
            This function extracts sequence Ids if available in the alignment file and
            stores the edited alignment file in ProcInput directory
            
            Returns - Dictionary element with gene name, taxon name and sequence Ids.
            
        """
        
        os.chdir("Input")
        file_list = glob.glob('*.nex')

        save_path = 'ProcInput'
        
        idDataDict = dict()
        
        print("Running ConCat-build|    %s|    Creating database ID-less files for Concatenation" %(time.strftime("%c")))
        
        totalLength = len(file_list)
        counter = 0

        for filename in file_list:
            if float(counter)/totalLength*100 < 10:
                print("Extract ID's        |    %s|    %.2f percent extraction completed  |    Extracting from and transferring %s file" %(time.strftime("%c"), float(counter)/totalLength*100, filename))
            elif 10 <= float(counter)/totalLength*100 < 100:
                print("Extract ID's        |    %s|   %.2f percent extraction completed  |    Extracting from and transferring %s file" %(time.strftime("%c"), float(counter)/totalLength*100, filename))
            else:
                print("Extract ID's        |    %s|  %.2f percent extraction completed  |    Extracting from and transferring %s file" %(time.strftime("%c"), float(counter)/totalLength*100, filename))
            counter = counter + 1
            completeName = os.path.join(save_path, filename.split('.')[0] + ".nex")
            fp = open(completeName, 'w')
            flist = [filename]

            nexi =  [(fname, Nexus.Nexus(fname)) for fname in flist]
            combined = Nexus.combine(nexi)
            
            try:
                combined.matrix = {k.split('|')[0]: v for k, v in combined.matrix.items()}
            except KeyError:
                pass
        
            idOneDict = dict()
            for i, IDs in enumerate(combined.taxlabels):
                try:
                    idOneDict[IDs[0: IDs.find('|')]] = (IDs[IDs.find('|') + 1: len(IDs)])
                    combined.taxlabels[i] = IDs.split('|')[0]
                except KeyError:
                    pass

            idDataDict[filename] = (idOneDict)
            combined.write_nexus_data(fp)
            fp.close()
        os.chdir("..")
        
        print("Extract ID's        |    %s|   ---------- ID extraction and file transfer done ----------" %(time.strftime("%c")))
        
        return idDataDict



    ##############################################################################
    #                       Main Program                                         #
    ##############################################################################
    
    

    def runRNAOperation(self, procRNAdict, combined):
        Flag = True
        for key, val in procRNAdict.items():
            if val[-1] -1 > len(combined.matrix[list(combined.matrix.keys())[0]]):
                print("Running ConCat-build|    %s|    Please check the RNA structure manual entry in alignment files" %(time.strftime("%c")))
                print("Running ConCat-build|    %s|    Looks like the length of RNA structure is greater than the length of alignment." %(time.strftime("%c")))
                print("Running ConCat-build|    %s|    It usually happens when you put wrong RNA structure starting position which causes final RNA structure position in the respective alignment file to exceed the total length of alignment" %(time.strftime("%c")))
                print("Running ConCat-build|    %s|    Skipping RNA structure matrix construction" %(time.strftime("%c")))
                
                #print("Please check the RNA structure manual entry in alignment files \n Looks like the length of RNA structure is greater than the length of alignment. It usually happens when you put wrong RNA structure starting position which causes final RNA structure position in the respective alignment file to exceed the total length of alignment \nSkipping RNA structure matrix construction\n")
                Flag = False
            if Flag:
                combined.charsets.update(procRNAdict)
    
        return combined

    def withoutTaxEdit(self,
                       combined,
                       usr_inpT,
                       RNAstrucData,
                       runRNA,
                       runShanon
                       ):
        
        
        #print("Checking for missing taxa in alignment files \n")
        print("Running ConCat-build|    %s|    Checking for missing taxa in alignment files" %(time.strftime("%c")))
        combined = self.missingScan(combined)
        entropyStore = dict()
        if runShanon == True:
            
            entropyDetails = self.entropyCal(combined)
            entropyGenes = entropyDetails[0]
            entropyStore = entropyDetails[1]
        
        elif runShanon == False:
            entropyStore = dict()
        else:
            pass
        
        if runRNA == True:
            print("Running ConCat-build|    %s|    Constructing RNA structure matrix" %(time.strftime("%c")))
            #print("Constructing RNA structure matrix \n")
            try:
                self.RNAfoldConsensus()
                procRNAdict = self.RNAtoNexus(combined, RNAstrucData)
                try:
                    combined = self.runRNAOperation(procRNAdict, combined)
                except:
                    print("Running ConCat-build|    %s|    Cant find RNA data" %(time.strftime("%c")))
                    #print("Cant find RNA data \n")
            except OSError:
                print("Running ConCat-build|    %s|    RNAfold not found. Skipping this step" %(time.strftime("%c")))
                #print("RNAfold not found \n Skipping this step")
    
    
        if runShanon == True:
            combined.charsets.update(entropyGenes)
        else:
            pass
        
        if usr_inpT == 1:
            
            print("Running ConCat-build|    %s|    ---------------- Updating taxsets ---------------" %(time.strftime("%c")))
            
            os.chdir('Input')
            listIDs = BaseHandle(2).fileOpenID()
            os.chdir('..')
            for idKey, val in listIDs.items():
                l = []
                for inval in val:
                    data = [val for i, val in enumerate(inval.split('|')) if i>=1]
                    l.append(inval)
                listIDs[idKey] = (l)
                listIDs['Database_IDs_' + idKey] = listIDs.pop(idKey)
                
            combined.taxsets.update(listIDs)
    
            print("Running ConCat-build|    %s|    ----------------- Taxsets updated ---------------" %(time.strftime("%c")))
    
            
        else:
            pass
        
        return [combined, entropyStore]
        

    def withTaxEdit(self,
                    listExclude,
                    combined,
                    usr_inpT,
                    RNAstrucData,
                    runRNA,
                    runShanon
                    ):
        

        remTaxList = []
        remTaxDict = dict()

        for delTaxa in listExclude:
            if delTaxa in combined.taxlabels:
                remTaxList.append(delTaxa)
                if delTaxa in combined.matrix: del combined.matrix[delTaxa]
                if delTaxa in combined.taxlabels: combined.taxlabels.remove(delTaxa)
            
            else:
                print("Running ConCat-build|    %s|    %s not found in the alignment. Check for spelling mistakes in taxa editing input file" %(time.strftime("%c"), delTaxa))
                #print("%s not found in the alignment. Check for spelling mistakes in taxa editing input file" % delTaxa)
        combined = self.missingScan(combined)
        if usr_inpT == 1:
            print("Running ConCat-build|    %s|    ---------------- Updating taxsets ---------------" %(time.strftime("%c")))
            os.chdir('Input')
            listIDs = BaseHandle(2).fileOpenID()
            os.chdir('..')
            for idKey, val in listIDs.items():
                l = []
                for inval in val:
                    if inval.split('|')[0] in combined.taxlabels:
                        data = [val for i, val in enumerate(inval.split('|')) if i>=1]
                        l.append(inval)
                listIDs[idKey] = (l)
                listIDs['Database_IDs_' + idKey] = listIDs.pop(idKey)
        
            combined.taxsets.update(listIDs)
    
            print("Running ConCat-build|    %s|    ----------------- Taxsets updated ---------------" %(time.strftime("%c")))
    
        else:
            pass
    
        positions = Nexus.Nexus.gaponly(combined, include_missing = True)
        entropyStore = dict()
        if runShanon == True:
            
            entropyDetails = self.entropyCal(combined)
            entropyGenes = entropyDetails[0]
            entropyStore = entropyDetails[1]
            entropyGenes = self.entropyDictUpdate(entropyGenes, positions)
        else:
            pass
    
        combined = self.charEdit(combined, positions)
    
        finalAlignment = MultipleSeqAlignment(self.combineToRecord(combined))
        align = BaseHandle(2).autoGapRemove(finalAlignment)
    
        combined = self.msaToMatrix(align, combined)

        if runRNA == True:
            print("Running ConCat-build|    %s|    Constructing RNA structure matrix" %(time.strftime("%c")))
            #print("Constructing RNA structure matrix \n")
            try:
                self.RNAfoldConsensus()
                procRNAdict = self.RNAtoNexus(combined, RNAstrucData)
                procRNAdict = self.rnaSetUpdate(positions, procRNAdict)
                try:
                    combined = self.runRNAOperation(procRNAdict, combined)
                except:
                    print("Running ConCat-build|    %s|    Cant find RNA data" %(time.strftime("%c")))
                    #print("Cant find RNA data \n")
            except OSError:
                print("Running ConCat-build|    %s|    RNAfold not found. Skipping this step" %(time.strftime("%c")))
                #print("RNAfold not found \n Skipping this step")
        else:
            pass

        if runShanon == True:
            combined.charsets.update(entropyGenes)
        else:
            pass

        try:
            remTaxDict["REMOVED_TAX"] = ([remTaxList])
        except:
            pass
    
        return [combined, remTaxDict, entropyStore]
    

    def NexusHandle(self,
                    combined,
                    usr_inpT,
                    RNAstrucData,
                    runRNA,
                    runShanon,
                    includeTax,
                    excludeTax
                    ):

        if includeTax.isatty() == False:
            print("Running ConCat-build|    %s|    Removing taxon absent in include taxon file supplied by user" %(time.strftime("%c")))
            #print("Removing taxon absent in include taxon file supplied by user \n")
            
            listInclude = []
            for lines in includeTax:
                listInclude.append(lines.rstrip('\n'))

            listExclude = [taxa for taxa in combined.taxlabels if taxa not in listInclude]
            
            
            outWithTaxa = self.withTaxEdit(listExclude,
                                           combined,
                                           usr_inpT,
                                           RNAstrucData,
                                           runRNA,
                                           runShanon
                                           )
            return outWithTaxa

        elif excludeTax.isatty() == False:
            
            #print("Removing taxon \n")
            print("Running ConCat-build|    %s|    Removing taxon" %(time.strftime("%c")))
            listExclude = []
            for lines in excludeTax:
                listExclude.append(lines.rstrip('\n'))

            outWithTaxa = self.withTaxEdit(listExclude,
                                           combined,
                                           usr_inpT,
                                           RNAstrucData,
                                           runRNA,
                                           runShanon
                                           )
            return outWithTaxa

            
        else:
            combinedRet = self.withoutTaxEdit(combined,
                                           usr_inpT,
                                           RNAstrucData,
                                           runRNA,
                                           runShanon
                                           )
            remTaxDict = dict()
            return [combinedRet[0], remTaxDict, combinedRet[1]]


