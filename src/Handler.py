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
import shutil
import subprocess
import math
import collections
import os.path

from Bio.Alphabet.IUPAC import IUPACAmbiguousDNA
from collections import Counter
from sets import Set
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from Bio import SeqIO
from Bio import AlignIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio.Nexus import Nexus



class myDict(dict):
    
    def __init__(self):
        self = dict()
    
    def add(self, key, value):
        self[key] = value


class BaseHandle:
    
    def __init__(self, file_format):
        self.file_format = file_format

    ################################################################################
    #                               ID Import                                      #
    ################################################################################


    def fileOpenConc(self):
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

    ################################################################################
    # Quick Record Import. This program creates a list of records for all files    #
    ################################################################################

    
    def fileOpenID(self):
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


    #####################################################################################################
    # This function imports all the accession IDs linked with the species name in alignment matrix.     #
    #####################################################################################################

    
    def storeId(self):
        idDict = self.fileOpenID()
        print idDict
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



    #############################################################################################################################
    #           This function is meant to check the spelling errors in the species name in alignment files                      #
    #############################################################################################################################


    def fuzyName(self):
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
                        print "Found " + i + " in gene " + keyList[counter-1] + " but gene " + keyList[n] + " (has " + j+")\n"
            
                n = n + 1
            counter = counter + 1


    def selectRNAgene(self):
        listRNA = glob.glob('*.*')
        
        loopCount = 0
        while loopCount < len(listRNA):
            print loopCount + 1, '-->', listRNA[loopCount]
            loopCount = loopCount + 1
        
        for name in listRNA:
            print "Would you like to add ", name.split('.')[0], " gene for RNA processing? \n Enter Y/N \n"
            usr_inp = raw_input('\n')
            if usr_inp == 'y' or usr_inp == 'Y' or usr_inp == 'Yes' or usr_inp == 'yes':
                newListRNA.append(name)
                print name.split('.')[0], " gene has been added to the list \n"
    
            elif usr_inp == 'n' or usr_inp == 'N' or usr_inp == 'No' or usr_inp == 'no':
                continue
    
        for filename in newListRNA:
            shutil.copy(filename, 'RNAuserData')
    
    

    ######################################################################################################################
    # This function creates the file with RNA secondry strcuture and its total energy value for all the gene and species #
    ######################################################################################################################



    def RNAfoldPred(self):
        extList = ["*.fasta", "*.nexus", "*.phy", "*.phy", "*.phy"]
        fileList = glob.glob(extList[self.file_format-1])
        fp = open('ResRNA.txt', 'w')
        geneListOut = []
    
        for filename in fileList:
            geneName = filename.split(".")
            geneListOut.append(geneName[0])
            fp.write("Begin; \n")
            fp.write("[" + geneName[0] + "]\n")
            fp.write(subprocess.check_output("RNAfold < %s" % filename, shell = True))
            fp.write("End;\n\n")
    
        fp.write("Close;")
        fp.close()
        dictWidGeneName = rnaExtract(geneListOut)
        rnaMatrix(dictWidGeneName)
            

    def rnaExtract(self, geneName):
        listData = []
        data=[]
        flag=False
        outFlag = 'True'
        dictWidGeneName = {}
        with open('ResRNA.txt','r') as f:
            while outFlag == 'True':
                fp = open('partRNA.txt', 'w')
                for line in f:
                    if line.startswith('Begin;'):
                        flag = True
                    elif line.startswith('Close;'):
                        outFlag = 'False'
                
                    if flag:
                        fp.write(line)
                
                    if line.strip().endswith('End;'):
                        flag=False
                        break;
            
                fp.close()
            
                d = open('partRNA.txt', 'r')
                data = d.read()
                d.close()
                dictWidSpKey = {}
                x = 0
                data1 = data.split()
                while x < len(data1) - 3:
                    if '>' in data1[x + 2]:
                        dictWidSpKey[data1[x + 2].strip('\n')] = (data1[x + 5].split())
                
                    x = x + 4
            
                listData.append(dictWidSpKey)
                if outFlag == 'False':
                    break;
    
        counter = 0
        newGeneList = []
        for i in geneName:
            if i not in newGeneList:
                newGeneList.append(i)

        while counter < len(newGeneList):
            dictWidGeneName[newGeneList[counter]] = (listData[counter])
            counter = counter + 1

        return dictWidGeneName
        
        
    ################################################################################
    #                       Creates Matrix of RNA energy values                    #
    ################################################################################
        

    def rnaMatrix(self, dictWidGeneName):
        newdict = {}
        subkeyList = []
        for keys in dictWidGeneName:
            for subkey in dictWidGeneName[keys]:
                subkeyList.append(subkey)
    
        newSubList = []
        for i in subkeyList:
            if i not in newSubList:
                newSubList.append(i)
    
        print '\t' + '\t'.join(newSubList)
        geneKeys = dictWidGeneName.keys()
        for u in geneKeys:
            def f(v):
                try:
                    return str(dictWidGeneName[u][v])
                except KeyError:
                    return "NA"
            print u + '\t' + '\t'.join(f(v) for v in newSubList)
    


    ################################################################################
    #                       Concatenation Program                                  #
    ################################################################################


    def concatenate(self, records1, records2):
        numRecords1 = len(records1)
        numRecords2 = len(records2)

        n1 = 1
        n2 = 1
        i1 = 1
        i2 = 1

        list1 = [""]
        list2 = [""]
        item_list1 = [""]
        item_list2 = [""]

        while n1<=numRecords1:
            list1.append(records1[n1-1].id)
            n1 = n1 + 1

        while n2<=numRecords2:
            list2.append(records2[n2-1].id)
            n2 = n2 + 1

        for item in list1:
            if item in list2:
                continue
            else:
                item_list1.append(item)

        for item in list2:
            if item in list1:
                continue
            else:
                item_list2.append(item)

        lenghtOFdiff1 = len(item_list1)

        lenghtOFdiff2 = len(item_list2)


        #This section helps in substituting "?" character for the missing species in alignment files

        if lenghtOFdiff1 >= 2:
            while i1 < lenghtOFdiff1:
                dummy = SeqRecord(Seq("?"*len(records2[1])), id = item_list1[i1], description="")
                records2.append(dummy)
                i1 = i1 + 1


        if lenghtOFdiff2 >= 2:
            while i2 < lenghtOFdiff2:
                dummy = SeqRecord(Seq("?"*len(records1[1])), id = item_list2[i2], description="")
                records1.append(dummy)
                i2 = i2 + 1

        msa1 = MultipleSeqAlignment(records1)
        msa2 = MultipleSeqAlignment(records2)

        #Sorting and Concatenation

        msa1.sort()
        msa2.sort()

        complete = msa1 + msa2

        return(complete)


    ################################################################################
    #                       Produces Output result file                            #
    ################################################################################


    def alignOutput(self, combine):
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


    ################################################################################
    #                               Taxon Removal                                  #
    ################################################################################


    def taxaRemove(self, combine):
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


    def geneList(self, file_format):
        initList = ['*.fasta', '*.nex', '*.phy', '*.phy', '*.phy']
        print "The available files are \n"
        print glob.glob(initList[file_format-1])
    
        list_files = glob.glob(initList[file_format-1])
        fileCount = len(list_files)
    
        loopCount = 0
    
        print "Select from files listed below \n"
    
        while loopCount < fileCount:
            print loopCount + 1, "->", list_files[loopCount], '\n'
            loopCount = loopCount + 1
    
        usr_inp = input('\n')
    
        return list_files[usr_inp-1]


    def taxonList(self, combine):
        idlist = []
        for data in combine:
            idlist.append(data.id)
    
        loopCount = 0
        print "Select from taxa listed below \n"
    
        while loopCount < len(idlist):
            print loopCount + 1, "-->", idlist[loopCount], '\n'
            loopCount = loopCount + 1
    
        usr_inp = input('\n')
    
        return idlist[usr_inp-1]



    ################################################################################
    #   Remove positions that has gap as well as missing data for all taxa         #
    ################################################################################


    def autoGapRemove(self, combine):
        
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


    ###################################################################################
    #                                  NEXML Output                                   #
    ###################################################################################


    def nexML(self, filename):
        fp = open('Results.xml', 'w')
        handleXML = open(filename, 'rU')
        recordsXML = list(SeqIO.parse(handleXML, "nexus"))
        SeqIO.write(recordsXML, fp, "seqxml")
        fp.close()
        handleXML.close()

    ###################################################################################################
    #   RY-coding program: It replaces A & G to R and C & T to Y at the third codon position          #
    ###################################################################################################

    def RYcoding(self, file, position, msaObject):
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
                newSeqData = ReplaceThird(self, seqData, position)
                newSeqList.append(newSeqData)


        counter = 0
        while counter < len(newSeqList):
            data.append(SeqRecord(Seq(newSeqList[counter], generic_dna), id = records[counter].id, name = records[counter].name, description = records[counter].description))
            counter = counter + 1

        newmsa = MultipleSeqAlignment(data)

        return newmsa




#############################################################################################################################
#                   taxanomyClass is designed to add or remove taxanomic classes from the sequence ID's                     #
#############################################################################################################################


class taxanomyClass:
    def __init__(self, taxDict, combined):
        self.taxDict = taxDict
        self.combined = combined
    
    def addTaxanomy(self):
        combined = self.combined
        taxDict = self.taxDict
        for key, data in combined.matrix.items():
            taxID = [taxDict[inkey] for inkey in taxDict if key.split('_')[0] + '_' + key.split('_')[1]  == inkey]
            combined.matrix[key + '_' + taxID[0][0]] = combined.matrix.pop(key)
        for i, labels in enumerate(combined.taxlabels):
            taxID = [taxDict[inkey] for inkey in taxDict if labels.split('_')[0] + '_' + labels.split('_')[1]  == inkey]
            combined.taxlabels[i] = labels + '_' + taxID[0][0]

        return combined

    def remTaxanomy(self):
        combined = self.combined
        taxDict = self.taxDict
        for key, data in combined.matrix.items():
            combined.matrix[key.rstrip('_' + key.split('_')[len(key.split('_'))-1])] = combined.matrix.pop(key)
        combined.taxsets.clear()
        return combined


#############################################################################################################################
#                                       Nexus file handling operations                                                      #
#############################################################################################################################


class NexusHandler:
    
    def __init__(self, filename):
        self.filename = filename

    
    ################################################################################
    #                       Creates Consensus RNA structure                        #
    ################################################################################


    def RNAfoldConsensus(self):
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
            fp.write("[ %s ]\n" % name.split('.')[0])
            fp.write(subprocess.check_output("RNAalifold < %s" % name, shell = True))
            fp.write('\n\n')
            
        fp.close()


    def fileOpenConcNex(self):
        fileList = glob.glob("*.nex")
        recordList = []
        for filename in fileList:
            handle = open(filename, "rU")
            record = list(SeqIO.parse(handle, "nexus"))
            recordList.append(record)
            
        return recordList

    


    ################################################################################
    #                       Writes RNA data to Nexus                               #
    ################################################################################


    def RNAtoNexus(self, combined, RNAstrucData):
        
        rnaDict = {}

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
            
            posStem = [x.start() for x in re.finditer('\(', rnaData)]
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
                        print "Error in %s ConCat Block Taxon name defined in RNA_Struc variable. Taxa not found in alignment. Please check the spelling \n" % key
            
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
                        print "Error in %s ConCat Block Taxon name defined in RNA_Struc variable. Taxa not found in alignment. Please check the spelling \n" % key


                posLoop = [x.start() for x in re.finditer('\.', rnaData)]
                rnaList.append([x + startPos for x in posLoop])
            
                posStem = [x.start() for x in re.finditer('\(', rnaData)]
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


	###########################################################################
	#			Updates RNA Structure Matrix								  #	
	###########################################################################


    def rnaSetUpdate(self, positions, procRNAdict):
        for key in procRNAdict:
            for i, val in enumerate(positions):
                x1=[x-1 for x in procRNAdict[key] if x > val-i]
                x2=[x for x in procRNAdict[key] if x < val-i]
                procRNAdict[key] = (x2+x1)

        return procRNAdict
    

    ############################################################################################
    #    Calculates entropy of alignment. It helps in identifying the alignment quality        #
    ############################################################################################

    def ProtEntropy(self, msa):
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
                        'a': ['D', 'E'],
                        'b': ['R', 'K'],
                        'i': ['I', 'V'],
                        'l': ['L', 'M'],
                        'f': ['F', 'W', 'Y'],
                        'n': ['N', 'Q'],
                        's': ['S', 'T']
                                             
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
                                             
        dna = set("ATGC-N?")
        def checkDNA(seq, alphabet=dna):
            "Checks if string is a DNA sequence"
            leftover = set(seq.upper()) - alphabet
            return not leftover

        entropyDict = {}
        msaObject = MultipleSeqAlignment(self.combineToRecord(combined))
        for i, mkeys in enumerate(combined.charsets):

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
        return entropyDict



    def combineToRecord(self, combined):
        matrixNew = combined.matrix
        listID = combined.taxlabels
        records = []
        
        for ids in listID:
            records.append(matrixNew[ids])
            
        for i, sequence in enumerate(records):
            newrecord = SeqRecord(records[i], id = listID[i])
            records[i] = newrecord

        return records

    ####################################################################################################
    #               This function edits the character sets after gap removal                           #
    ####################################################################################################

    def charEdit(self, combined, positions):
        listGap=[]
        n = 0
        while n < len(combined.charsets):
            list1 = combined.charsets[combined.charsets.keys()[n]]
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
            combined.charsets[combined.charsets.keys()[p]] = (listGap[p])
            p = p + 1

        return combined

    def msaToMatrix(self, msa, combined):
        for i, values in enumerate(msa):
            combined.matrix[values.id] = msa[i].seq
        
        return combined
    

    ##############################################################################
    #                       Entropy dictionary Update                            #
    ##############################################################################


    def entropyDictUpdate(self, entropyDict, positions):
        listGap = []
        n = 0
        while n < len(entropyDict):
            list1 = entropyDict[entropyDict.keys()[n]]
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
            entropyDict[entropyDict.keys()[p]] = (listGap[p])
            p = p + 1

        for key, val in entropyDict.items():
            entropyDict[key+"_Bad_Alignment"] = entropyDict.pop(key)

        return entropyDict


    ##############################################################################
    #                       Checks if the sequence is missing data               #
    ##############################################################################


    def validate(self, seq, alphabet = set("?")):
        leftover = set(seq.upper()) - alphabet
        return not leftover


    ##############################################################################
    #                  Cretes a list of missing data for TaxaLabels              #
    ##############################################################################


    def missingScan(self, combined):
        MSA = MultipleSeqAlignment(self.combineToRecord(combined))

        missingList = {}
        for i,key in enumerate(combined.charsets):
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
        
        for key, val in missingList.items():
            missingList["Missing_"+key] = missingList.pop(key)
    
        combined.taxsets.update(missingList)
        
        return combined



    def managePipes(self):
        os.chdir("Input")
        file_list = glob.glob('*.nex')

        save_path = 'ProcInput'
        
        idDataDict = dict()

        for filename in file_list:
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
        
        return idDataDict

    ##############################################################################
    #                       Main Program                                         #
    ##############################################################################

    def runRNAOperation(self, procRNAdict, combined):
        Flag = True
        for key, val in procRNAdict.items():
            if val[-1] -1 > len(combined.matrix[combined.matrix.keys()[0]]):
                print "Please check the RNA structure manual entry in alignment files \n Looks like the length of RNA structure is greater than the length of alignment. It usually happens when you put wrong RNA structure starting position which causes final RNA structure position in the respective alignment file to exceed the total length of alignment \nSkipping RNA structure matrix construction\n"
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
        
        
        print "Checking for missing taxa in alignment files \n"
        combined = self.missingScan(combined)
        
        if runShanon == True:
            print "Checking alignment quality [Shannons Entropy] \n"
            entropyGenes = self.entropyCal(combined)
        else:
            pass
        
        if runRNA == True:
            print "Constructing RNA structure matrix \n"
            try:
                self.RNAfoldConsensus()
                procRNAdict = self.RNAtoNexus(combined, RNAstrucData)
                try:
                    combined = self.runRNAOperation(procRNAdict, combined)
                except:
                    print "Cant find RNA data \n"
            except OSError:
                print "RNAfold not found \n Skipping this step"
    
    
        if runShanon == True:
            combined.charsets.update(entropyGenes)
        else:
            pass
        
        if usr_inpT == 1:
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
            
        else:
            pass
        
        return combined
        

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
                combined = self.missingScan(combined)
            else:
                raise ValueError("% not found in the alignment. Check for spelling mistakes in taxa editing input file" % delTaxa)
    
        if usr_inpT == 1:
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
    
        else:
            pass
    
        positions = Nexus.Nexus.gaponly(combined, include_missing = True)
        
        if runShanon == True:
            print "Checking alignment quality [Shannons Entropy] \n"
            entropyGenes = self.entropyCal(combined)
            entropyGenes = self.entropyDictUpdate(entropyGenes, positions)
        else:
            pass
    
        combined = self.charEdit(combined, positions)
    
        finalAlignment = MultipleSeqAlignment(self.combineToRecord(combined))
        align = BaseHandle(2).autoGapRemove(finalAlignment)
    
        combined = self.msaToMatrix(align, combined)

        if runRNA == True:
            print "Constructing RNA structure matrix \n"
            try:
                self.RNAfoldConsensus()
                procRNAdict = self.RNAtoNexus(combined, RNAstrucData)
                procRNAdict = self.rnaSetUpdate(positions, procRNAdict)
                try:
                    combined = self.runRNAOperation(procRNAdict, combined)
                except:
                    print "Cant find RNA data \n"
            except OSError:
                print "RNAfold not found \n Skipping this step"
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
    
        return [combined, remTaxDict]
    

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
            print "Removing taxon \n"
            listExclude = []
            for lines in includeTax:
                listExclude.append(lines.rstrip('\n'))
            
            outWithTaxa = self.withTaxEdit(listExclude,
                                           combined,
                                           usr_inpT,
                                           RNAstrucData,
                                           runRNA,
                                           runShanon
                                           )
            return outWithTaxa

        elif excludeTax.isatty() == False:
            print "Removing taxon absent in include taxon file supplied by user \n"
            listInclude = []
            for lines in excludeTax:
                listInclude.append(lines.rstrip('\n'))
            
            listExclude = [taxa for taxa in combined.taxlabels if taxa not in listInclude]
            self.withTaxEdit(listExclude, combined, usr_inpT, RNAstrucData)
            outWithTaxa = self.withTaxEdit(listExclude, combined, usr_inpT, RNAstrucData)
            return outWithTaxa

            
        else:
            combined = self.withoutTaxEdit(combined,
                                           usr_inpT,
                                           RNAstrucData,
                                           runRNA,
                                           runShanon
                                           )
            remTaxDict = dict()
            return [combined, remTaxDict]


