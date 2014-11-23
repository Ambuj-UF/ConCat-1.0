################################################################################################################
# Alignment editing program. It removes indels and stop codons and provides an option to adjust the alignment  #
# according to user defined protein or species sequence                                                        #
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
import os
import platform
import subprocess

try:
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from Bio.Seq import _translate_str
    from Bio.Data import CodonTable
    from Bio.Alphabet import generic_dna
    from Bio.Alphabet import IUPAC
except ImportError, e:
    sys.exit("Biopython not found")


def _spliter(str, num):
    '''Splits the string object'''
    return [ str[start:start+num] for start in range(0, len(str), num) ]



def _groupy(L):
    first = last = L[0]
    for n in L[1:]:
        if n - 1 == last:
            last = n
        else:
            yield first, last
            first = last = n
    yield first, last


def _translator(recordData, ign, omit, table):
    proteinSeqList = list()
    recordsFunc = recordData
    
    for i, rec in enumerate(recordsFunc):
        counter = dict()
        seqT = _translate_str(str(rec.seq), table)
        
        if ign == False:
            if "*" in seqT:
                counter['one'] = seqT.count('*')
                seqT = _translate_str(str(rec.seq[1:len(rec.seq)] + Seq("N", generic_dna)), table)
                if "*" in seqT:
                    counter['two'] = seqT.count('*')
                    seqT = _translate_str(str(rec.seq[2:len(rec.seq)] + Seq("NN", generic_dna)), table)
                    if "*" in seqT:
                        counter['three'] = seqT.count('*')
                        if omit == False:
                            if min(counter, key=counter.get) == 'one':
                                seqT = _translate_str(str(rec.seq), table)
                            elif min(counter, key=counter.get) == 'two':
                                seqT = _translate_str(str(rec.seq[1:len(rec.seq)] + Seq("N", generic_dna)), table)
                                recordsFunc[i].seq = recordsFunc[i].seq[1:len(rec.seq)] + Seq("N", generic_dna)
                            elif min(counter, key=counter.get) == 'three':
                                seqT = _translate_str(str(rec.seq[2:len(rec.seq)] + Seq("NN", generic_dna)), table)
                                recordsFunc[i].seq = recordsFunc[i].seq[2:len(rec.seq)] + Seq("NN", generic_dna)
                    
                    else:
                        seqT = _translate_str(str(rec.seq[2:len(rec.seq)] + Seq("NN", generic_dna)), table)
                        recordsFunc[i].seq = recordsFunc[i].seq[2:len(rec.seq)] + Seq("NN", generic_dna)
                else:
                    seqT = _translate_str(str(rec.seq[1:len(rec.seq)] + Seq("N", generic_dna)), table)
                    recordsFunc[i].seq = recordsFunc[i].seq[1:len(rec.seq)] + Seq("N", generic_dna)
            
            else:
                pass
        
        for j, obj in enumerate(seqT):
            if '*' in obj:
                seqT = seqT[:j] + 'Z' + seqT[j+1:]
        
        proteinSeqList.append(SeqRecord(Seq(seqT, IUPAC.protein), id=rec.id, name=rec.name, description=rec.description))
    
    
    with open('translated.fas', 'w') as fp:
        SeqIO.write(proteinSeqList, fp, 'fasta')
    
    
    return recordsFunc




def _alignP(pkg, arguments=None):
    if pkg == 'muscle':
        if 'Darwin' in platform.system():
            subprocess.call("./src/muscle/muscle -in translated.fas -out tAligned.fas", shell=True)
        else:
            subprocess.call("./src/muscle/muscleLinux -in translated.fas -out tAligned.fas", shell=True)
    
    else:
        arguments = arguments.replace('[', '').replace(']', '')
        subprocess.call("./src/mafft/mafft.bat %s translated.fas > tAligned.fas" %arguments, shell=True)


def _cleanAli(recordNuc, omit, fileName):
    handleP = open('tAligned.fas', 'rU')
    records = list(SeqIO.parse(handleP, 'fasta'))
    
    store = list()
    for i, rec in enumerate(records):
        nucData = [x.seq for x in recordNuc if x.id in rec.id]
        nucSeqData = _spliter(nucData[0], 3)
        sequence = Seq("", generic_dna); pos = 0
        for j, amino in enumerate(rec.seq):
            if amino == '-':
                sequence = sequence + Seq("---", generic_dna)
            elif amino == 'Z':
                sequence = sequence + Seq("NNN", generic_dna)
                pos = pos + 1
            else:
                try:
                    sequence = sequence + nucSeqData[pos]
                    pos = pos + 1
                except:
                    if rec.id not in store:
                        store.append(rec.id)
        
        
        records[i].seq = Seq(str(sequence), generic_dna)

    records = [x for x in records if x.id not in store]
    if store != []:
        print("Failed to align following sequences: %s" %store)
    
    if omit == False:
        with open("Input/" + fileName.split('.')[0] + ".nex", 'w') as fp:
            SeqIO.write(records, fp, "nexus")
    else:
        with open("Input/" + fileName.split('.')[0] + "_omited.nex", 'w') as fp:
            SeqIO.write(records, fp, "nexus")

    os.remove('translated.fas')
    os.remove('tAligned.fas')


def cdsAlign(inputFile, pkg='muscle', omit=True, ign=False, CT=None):
    
    codonTables = ['Ascidian Mitochondrial', 'SGC9', 'Coelenterate Mitochondrial', 'Protozoan Mitochondrial', 'Vertebrate Mitochondrial', 'Plant Plastid', 'Thraustochytrium Mitochondrial', 'Blepharisma Macronuclear', 'Mold Mitochondrial', 'Invertebrate Mitochondrial', 'Standard', 'Trematode Mitochondrial', 'Scenedesmus obliquus Mitochondrial', 'Euplotid Nuclear', 'Yeast Mitochondrial', 'Spiroplasma', 'Alternative Flatworm Mitochondrial', 'Ciliate Nuclear', 'SGC8', 'Alternative Yeast Nuclear', 'Hexamita Nuclear', 'SGC5', 'SGC4', 'SGC3', 'SGC2', 'SGC1', 'SGC0', 'Flatworm Mitochondrial', 'Dasycladacean Nuclear', 'Chlorophycean Mitochondrial', 'Mycoplasma', 'Bacterial', 'Echinoderm Mitochondrial']
    
    if CT == None:
        table = CodonTable.ambiguous_dna_by_id[1]
    elif CT != None and CT in codonTables:
        table = CodonTable.ambiguous_generic_by_name[CodonTable]
    else:
        table = CodonTable.ambiguous_generic_by_name['Standard']

    
    handle = open("Align/" + inputFile, 'rU')
    records = list(SeqIO.parse(handle, 'fasta'))
    for j, rec in enumerate(records):
        if 'TAA' in rec.seq[-3:] or 'TGA' in rec.seq[-3:] or 'TAG' in rec.seq[-3:]:
            records[j].seq = rec.seq[0:-3]

    if omit == True:
        badQuality = list()
        fdata = open("Align/" + inputFile.split('.')[0] + '.log', 'r').readlines()
        for lines in fdata:
            badQuality.append(lines.split(' ')[0])
        
        newRecords = list()
        for rec in records:
            if rec.id.split('|')[1] not in badQuality:
                newRecords.append(rec)
        
        records = newRecords

    records = _translator(records, ign, omit, table)
    _alignP(pkg)
    _cleanAli(records, omit, inputFile)


def mrnaAlign(inputFile, pkg, arguments=None):
    if pkg != 'muscle' and arguments == None:
        pkg = 'muscle'
    
    if pkg == 'muscle':
        if 'Darwin' in platform.system():
            subprocess.call("./src/muscle/muscle -in %s -out %s" %("Align/"+inputFile, "Input/"+inputFile), shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
            handle = open("Input/"+inputFile, 'rU')
            record = list(SeqIO.parse(handle, 'fasta'))
            with open("Input/" + inputFile.split('.')[0] + ".nex", 'w') as fp:
                SeqIO.write(record, fp, 'nexus')
            os.remove("Input/"+inputFile)
        else:
            subprocess.call("./src/muscle/muscleLinux -in %s -out %s" %("Align/"+inputFile, "Input/"+inputFile), shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
            handle = open("Input/"+inputFile, 'rU')
            record = list(SeqIO.parse(handle, 'fasta'))
            with open("Input/" + inputFile.split('.')[0] + ".nex", 'w') as fp:
                SeqIO.write(record, fp, 'nexus')
            os.remove("Input/"+inputFile)
    else:
        arguments = arguments.replace('[', '').replace(']', '')
        subprocess.call("./src/mafft/mafft.bat %s %s > %s" %(arguments,"Align/"+inputFile, "Input/"+inputFile), shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        handle = open("Input/"+inputFile, 'rU')
        record = list(SeqIO.parse(handle, 'fasta'))
        with open("Input/" + inputFile.split('.')[0] + ".nex", 'w') as fp:
            SeqIO.write(record, fp, 'nexus')
        os.remove("Input/"+inputFile)













