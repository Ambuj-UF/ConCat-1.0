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
import glob
import shutil
import subprocess
import platform
import argparse
import textwrap
import warnings
from time import sleep
from src.fetchSeq import cdsImport, mrnaImport, oneGeneCdsImport, oneGeneMrnaImport
from src.aligner import cdsAlign, mrnaAlign
from Bio import AlignIO, SeqIO
from StringIO import StringIO
from Bio.Alphabet import IUPAC, Gapped


parser = argparse.ArgumentParser(prog='ConCat-Align',
                                 version= 'ConCat-1.0',
                                 formatter_class=argparse.RawDescriptionHelpFormatter,
                                 description=textwrap.dedent('''\
    ----------------------------------------------------------------------------------------------------------
    \t\t\t\t Welcome to ConCat-1.0 alignment handling program
    \t\t\t ConCat-Align conducts sequence alignment using Muscle/Mafft 
    \t\t\t sequence alignment programs.
    \t\t\tWritten by Ambuj Kumar, University of Florida
    \t      ConCat sequence alignment concatenation program was designed at Kimball-Braun
    \t\tlab group to handle large amount of sequence dataset for phylogeny study.
    
    ----------------------------------------------------------------------------------------------------------
    
    '''))


parser.add_argument('-cds', type=str, default=None,
                    help='Takes gene name via file for CDS import')

parser.add_argument('-mrna', type=str, default=None,
                    help='Takes gene name via file for mRNA import')

parser.add_argument('-gcds', type=str, default=None,
                    help='Takes gene name via file for CDS import')

parser.add_argument('-gmrna', type=str, default=None,
                    help='Takes gene name via file for CDS import')

parser.add_argument('-orgn', type=str, default=None,
                    help='Takes group name to extract sequence data')

parser.add_argument('-pkg', type=str, default='muscle',
                    choices=['muscle', 'mafft'], required=True,
                    help='User defined program selection')

parser.add_argument('-args', type=str,
                    help='Arguments to run MAFFT. EXAMPLE: "--retree 2 --maxiterate 10"')

parser.add_argument('-argf', type=str, default=None,
                    help='Takes argument file as input')


argmnts = parser.parse_args()


if argmnts.pkg == 'mafft' and not argmnts.args:
    parser.error('-args argument is required in "mafft" mode.')

if argmnts.cds == True and argmnts.mrna == True:
    parser.error('-cds argument and -cds cannot be used together')

if argmnts.cds == True and not argmnts.orgn:
    parser.error('-orgn argument is required in "-cds" mode.')

if argmnts.mrna == True and not argmnts.orgn:
    parser.error('-orgn argument is required in "-mrna" mode.')


def warnfxn():
    warnings.warn("deprecated", DeprecationWarning)


with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    warnfxn()


def remDuplicate():
    files = glob.glob('Align/*.fas')
    for file in files:
        handle = open(file, 'rU')
        records = list(SeqIO.parse(handle, 'fasta'))
        store = list()
        print("Removing duplicate sequences from file %s" %file)
        for x in records:
            for y in [z for z in records if z.id != x.id]:
                if x.id.split('|')[0] == y.id.split('|')[0]:
                    if x.seq.count('-') > y.seq.count('-'):
                        store.append(x.id)
                    elif x.seq.count('-') <= y.seq.count('-'):
                        store.append(y.id)
    
        newRec = [x for x in records if x.id not in store]
        with open(file, 'w') as fp:
            SeqIO.write(newRec, fp, 'fasta')


def main():
    if argmnts.cds != None:
        genes = [x for x in open(argmnts.cds, 'r').readlines() if x != '' and x != '\n']
        for geneName in genes:
            try:
                warnings.filterwarnings("ignore")
                cdsImport(geneName.rstrip('\n'), argmnts.orgn)
                remDuplicate()
                cdsAlign(geneName.rstrip('\n') + ".fas")
            except:
                print("Failed to import %s sequence. HTTP Connection error!!\n" %geneName)
                print("Resuming parser in 5 seconds\n")
                sleep(5)
                try:
                    warnings.filterwarnings("ignore")
                    cdsImport(geneName.rstrip('\n'), argmnts.orgn)
                    remDuplicate()
                    cdsAlign(geneName.rstrip('\n') + ".fas")
                except:
                    print("Failed to import %s sequence. Skipping %s" %geneName)
                    continue


    elif argmnts.mrna != None:
        genes = [x for x in open(argmnts.mrna, 'r').readlines() if x != '' and x != '\n']
        for geneName in genes:
            try:
                warnings.filterwarnings("ignore")
                mrnaImport(geneName.rstrip('\n'), argmnts.orgn)
                remDuplicate()
                mrnaAlign(geneName.rstrip('\n') + ".fas")
            except:
                print("Failed to import %s sequence. HTTP Connection error!!\n" %geneName)
                print("Resuming parser in 5 seconds\n")
                sleep(5)
                try:
                    warnings.filterwarnings("ignore")
                    mrnaImport(geneName.rstrip('\n'), argmnts.orgn)
                    remDuplicate()
                    mrnaAlign(geneName.rstrip('\n') + ".fas")
                except:
                    print("Failed to import %s sequence. Skipping %s" %geneName)
                    continue


    elif argmnts.gcds != None:
        genes = [x for x in open(argmnts.gcds, 'r').readlines() if x != '' and x != '\n']
        geneRecord = list()
        for geneName in genes:
            try:
                warnings.filterwarnings("ignore")
                geneRecord.append(oneGeneCdsImport(geneName.rstrip('\n'), argmnts.orgn))
            except:
                print("Failed to import %s sequence. HTTP Connection error!!\n" %geneName)
                print("Resuming parser in 5 seconds\n")
                sleep(5)
                try:
                    warnings.filterwarnings("ignore")
                    geneRecord.append(oneGeneCdsImport(geneName.rstrip('\n'), argmnts.orgn))
                except:
                    print("Failed to import %s sequence. Skipping %s" %geneName)
                    continue

            os.remove("Align/" geneName + ".log")
    
        with open("Align/" + argmnts.orgn + ".fas", 'w') as fp:
            SeqIO.write(geneRecord, fp, "fasta")
        
        cdsAlign(argmnts.orgn + ".fas")
        shutil.move("Input/" + argmnts.orgn + ".nex", argmnts.orgn + ".nex")
        
            
            
    elif argmnts.gmrna != None:
        genes = [x for x in open(argmnts.gmrna, 'r').readlines() if x != '' and x != '\n']
        geneRecord = list()
        for geneName in genes:
            try:
                warnings.filterwarnings("ignore")
                geneRecord.append(oneGeneMrnaImport(geneName.rstrip('\n'), argmnts.orgn))
            except:
                print("Failed to import %s sequence. HTTP Connection error!!\n" %geneName)
                print("Resuming parser in 5 seconds\n")
                sleep(5)
                try:
                    warnings.filterwarnings("ignore")
                    geneRecord.append(oneGeneCdsImport(geneName.rstrip('\n'), argmnts.orgn))
                except:
                    print("Failed to import %s sequence. Skipping %s" %geneName)
                    continue
    
            os.remove("Align/" geneName + ".log")

        with open("Align/" + argmnts.orgn + ".fas", 'w') as fp:
            SeqIO.write(geneRecord, fp, "fasta")
                
        cdsAlign(argmnts.orgn + ".fas")
        shutil.move("Input/" + argmnts.orgn + ".nex", argmnts.orgn + ".nex")


    else:
        files = glob.glob("Align/*.*")
        files.remove('Align/README.txt')

        if argmnts.pkg == 'muscle':
            for filename in files:
                fname = filename.replace('Align/', '').split('.')[0] + '.fas'
                if 'Darwin' in platform.system():
                    subprocess.call("./src/muscle/muscle -in %s -out Output/%s -verbose -refine" %(filename, fname), shell=True)
                else:
                    subprocess.call("./src/muscle/muscleLinux -in %s -out Output/%s -verbose -refine" %(filename, fname), shell=True)


        elif argmnts.pkg == 'mafft':
            if argmnts.argf != None:
                with open(argmnts.argf) as f:
                    data = f.readlines()
                    dataDict = dict()
                    for lines in data:
                        dataDict[lines.split('=')[0].lstrip(' ').rstrip(' ')] = (lines.split('=')[1].lstrip(' ').rstrip(' '))
                            
                for filename in files:
                    fname = filename.replace('Align/', '').split('.')[0] + '.fas'
                    try:
                        subprocess.call("./src/mafft/mafft.bat %s %s > Output/%s" %(dataDict[filename.replace('Align/', '')], filename, fname), shell=True)
                    except:
                        print "Error in argument passed for %s in %s file" %(filename.replace('Align/', ''), argmnts.argf)
                        continue


            else:
                for filename in files:
                    fname = filename.replace('Align/', '').split('.')[0] + '.fas'
                    arguments = argmnts.args.replace('[', '').replace(']', '')
                    subprocess.call("./src/mafft/mafft.bat %s %s > Output/%s" %(arguments, filename, fname), shell=True)

        os.chdir('Output')
        files = glob.glob('*.fas')

        for filename in files:
            try:
                alignment = AlignIO.read(open(filename), "fasta", alphabet=Gapped(IUPAC.protein))
                g = open('../Input/' + filename.split(".")[0] + '.nex', 'w')
                g.write(alignment.format("nexus"))
                g.close()
        
            except ValueError:
                continue
    

        os.chdir('..')




if __name__ == "__main__":
    main()





