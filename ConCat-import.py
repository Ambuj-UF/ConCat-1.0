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
import subprocess
import platform
import argparse
import textwrap
import warnings
from src.fetchSeq import cdsImport, mrnaImport
from src.aligner import cdsAlign, mrnaAlign
from Bio import AlignIO
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


parser.add_argument('-impcds', type=str, default=None,
                    help='Takes gene name via file for CDS import')

parser.add_argument('-impmrna', type=str, default=None,
                    help='Takes gene name via file for mRNA import')

parser.add_argument('-group', type=str, default=None,
                    help='Takes group name to extract sequence data')

parser.add_argument('-pkg', type=str, default='muscle',
                    choices=['muscle', 'mafft'], required=True,
                    help='User defined program selection')

parser.add_argument('-args', type=str,
                    help='Arguments to run MAFFT. EXAMPLE: "--retree 2 --maxiterate 10"')

parser.add_argument('-sep', action='store_true', default=False,
                    help='Include if you want to run different models for different alignment files')

parser.add_argument('-argf', type=str,
                    help='Takes argument file as input')


argmnts = parser.parse_args()


if argmnts.pkg == 'mafft' and argmnts.sep == False and not argmnts.args:
    parser.error('-CMND argument is required in "mafft" mode.')

if argmnts.sep == True and not argmnts.argf:
    parser.error('-sep argument is required in "-argf" mode.')

if argmnts.impcds == True and argmnts.impmrna == True:
    parser.error('-impcds argument and -impcds cannot be used together')

if argmnts.impcds == True and not argmnts.group:
    parser.error('-group argument is required in "-impcds" mode.')

if argmnts.impmrna == True and not argmnts.group:
    parser.error('-group argument is required in "-impmrna" mode.')


def warnfxn():
    warnings.warn("deprecated", DeprecationWarning)


with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    warnfxn()



def main():
    if argmnts.impcds != None:
        genes = open(argmnts.impcds, 'r').readlines()
        for geneName in genes:
            warnings.filterwarnings("ignore")
            cdsImport(geneName.rstrip('\n'), argmnts.group)
            cdsAlign(geneName.rstrip('\n') + ".fas")

    elif argmnts.impmrna != None:
        genes = open(argmnts.impmrna, 'r').readlines()
        for geneName in genes:
            warnings.filterwarnings("ignore")
            mrnaImport(geneName.rstrip('\n'), argmnts.group)
            mrnaAlign(geneName.rstrip('\n') + ".fas")

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
            if argmnts.sep == True:
                with open(argmnts.argf) as f:
                    data = f.readlines()
                    dataDict = dict()
                    for lines in data:
                        dataDict[lines.split(' = ')[0]] = (lines.split(' = ')[1])
                            
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





