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





import argparse
import textwrap
from Bio.Nexus import Nexus
from Bio import SeqIO
import sys
from src.Functions import fastEvol, Convert
from Bio.AlignIO import MultipleSeqAlignment
from src.Handler import NexusHandler


parser = argparse.ArgumentParser(prog='FastEvol',
                                 version= 'ConCat-1.0',
                                 formatter_class=argparse.RawDescriptionHelpFormatter,
                                 description=textwrap.dedent('''\
    ----------------------------------------------------------------------------------------------------------
    \t\t\t\t Welcome to ConCat-1.0 alignment handling program
    \t\t\t FastEvol removes all the fast evolving sites from the input alignment.
    \t\t\t
    \t\t\t\tWritten by Ambuj Kumar, University of Florida
    \t\t      ConCat sequence alignment concatenation program was designed at Kimball-Braun
    \t\t\tlab group to handle large amount of sequence dataset for phylogeny study.
    
    ----------------------------------------------------------------------------------------------------------
    
    '''))

parser.add_argument('-fevol', action='store_true', default=False,
                    help='remove fast evolving sites')

parser.add_argument('-rembin', action='store_true', default=False,
                    help='remove sections from alignment through bin selection')

parser.add_argument('-i', type=str, required=True,
                    help='Enter input alignnmnet file')

parser.add_argument('-fin', type=str, required=True, default = 'nexus',
                    help='Enter input alignnmnet file format')

parser.add_argument('-o', type=str, required=True,
                    help='Enter output alignnmnet file')

parser.add_argument('-fout', type=str, required=True, default = 'nexus',
                    help='Enter output alignnmnet file format')

parser.add_argument('-auto', action='store_true', default=False,
                    help='Include if you want to run remove the list of fast evoling sites in Fast_Evolving_Sites.txt file')

parser.add_argument('-OV', type=float, default=None,
                    help='Enter the OV cutoff value for detecting fast evolving sites')


args = parser.parse_args()

def main():
    
    if args.fevol == True:
        
        """
            Removes fast evolving sites from the concatenated alignment and updates all the charset data.
            """
        
        file = [args.i]
        nexi =  [(fname, Nexus.Nexus(fname)) for fname in file]
        msaObject = MultipleSeqAlignment(NexusHandler(1).combineToRecord(nexi[0][1]))
        posMatrix = []

        if args.auto == True:
            with open('Fast_Evolving_Sites', 'r') as fr:
                data = fr.readlines()
                for lines in data:
                    posMatrix.append(int(lines.rstrip('\n')))

        else:
            data = fastEvol(nexi[0][1], args.OV)
            if data != []:
                for val in data:
                    posMatrix.append(int(val[0].split('_')[1]))
            else:
                sys.exit('Zero fast evolvong site found. Program Terminated \n')

        cycles = 0
        posMatrix.sort()
        msanewObject = msaObject[:, :posMatrix[0]]
        for i, val in enumerate(posMatrix):
            if i > 0:
                if posMatrix[i] == posMatrix[-1]:
                    msanewObject = msanewObject + msaObject[:, val:len(msaObject[1])]
                    break
                else:
                    msanewObject = msanewObject + msaObject[:, val:posMatrix[i+1]]
    
        msaObject = msanewObject
        if args.fout == 'nexus':
            combined = nexi[0][1]
            combined = NexusHandler(1).msaToMatrix(msaObject, combined)

            def setUpdate(sets, positions):
                for key in sets:
                    for i, val in enumerate(positions):
                        x1=[x-1 for x in sets[key] if x > val-i]
                        x2=[x for x in sets[key] if x < val-i]
                        sets[key] = (x2+x1)
                    
                return sets

            combined.charsets = setUpdate(combined.charsets, posMatrix)
            combined.write_nexus_data(filename=open(args.o, 'w'))
                                 
        else:
            with open(args.o, 'w') as fp:
                SeqIO.write(msaObject, fp, args.fout)

    elif args.rembin == True:
        file = [args.i]
        nexi =  [(fname, Nexus.Nexus(fname)) for fname in file]
        msaObject = MultipleSeqAlignment(NexusHandler(1).combineToRecord(nexi[0][1]))



if __name__ == "__main__":
    main()













