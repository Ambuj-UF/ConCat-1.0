################################################################################################################
# Basic alignment Handling                                                                                     #
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

try:
    import Bio
except ImportError, e:
    sys.exit('BioPython not found on your system. Program Terminated')

try:
    import argparse
except ImportError, e:
    sys.exit('Requires python 2.7 and above to run ConCat-1.0. Program Terminated')


import textwrap
from src.Nex import richNexusCall
from src.Functions import Convert, ConvertAll



parser = argparse.ArgumentParser(prog='ConCat',
                                 version= 'ConCat-1.0',
                                 formatter_class=argparse.RawDescriptionHelpFormatter,
                                 description=textwrap.dedent('''\
                                                            ----------------------------------------------------------------------------------------------------------
                                                            \t\t\t\t Welcome to ConCat-1.0 
                                                            \t\t\t A sequence alignment Concatenation Program
                                                            \t\t\tWritten by Ambuj Kumar, University of Florida
                                                            \t      ConCat sequence alignment concatenation program was designed at Kimball-Braun
                                                            \t\tlab group to handle large amount of sequence dataset for phylogeny study.
                                                            
                                                            ----------------------------------------------------------------------------------------------------------
                                                            
                                                            COMMAND LINE OPERATION
                                                            
                                                            python ConCat.py -ftype [input filetype[defult is nexus]]
                                                                    -otype [output filetype[defult is nexus]] -RNA -rcv -pipe
                                                                    -shannon -RY -inc [filename] -exc [filename] -addT [value]
                                                                    -remT [integer value]
                                                                    
                                                            ----------------------------------------------------------------------------------------------------------
                                                            
                                                            '''))


parser.add_argument('-ftype', type=str, default='nexus',
                    choices=['fasta', 'nexus', 'phylip', 'phylip-interleived', 'phylip-relaxed'],
                    help='Enter the input file format for Concatenation. Default is nexus.')

parser.add_argument('-otype', type=str, default='nexus',
                    choices=['fasta', 'nexus', 'phylip', 'phylip-interleived', 'phylip-relaxed'],
                    help='Enter the output file format for Concatenation. Default is nexus.')

parser.add_argument('-RNA', action='store_true', default=False,
                    help='Include if you want to run RNAfold structure prediction')

parser.add_argument('-inc', nargs='?', type=argparse.FileType('r'), default=sys.stdin,
                    help='Enter the text file with taxon names to be included in the alignment file')

parser.add_argument('-exc', nargs='?', type=argparse.FileType('r'), default=sys.stdin,
                    help='Enter the tet file with taxon names to be removed from the alignment file')

parser.add_argument('-shannon', action='store_true', default=False,
                    help='Include if you want to run Shannons Entropy Calculation for the alignments')

parser.add_argument('-rcv', action='store_true', default=False,
                    help='Include if you want to run RCV Calculation for the alignments')

parser.add_argument('-addT', type=str, default=None,
                    choices=['Class', 'Family', 'Order', 'Phylum', 'Kingdom'],
                    help='Enter the txaonomy classes to add in the final alignment taxa name. Sperate Multiple values using comma (-addT Class,Family,Order)')

parser.add_argument('-remT', type=str, default=None,
                    choices=['Class', 'Family', 'Order', 'Phylum', 'Kingdom'],
                    help='Enter the numer of txaonomy classes to remove from the end of taxon namees')

parser.add_argument('-pipe', action='store_true', default=False,
                    help='Include if you have IDs stored in taxon names in the alignment file')

parser.add_argument('-RY', nargs='?', type=argparse.FileType('r'), default=sys.stdin,
                    help='Enter the text file with alignment file names and positions for RY coding')



args = parser.parse_args()

def main():
    if args.ftype == 'nexus' and args.otype == 'nexus':
        richNexusCall(args.RNA,
                      args.inc,
                      args.exc,
                      args.shannon,
                      args.rcv,
                      args.addT,
                      args.remT,
                      args.pipe,
                      args.RY
                      )
    
    else:
        ConvertAll(input)
        richNexusCall(args.RNA,
                      args.inc,
                      args.exc,
                      args.shannon,
                      args.rcv,
                      args.addT,
                      args.remT,
                      args.pipe,
                      args.RY
                      )
                      
        Convert('nexus', args.otype, 'Combined.nex')




if __name__ == "__main__":
    main()