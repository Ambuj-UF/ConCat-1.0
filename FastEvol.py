import argparse
import textwrap
from Bio.Nexus import Nexus
from Bio import SeqIO
from Functions import fastEvol, Convert
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



args = parser.parse_args()

def main():
    file = [args.i]
    nexi =  [(fname, Nexus.Nexus(fname)) for fname in file]
    msaObject = MultipleSeqAlignment(NexusHandler(1).combineToRecord(nexi[0][1]))
    posMatrix = []

    if args.auto == True:
        with open('Fast_Evolving_Sites.txt', 'r') as fr:
            data = fr.readlines()
            for lines in data:
                posMatrix.append(int(lines.rstrip('\n'))

    else:
        data = Functions.fastEvol(nex[1][0], args.OV)
        for val in data:
            posMatrix.append(int(val[0].split('_')[1]))

    cycles = 0

    while cycles < len(posMatrix)-1:
        if cycles == 0:
            Position = posMatrix[1]
            termPos =  Position - len(msaObject[1])
            varFirst = Position-1
            msaObject = msaObject[:, :varFirst] + msaObject[:, termPos:]
    
        else:
            Position = posMatrix[cycles+1]
            termPos = (Position-cycles) - len(msaObject[1])
            varFirst = Position-cycles-1
            msaObject = msaObject[:, :varFirst] + msaObject[:, termPos:]
        cycles = cycles + 1

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


                                 
                                 
                                 
if __name__ == "__main__":
    main()













