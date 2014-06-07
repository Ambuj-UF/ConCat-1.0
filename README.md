ConCat-1.0
==========
COMMAND LINE OPERATION

python ConCat.py -ftype [input filetype[defult is nexus]] -otype [output filetype[defult is nexus]] -RNA -rcv -pipe -shannon -RY -inc [filename] -exc [filename] -addT [value] -remT [integer value]

DEFAULT OPERATION

python ConCat.py

Help:

python ConCat.py -h

Requirements:

BioPython and python 2.7 and above

use -RNA option if RNAfold installed on the system

Prepare your own Taxonomy.csv file before using -remT and addT functions

use -pipe option to run analysis on datasets that contaons database Id in taxon name:

example: Homo_sapiens|gi|571026644|ref|NM_014453.3| 

ConCat Block:

ConCat provide users to define the alignment file type, define files to run RNA structure mapping and supply user defined RNA structure through ConCat block option.

ConCat block architecture:

#NEXUS

begin ConCat;
  Ali_Type = DNA; 
  RNA_Type = True; 
  RNA_Struc = (((((((...((...,6;
end;

begin data;
dimensions ntax=32 nchar=912;
format datatype=protein missing=? gap=-;
matrix

Homo_sapiens CCGAACAATTCTGCGCGAGGTAGGGAGGCCATGGCG....................


Ali_Type: This variable is used for defining the alignment type. It can either be DNA, Codons, Introns, Proteins or any other alignment type. Program creates a partition file for RaxML using information supplied in Ali_Type variable.
RNA_Type = True/None. If True then the ConCat program uses the alignment file for RNA structure prediction.
RNA_Struc = This variable takes RNA structure as input. If you plan to enter RNA structure then avoid using RNA_Type function. Input contains RNA structure followed by comma and RNA structure starting position. By default the starting position is set to 0 (which is the starting position of alignment).







