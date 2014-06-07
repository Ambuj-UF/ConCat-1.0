ConCat-1.0
==========


=======
  -h, --help            show this help message and exit
  
  -v, --version         show program's version number and exit
  
  -ftype {fasta,nexus,phylip,phylip-interleived,phylip-relaxed}
                        Enter the input file format for Concatenation. Default
                        is nexus.
                        
  -otype {fasta,nexus,phylip,phylip-interleived,phylip-relaxed}
                        Enter the output file format for Concatenation.
                        Default is nexus.
                        
  -spell                Include if you want to check spelling mistakes in
                        alignment files
                        
  -block                Include if you have ConCat block defined in alignment
                        files
                        
  -pipe                 Include if you have IDs stored in taxon names in the
                        alignment file
                        
  -RNA                  Include if you want to run RNAfold structure
                        prediction
                        
  -inc [INC]            Enter the text file with taxon names to be included in
                        the alignment file
                        
  -exc [EXC]            Enter the tet file with taxon names to be removed from
                        the alignment file
                        
  -shannon              Include if you want to run Shannons Entropy
                        Calculation for the alignments
                        
  -rcv                  Include if you want to run RCV Calculation for the
                        alignments
                        
  -addT {Class,Family,Order,Phylum,Kingdom}
                        Enter the txaonomy classes to add in the final
                        alignment taxa name. Sperate Multiple values using
                        comma (-addT Class,Family,Order)
                        
  -remT REMT            Enter the numer of txaonomy classes to remove from the
                        end of taxon namees
                        
  -RY [RY]              Enter the text file with alignment file names and
                        positions for RY coding
                        
  -convert              Converts fasta and phylip alignment files to nexus
                        alignment

=======


DEFAULT OPERATION

python ConCat.py

Help:

python ConCat.py -h

Requirements:

BioPython and python 2.7 and above

Handling different file formats: Importance of ConCat lies in its ability to store and generate rich annotations data in nexus output file. It allows user to define several features of alignemnt file through ConCat block as discussed above. Thus, it is important for ConCat to use nexus files as input for better functionality. To make this easier for users, ConCat -convert function scans for the files with user defined file format in Input Directory and converts it into nexus file format. 

Just type:

python ConCat.py -convert -ftype filetype [fasta, phylip, relaxed-phylip, interlieved-phylip]

Now proceed to analysis. Don't bother removing old fileformat input files from Input directory. Run -convert function and your good to go.


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





Some commonly used operation:

Run simultaneous fasta to nexus input file format conversion and analysis with RCV calculation, Alignment Entropy calculation and spellin check function set as True for a set of alignment files that has sequence IDs stored in taxon name:

python ConCat.py -CA -rcv -shannon -spell -pipe








