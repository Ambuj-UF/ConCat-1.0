#ConCat-1.0
###Read Manual for clear instructions


####Quick instructions

###Few basic arguments

```

  -h, --help    ->                 show this help message and exit
  
  -v, --version    ->              show program's version number and exit
      
  -ftype {fasta,nexus,phylip,phylip-interleived,phylip-relaxed}
                      ->           Enter the input file format for Concatenation. Default
                                is nexus.
                        
  -otype {fasta,nexus,phylip,phylip-interleived,phylip-relaxed}
                      ->          Enter the output file format for Concatenation.
                                Default is nexus.
                        
  -spell              ->              Include if you want to check spelling mistakes in
                                alignment files
                        
  -block              ->              Include if you have ConCat block defined in alignment
                                files
                        
  -pipe               ->              Include if you have IDs stored in taxon names in the
                                alignment file
                        
  -RNA                ->              Include if you want to run RNAfold structure
                                prediction
                        
  -inc [INC]          ->              Enter the text file with taxon names to be included in
                                the alignment file
                        
  -exc [EXC]          ->              Enter the tet file with taxon names to be removed from
                                the alignment file
                        
  -shannon            ->              Include if you want to run Shannons Entropy
                                Calculation for the alignments
                        
  -rcv                ->              Include if you want to run RCV Calculation for the
                                alignments
                        
  -addT {Class,Family,Order,Phylum,Kingdom}
                      ->              Enter the txaonomy classes to add in the final
                                alignment taxa name. Sperate Multiple values using
                                comma (-addT Class,Family,Order)
                          
  -remT REMT          ->              Enter the numer of txaonomy classes to remove from the
                                end of taxon namees
                        
  -RY [RY]            ->              Enter the text file with alignment file names and
                                positions for RY coding
                        
  -convert            ->              Converts fasta and phylip alignment files to nexus
                                alignment

```

##Requirements:

BioPython and python 2.7 and above

Handling different file formats: Importance of ConCat lies in its ability to store and generate rich annotation data in nexus output file. It allows user to define several features of alignemnt file through ConCat block. Thus, it is important for ConCat to use nexus files as input for extended functionality. To make this easier, ConCat -convert function scans for the files with user defined file format in 'Input' Directory and converts it into nexus file format. 


##Instructions
####Running ConCat with default setup.

In default setup ConCat assumes that the alignment files are placed in "Input" directory and assumes nexus as input and output alignment format.
```
python ConCat.py
```

####Help command
```
python ConCat.py -h
```


###Quick start

1. Put all your alignment files in "Input directory"
2. run "python ConCat.py" on command line

`
Input-Output format supported : [fasta, phylip, relaxed-phylip, interlieved-phylip]
`

File formats can be supplied through arguments -ftype and -otype

```
python ConCat.py -ftype fasta -otype phylip-relaxed
```

Use -RY argument to perform RY coding 

-RY argument requires and input file that contains filename and RY coding position separated by comma. For RY coding at 3rd position it typically looks like

File1.nex, 3

File2.nex, 3

File3.nex, 3


###File format conversion:

```
python ConCat.py -convert -ftype filetype [fasta, phylip, relaxed-phylip, interlieved-phylip]
```


use -pipe argument to run analysis on datasets that contains database Id in taxon name:
`
example: Homo_sapiens|gi|571026644|ref|NM_014453.3| 
`

```
python ConCat.py -pipe
```

###RNA structure mapping

use -RNA option if RNAfold is installed on the system. -RNA option can be initiated with -block argument.

-block argument allows ConCat to read ConCat block present in input nexus alignment files.

```
python ConCat.py -ftype fasta -otype nexus -block -RNA
```
**This function will perform RNA mapping for the alignment files that has RNA_Type = True assigned in its ConCat block.** 

####ConCat Block Architecture:

ConCat allows users to define the alignment file type, define files to run RNA structure mapping and supply user defined RNA structure through ConCat block option.

ConCat block architecture:

```
begin ConCat;

  Ali_Type = DNA; 

  RNA_Type = True; 

  RNA_Struc = Species_Name,(((((((...((...,6;

end;

begin data;

dimensions ntax=32 nchar=912;

format datatype=protein missing=? gap=-;

matrix

Homo_sapiens CCGAACAATTCTGCGCGAGGTAGGGAGGCCATGGCG....................
```

Ali_Type: This variable is used for defining the alignment type. It can either be DNA, Codons, Introns, Proteins or any other alignment type. Program creates a partition file for RaxML using information supplied in Ali_Type variable.

RNA_Type = True/None. If True then the ConCat program uses the alignment file for RNA structure prediction.

RNA_Struc = This variable takes RNA structure as input. If you plan to enter RNA structure then avoid using RNA_Type function. Input contains RNA structure followed by comma and RNA structure starting position. By default the starting position is set to 0 (which is the starting position of alignment).

NOTE: '.' character is not allowed in the filenames...

```
file.name.nex - Not Allowed
file_name.nex - Allowed
file-name.nes - Allowed
```

####Few commonly used operation:

Run simultaneous fasta to nexus input file format conversion and analysis with RCV calculation, Alignment Entropy calculation and spellin check function set as True for a set of alignment files that has sequence IDs stored in taxon name:
```
python ConCat.py -CA -rcv -shannon -spell -pipe -ftype fasta
```

#ConCat-Align
ConCat-Align module conducts batch sequence alignment of raw input sequence file. Direct implementation of Muscle and Mafft sequence alignment programs can be performed through ConCat-1.0.
Steps:
```
1. Place all your raw sequence data in Align directory. 
2. Run python ConCat-Align script in terminal
3. Choose between Muscle[muscle] and Mafft[mafft] program through -pkg. 
```

####Special Note for using Mafft program:  
*We provide two options for passing arguments to run Mafft.*

You can directly pass the argument through -args for conducting batch alignment with same argument.

```
python ConCat-Align.py -pkg mafft -args "--maxiterate 1000 --localpair"
```

Enter arguments within " ".


ConCat also provide an option to pass separate arguments for each alignment files. Add alignment file name and the corresponding arguments in text file and run.
```
python ConCat-Align.py -pkg mafft -sep -argf argumrntFileName.txt
```

Architecture of argument file is shown below:

```
test.afa = --maxiterate 1000 --localpair

test1.afa = --globalpair --maxiterate 1000
```

There are many more argument for alignment pre and post processing modules in ConCat. I am still working on the documentation section. 


