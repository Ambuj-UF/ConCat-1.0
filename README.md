#ConCat v1.0

###BioPython is not a prerequisite for ConCat anymore
##Coevolution module under development. Check Coevolver folder for currently added coevolution functionalities. Suggestions are invited. 

###Read Manual for detail instructions

##Requirements:

BioPython and python 2.7 and above

Handling different file formats: Importance of ConCat lies in its ability to store and generate rich annotation data. It allows user to define several features in alignment file through ConCat block. Thus, it is important for ConCat to use nexus files as input for extended functionality. To make this easier, ConCat -convert function scans alignment files in any file format and converts it into nexus file format. 

Alignment file present in Input directory contains database ID's embedded in  the taxon names. Use -pipe argument to concatenate this dataset.

``` 
python ConCat-build.py -pipe
```

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

Ali_Type: This variable is used to define alignment type. It can either be DNA, Codons, Introns, Proteins or any other alignment type. Program creates a partition file for RaxML using information supplied in Ali_Type variable.

RNA_Type = True/None. If True then the ConCat program uses the alignment file for RNA structure prediction.

RNA_Struc = This variable takes RNA structure as input. If you plan to enter RNA structure then avoid using RNA_Type function. Input contains RNA structure followed by comma and RNA structure starting position. By default the starting position is set to 0 (which is the starting position of alignment).

NOTE: '.' character is not allowed in filename...

```
file.name.nex - Not Allowed
file_name.nex - Allowed
file-name.nes - Allowed
```

#ConCat-import
ConCat-import module conducts batch sequence alignment of raw input sequence data. Direct implementation of Muscle and Mafft sequence alignment programs can be performed through ConCat-1.0.

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

Enter arguments within double quotes.


ConCat also provide an option to pass separate arguments for each alignment files. Add alignment file name and the corresponding arguments in text file and run.
```
python ConCat-Align.py -pkg mafft -sep -argf argumrntFileName.txt
```

Architecture of argument file is shown below:

```
test.afa = --maxiterate 1000 --localpair

test1.afa = --globalpair --maxiterate 1000
```

There are many more arguments that allows user to fetch CDS and mRNA sequences from NCBI database.

#ConCat-build

This module performs ultrafast alignment concatenation and processing. Read Manual for detail instructions.

#ConCat-analyze

ConCat-analyze module allows user to perform various sequence alignment based analysis to retrieve biologically significant information. Read Manual for the detailed instructions.


