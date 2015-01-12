# -*- coding: utf-8 -*-
# Copyright 2014 by Ambuj Kumar, Kimball-Braun lab group, University of Florida.
# All rights reserved. This code is part of the ConCat distribution and governed
# by its license. Please see the LICENSE file that should have been included
# as part of this package.
#
# Bug reports welcome: ambuj@ufl.edu 


"""
    Detects intra and inter molecular coevolving sites and proteins
    
    Program runs in three different mode - Inter, Intra and Mega
    
    
                                        Mega mode
    
    In Mega mode program calculates pairwise distance between amino acid/cds sequences 
    of all the taxon present in the input dataset and performs pearson correlation
    analysis to detect coevolution statistics.
    
    
    
    Based on
    
    Hua Zhou, Eric Jakobsson (2013)
    Predicting Protein-Protein Interaction by the Mirrortree Method: Possibilities and Limitations.
    PLoS One 8(12): e81100.
    
    
    
                                    Inter and Intra mode
    
    In Intra and Inter mode, program estimates site coevolution statistics using CAPS method.
    Program also allows user to choose between standard pearson correlation coefficient method
    implimented in CAPS program and Gobels site correlation test algorithm. 
    
    Program allows user to restrict the execution to the selected list of alignment positions
    to reduce the execution runtime in the cases where the user wants to test correlation
    between specific alignment sites.
    
    It can also detect the presence/absence of interatomic interactions, such as the hydrogen bonds,
    at the coevolving sites.
    
    
    
    Based upon 
    
    'CAPS: coevolution analysis using protein sequences'
    Fares MA, McNally D. Bioinformatics. 2006 Nov 15;22(22):2821-2.
    
    Gobel U, Sander C, Schneider R and Valencia A. (1994) Proteins, 18, 309-317.
    
    Manish C.Saraf, Gregory L.Moore and Costas D.Maranas.
    Using multiple sequence correlation analysis to characterize functionally important protein regions.
    Protein Engineering vol. 16 no. 6 pp. 397-406, 2003
    
    
    
    
                "Site Correlation" - Inter and Intra Mode
    
    BLOSUM values are normalized by the time of divergence between
    sequences. BLOSUM values (Bek) are thus weighted for the transition 
    between amino acids e and k using the time (t) since the divergence 
    between sequences i and j:
    
                      θ_ek(ij) = Bek*(t^-1)(ij)
    
    The assumption made in equation 1 is that the different types of amino 
    acid transitions (slight to radical amino acid changes) in a particular 
    site follow a Poisson distribution along time. The greater the time since 
    the divergence between sequences i and j the greater the probability of 
    having a radical change. A linear relationship is thus assumed between the 
    BLOSUM values and time. Synonymous substitutions per site (dSij) are silent 
    mutations, as they do not affect the amino acid composition of the protein. 
    These mutations are therefore neutrally fixed in the gene. Assuming that 
    synonymous sites are not saturated or under constraints, dS is proportional 
    to the time since the two sequences compared diverged. Time (t) therefore is 
    measured as dS.
    
    The next step is the estimation of the mean θ parameter for each site (θC~) 
    of the alignment, so that:
    
                   θC~ = (1/T)sigma[1->T]θ_ek(S)
    
    Here S refers to each pairwise comparison, while T stands for the total number
    of pairwise sequence comparisons, and thus:
    
                          T = N*(N-1)/2
    
    Where N is the total number of sequences in the alignment.
    
    The variability of each pairwise amino acid transition compared to that of 
    the site column is estimated as:
    
                    Dˆ_ek = (θ_ek(ij) - θC~)^2
    
    The mean variability for the corrected BLOSUM transition values is:
    
               DC= (1/T)*sigma[1->T](θ_ek(S)−θC~)2
    
    The coevolution between amino acid sites (A and B) is estimated thereafter by 
    measuring the correlation in the pairwise amino acid variability, relative to 
    the mean pairwise variability per site, between them. This covariation is 
    measured as the correlation between their Dˆ ek values. 
    
    The mean correlation coefficient and its variance are then estimated. Correlation 
    coefficients are then tested for significance under a normal distribution.
    
    The statistical power of the test is optimised by analysing sites showing:
    
                         DC > phi − 2σ_phi
    
    Here, phi is the parametric value of DC and σ is the standard deviation of phi.
    
    Pair-wise comparisons including gaps in any or both sites at any sequence are 
    excluded from the analysis.
    
    
    
    Execution commands example
    
    >>> import Coevol
    >>> from Coevol import coevol
    
    Mega mode
    >>> fileList = ["file1.nex", "file2.nex", "file3.nex"]
    >>> output = coevol(fileList, method = "Mega")
    
    or
    
    >>> output = coevol(folder="folder_name", method="Mega")
    
    
    Inter mode
    
    Default run
    >>> fileList = ["file1.nex", "file2.nex"]
    >>> coevol(fileList, method = "Inter", corr="Pearson", rCut=0.9)
    
    With pdb coordinates and selected positions
    >>> list1 = [9, 10, 11, 15, 16, 17, 18, 19, 20, 42, 43, 44]
    >>> list2 = [600, 601, 602, 603, 604, 605, 624, 625, 626]
    
    >>> coevol(fileList, type="cds", method = "Inter", pdbfile="1T5C.pdb", chain1="A", chain2="A", atom_interval1="4-339", atom_interval2="75-224", select1=list1, select2=list2, corr="Gobel")
    
    Intra mode
    >>> coevol(fileList, method = "Intra", pdbfile="1T5C.pdb", chain1="A", atom_interval1="4-339", corr="Pearson")
    
    
    """




from __future__ import print_function


try:
    from Bio._py3k import zip
    from Bio._py3k import range
    from Bio._py3k import basestring
except:
    raise ImportError("Failed to import Biopython py3k module. Please check if Biopython is installed on your system")


from functools import reduce

import sys
import math
import glob
import multiprocessing
from multiprocessing import *


try:
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.Alphabet import generic_dna
    #from Bio import pairwise2
    #from Bio.SubsMat import MatrixInfo as matlist
    from Bio.Align import MultipleSeqAlignment
except:
    raise ImportError("Failed to import Biopython. Please check if Biopython is installed on your system")


__docformat__ = "restructuredtext en"



class CoevolError(Exception):
    pass



cero_fold = {
    'ATG': 'M',
    'TGG': 'W',
    'TAA': '-',
    'TAG': '-',
    'TGA': '-'
}


two_fold = {
    'TTT': 'F',
    'TTC': 'F',
    'TAT': 'Y',
    'TAC': 'Y',
    'CAT': 'H',
    'CAC': 'H',
    'CAA': 'Q',
    'CAG': 'Q',
    'AAA': 'K',
    'AAG': 'K',
    'GAT': 'D',
    'GAC': 'D',
    'GAA': 'E',
    'GAG': 'E',
    'AAT': 'N',
    'AAC': 'N',
    'TGT': 'C',
    'TGC': 'C'
}


four_fold = {
    'ATT': 'I',
    'ATC': 'I',
    'ATA': 'I',
    'GTA': 'V',
    'GTT': 'V',
    'GTC': 'V',
    'GTG': 'V',
    'CCA': 'P',
    'CCT': 'P',
    'CCC': 'P',
    'CCG': 'P',
    'ACA': 'T',
    'ACT': 'T',
    'ACC': 'T',
    'ACG': 'T',
    'GCA': 'A',
    'GCT': 'A',
    'GCC': 'A',
    'GCG': 'A',
    'GGA': 'G',
    'GGT': 'G',
    'GGC': 'G',
    'GGG': 'G'
}


six_fold = {
    'TTA': 'L',
    'TTG': 'L',
    'CTA': 'L',
    'CTT': 'L',
    'CTC': 'L',
    'CTG': 'L',
    'TCA': 'S',
    'TCT': 'S',
    'TCC': 'S',
    'TCG': 'S',
    'AGT': 'S',
    'AGC': 'S',
    'AGA': 'R',
    'AGG': 'R',
    'CGA': 'R',
    'CGT': 'R',
    'CGC': 'R',
    'CGG': 'R'
}

purines = {
    'A': 'Pur1',
    'B': 'Pur2'
}


pyrimidines = {
    'T': 'Pyr1',
    'C': 'Pyr2'
}





codonGroups = [cero_fold, two_fold, four_fold, six_fold]
codonDict = dict()

for sections in codonGroups:
    for key, val in sections.items():
        codonDict[key] = (val)



def _average(s): return sum(s) * 1.0 / len(s)

def _variance(avg, s): return map(lambda x: (x - avg)**2, s)

def _stDev(variance): return math.sqrt(_modulate(average(variance)))

def _sort_dict_by_val(d): return sorted(d.items(), key=lambda x: x[1])

def _sort_dict_by_key(d): return sorted(d.items(), key=lambda x: x[0])

def _sort_dict_by_val_doub(d): return sorted(d.items(), key=lambda x: x[1][0])

def _sort_list_by_key(l): return sorted([x for x in l], key=lambda x: x[0])

def _get_difference(x,y): return sum(amino_x != amino_y for amino_x, amino_y in zip(x,y))

def _spliter(str, num): return [ str[start:start+num] for start in range(0, len(str), num) ]




def Blossum(a1, a2):
    
    """Blosum matrix data return function"""
    
    if a1 == '?' or a1 == 'X':
        a1 = '-'
    if a2 == '?' or a2 == 'X':
        a2 = '-'

    data_Matrix = {
                   'A': [4,  0, -2, -1, -2,  0, -2, -1, -1, -1, -1, -2, -1, -1, -1,  1,  0,  0, -3, -2,  0],
                   'C': [0,  9, -3, -4, -2, -3, -3, -1, -3, -1, -1, -3, -3, -3, -3, -1, -1, -1, -2, -2,  0],
                   'D': [-2, -3,  6,  2, -3, -1, -1, -3, -1, -4, -3,  1, -1,  0, -2,  0, -1, -3, -4, -3,  0],
                   'E': [-1, -4,  2,  5, -3, -2,  0, -3,  1, -3, -2,  0, -1,  2,  0,  0, -1, -2, -3, -2,  0],
                   'F': [-2, -2, -3, -3,  6, -3, -1,  0, -3,  0,  0, -3, -4, -3, -3, -2, -2, -1,  1,  3,  0],
                   'G': [0, -3, -1, -2, -3,  6, -2, -4, -2, -4, -3,  0, -2, -2, -2,  0, -2, -3, -2, -3,  0],
                   'H': [-2, -3, -1,  0, -1, -2,  8, -3, -1, -3, -2,  1, -2,  0,  0, -1, -2, -3, -2,  2,  0],
                   'I': [-1, -1, -3, -3,  0, -4, -3,  4, -3,  2,  1, -3, -3, -3, -3, -2, -1,  3, -3, -1,  0],
                   'K': [-1, -3, -1,  1, -3, -2, -1, -3,  5, -2, -1,  0, -1,  1,  2,  0, -1, -2, -3, -2,  0],
                   'L': [-1, -1, -4, -3,  0, -4, -3,  2, -2,  4,  2, -3, -3, -2, -2, -2, -1,  1, -2, -1,  0],
                   'M': [-1, -1, -3, -2,  0, -3, -2,  1, -1,  2,  5, -2, -2,  0, -1, -1, -1,  1, -1, -1,  0],
                   'N': [-2, -3,  1,  0, -3,  0,  1, -3,  0, -3, -2,  6, -2,  0,  0,  1,  0, -3, -4, -2,  0],
                   'P': [-1, -3, -1, -1, -4, -2, -2, -3, -1, -3, -2, -2,  7, -1, -2, -1, -1, -2, -4, -3,  0],
                   'Q': [-1, -3,  0,  2, -3, -2,  0, -3,  1, -2,  0,  0, -1,  5,  1,  0, -1, -2, -2, -1,  0],
                   'R': [-1, -3, -2,  0, -3, -2,  0, -3,  2, -2, -1,  0, -2,  1,  5, -1, -1, -3, -3, -2,  0],
                   'S': [1, -1,  0,  0, -2,  0, -1, -2,  0, -2, -1,  1, -1,  0, -1,  4,  1, -2, -3, -2,  0],
                   'T': [0, -1, -1, -1, -2, -2, -2, -1, -1, -1, -1,  0, -1, -1, -1,  1,  5,  0, -2, -2,  0],
                   'V': [0, -1, -3, -2, -1, -3, -3,  3, -2,  1,  1, -3, -2, -2, -3, -2,  0,  4, -3, -1,  0],
                   'W': [-3, -2, -4, -3,  1, -2, -2, -3, -3, -2, -1, -4, -4, -2, -3, -3, -2, -3, 11,  2,  0],
                   'Y': [-2, -2, -3, -2,  3, -3,  2, -1, -2, -1, -1, -2, -3, -1, -2, -2, -2, -1,  2,  7,  0],
                   '-': [0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0]
    }

    tag = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y', '-']

    return data_Matrix[a1][tag.index(a2)]


def hydrophob():
    
    """
        Contains Hydrophobicity and molecular weight data
        
        Column1 => Hydrophobicity data at pH 7
        Column2 => Hydrophobicity data at pH 2
        Column3 => Molecular weight
        """
    
    data_Matrix = {
                    'F': [100, 92, 165.2],
                    'I': [99, 100, 131.2],
                    'W': [97, 84, 204.2],
                    'L': [97, 100, 131.2],
                    'V': [76, 79, 117.2],
                    'M': [74, 74, 149.2],
                    'Y': [63, 49, 181.2],
                    'C': [49, 52, 121.2],
                    'A': [41, 47, 89.1],
                    'T': [13, 13, 119.1],
                    'H': [8, -42, 155.2],
                    'G': [0, 0, 75.1],
                    'S': [-5, -7, 105.1],
                    'Q': [-10, -18, 146.2],
                    'R': [-14, -26, 174.2],
                    'K': [-23, -37, 146.2],
                    'N': [-28, -41, 132.1],
                    'E': [-31, 8, 147.1],
                    'P': [-46, -46, 115.1],
                    'D': [-55, -18, 133.1],
                    '-': [0, 0, 0]

                    }

    return data_Matrix


def z_table():

    data_Matrix = {
                    0:  [0.5,       0.504,	0.508,	0.512,	0.516,	0.52,	0.5239,	0.5279,	0.5319,	0.5359],
                    1:	[0.5398,	0.5438,	0.5478,	0.5517,	0.5557,	0.5596,	0.5636,	0.5675,	0.5714,	0.5753],
                    2:	[0.5793,	0.5832,	0.5871,	0.591,	0.5948,	0.5987,	0.6026,	0.6064,	0.6103,	0.6141],
                    3:	[0.6179,	0.6217,	0.6255,	0.6293,	0.6331,	0.6368,	0.6406,	0.6443,	0.648,	0.6517],
                    4:	[0.6554,	0.6591,	0.6628,	0.6664,	0.67,	0.6736,	0.6772,	0.6808,	0.6844,	0.6879],
                    5:	[0.6915,	0.695,	0.6985,	0.7019,	0.7054,	0.7088,	0.7123,	0.7157,	0.719,	0.7224],
                    6:	[0.7257,	0.7291,	0.7324,	0.7357,	0.7389,	0.7422,	0.7454,	0.7486,	0.7517,	0.7549],
                    7:	[0.758,     0.7611,	0.7642,	0.7673,	0.7704,	0.7734,	0.7764,	0.7794,	0.7823,	0.7852],
                    8:	[0.7881,	0.791,	0.7939,	0.7967,	0.7995,	0.8023,	0.8051,	0.8078,	0.8106,	0.8133],
                    9:	[0.8159,	0.8186,	0.8212,	0.8238,	0.8264,	0.8289,	0.8315,	0.834,	0.8365,	0.8389],
                    10: [0.8413,	0.8438,	0.8461,	0.8485,	0.8508,	0.8531,	0.8554,	0.8577,	0.8599,	0.8621],
                    11:	[0.8643,	0.8665,	0.8686,	0.8708,	0.8729,	0.8749,	0.877,	0.879,	0.881,	0.883],
                    12:	[0.8849,	0.8869,	0.8888,	0.8907,	0.8925,	0.8944,	0.8962,	0.898,	0.8997,	0.9015],
                    13:	[0.9032,	0.9049,	0.9066,	0.9082,	0.9099,	0.9115,	0.9131,	0.9147,	0.9162,	0.9177],
                    14:	[0.9192,	0.9207,	0.9222,	0.9236,	0.9251,	0.9265,	0.9279,	0.9292,	0.9306,	0.9319],
                    15:	[0.9332,	0.9345,	0.9357,	0.937,	0.9382,	0.9394,	0.9406,	0.9418,	0.9429,	0.9441],
                    16:	[0.9452,	0.9463,	0.9474,	0.9484,	0.9495,	0.9505,	0.9515,	0.9525,	0.9535,	0.9545],
                    17:	[0.9554,	0.9564,	0.9573,	0.9582,	0.9591,	0.9599,	0.9608,	0.9616,	0.9625,	0.9633],
                    18:	[0.9641,	0.9649,	0.9656,	0.9664,	0.9671,	0.9678,	0.9686,	0.9693,	0.9699,	0.9706],
                    19:	[0.9713,	0.9719,	0.9726,	0.9732,	0.9738,	0.9744,	0.975,	0.9756,	0.9761,	0.9767],
                    20: [0.9772,	0.9778,	0.9783,	0.9788,	0.9793,	0.9798,	0.9803,	0.9808,	0.9812,	0.9817],
                    21:	[0.9821,	0.9826,	0.983,	0.9834,	0.9838,	0.9842,	0.9846,	0.985,	0.9854,	0.9857],
                    22:	[0.9861,	0.9864,	0.9868,	0.9871,	0.9875,	0.9878,	0.9881,	0.9884,	0.9887,	0.989],
                    23:	[0.9893,	0.9896,	0.9898,	0.9901,	0.9904,	0.9906,	0.9909,	0.9911,	0.9913,	0.9916],
                    24:	[0.9918,	0.992,	0.9922,	0.9925,	0.9927,	0.9929,	0.9931,	0.9932,	0.9934,	0.9936],
                    25:	[0.9938,	0.994,	0.9941,	0.9943,	0.9945,	0.9946,	0.9948,	0.9949,	0.9951,	0.9952],
                    26:	[0.9953,	0.9955,	0.9956,	0.9957,	0.9959,	0.996,	0.9961,	0.9962,	0.9963,	0.9964],
                    27:	[0.9965,	0.9966,	0.9967,	0.9968,	0.9969,	0.997,	0.9971,	0.9972,	0.9973,	0.9974],
                    28:	[0.9974,	0.9975,	0.9976,	0.9977,	0.9977,	0.9978,	0.9979,	0.9979,	0.998,	0.9981],
                    29:	[0.9981,	0.9982,	0.9982,	0.9983,	0.9984,	0.9984,	0.9985,	0.9985,	0.9986,	0.9986],
                    30: [0.9987,	0.9987,	0.9987,	0.9988,	0.9988,	0.9989,	0.9989,	0.9989,	0.999,	0.999]
                }

    return data_Matrix



def li_synonymous(recordObj):
    
    """returns synonymous distance(ds) as divergence time"""
    
    syn_distance = 0
    inRecStore = list()
    L = list(); Ts = list(); Tv = list(); A = list(); B = list()
    counter = 0
    while counter <= 2:
        L.append(0)
        Ts.append(0)
        Tv.append(0)
        A.append(0)
        B.append(0)
        counter = counter + 1

    for rec in recordObj:
        inRecStore.append(rec.id)
        inRec = [x for x in recordObj if x.id not in inRecStore]
        for inrec in inRec:
            codonSeq1 = _spliter(rec.seq, 3)
            codonSeq2 = _spliter(inrec.seq, 3)
            for i, codon1 in enumerate(codonSeq1):
                L, Ts, Tv = _deg_transvers(codon1, codonSeq2[i], L, Ts, Tv)
            
            A, B = _calculate_A_B(Ts, Tv, L, A, B)
            
            if L[1] + L[2] > 0:
                syn_distance = syn_distance + (float((L[1]*A[1]) + (L[2]*A[2]))/(L[1] + L[2])) + B[2]
            else:
                syn_distance = syn_distance + B[2]
            
    return syn_distance



def _deg_transvers(codon1, codon2, L, Ts, Tv):
    if codon1 in cero_fold.keys():
        L[0] += 3
        for i, nuc in enumerate(codon1):
            if nuc != codon2[i]:
                if nuc in purines.keys() and codon2[i] in purines.keys()\
                    or nuc in pyrimidines.keys() and codon2[i] in pyrimidines.keys():
                    Ts[0] = Ts[0] + 1
                elif nuc in purines.keys() and codon2[i] in pyrimidines.keys()\
                    or nuc in pyrimidines.keys() and codon2[i] in purines.keys():
                    Tv[0] = Tv[0] + 1
    
    elif codon1 in two_fold.keys():
        L[0] += 2
        L[1] += 1
        for i, nuc in enumerate(codon1):
            if nuc != codon2[i]:
                if nuc in purines.keys() and codon2[i] in purines.keys()\
                    or nuc in pyrimidines.keys() and codon2[i] in pyrimidines.keys():
                    if i <= 1:
                        Ts[0] += 1
                    else:
                        Ts[1] += 1
                elif nuc in purines.keys() and codon2[i] in pyrimidines.keys()\
                    or nuc in pyrimidines.keys() and codon2[i] in purines.keys():
                    if i <= 1:
                        Tv[0] += 1
                    else:
                        Tv[1] += 1
    
    
    elif codon1 in six_fold.keys():
        L[0] += 1
        L[1] += 1
        L[2] += 1
        for i, nuc in enumerate(codon1):
            if nuc != codon2[i]:
                if nuc in purines.keys() and codon2[i] in purines.keys()\
                    or nuc in pyrimidines.keys() and codon2[i] in pyrimidines.keys():
                    if i == 0:
                        Ts[1] += 1
                    elif i == 1:
                        Ts[0] += 1
                    elif i == 2:
                        Ts[2] += 1
                elif nuc in purines.keys() and codon2[i] in pyrimidines.keys()\
                    or nuc in pyrimidines.keys() and codon2[i] in purines.keys():
                    if i == 0:
                        Ts[1] += 1
                    elif i == 1:
                        Ts[0] += 1
                    elif i == 2:
                        Ts[2] += 1
    
    elif "---" in codon1:
        pass
    
    return L, Ts, Tv


def _calculate_A_B(Ts, Tv, L, A, B):
    P = list(); Q = list()
    counter = 0
    while counter <= 2:
        P.append(0)
        Q.append(0)
        counter = counter + 1
    
    nexCounter = 0
    while nexCounter <= 2:
        if L[nexCounter] > 0:
            P[nexCounter] += float(Ts[nexCounter])/L[nexCounter]
            Q[nexCounter] += float(Tv[nexCounter])/L[nexCounter]
        
        razon_A = 0
        razon_B = 0
        
        razon_A += (1 - Q[nexCounter] - (2*P[nexCounter]))
        razon_B += (1 - (2*Q[nexCounter]))
        
        if razon_A < 0:
            razon_A = -razon_A
        if razon_B < 0:
            razon_B = -razon_B
        
        dividendo_A = 0
        dividendo_B = 0
        
        if razon_A > 0:
            dividendo_A += float(1)/razon_A
        else:
            dividendo_A += 1
        
        if razon_B > 0:
            dividendo_B += float(1)/razon_B
        else:
            dividendo_B += 1
        
        A[nexCounter] += ((0.5 * math.log(dividendo_A)) - (0.25*math.log(dividendo_B)))
        B[nexCounter] += (0.5 * math.log(dividendo_B))

        nexCounter = nexCounter + 1
    
    return A, B




def _poisson(distance, long):
    
    pdist = 0
    if 1 - float(distance)/long != 0:
        pdist = -log(1 - float(distance)/long)
    else:
        pdist = 0

    return pdist


def _poisson_dist(seq1, se2):
    
    """Poisson distance"""
    
    dist = 0
    gap = 0
    for i, amino in enumerate(seq1):
        if seq1[i] == seq2[i] and seq1[i] != '-' and seq2[i] != '-':
            dist = dist + 1
        if seq1[i] == '-' or seq2[i] == '-':
            gap = gap + 1

    return _poisson(dist, len(seq1) - gap)


def _relative_distance(distance, sorted_dist):
    for key, val in distance.items():
        if sorted_dist[0] != 0:
            distance[key] = float(val)/sorted_dist[0]

    return distance


def _optimize(record):
    
    """Optimize distance as divergence time"""
    
    distance = dict()
    seqNameVec = list()
    for i, rec in records:
        seqNameVec.append(rec.id)
        inSeqRecord = [x for x in records if x.id not in seqNameVec]
        for inrec in inSeqRecord:
            distance[str(rec.id) + "-" + str(inrec.id)] = _poisson_dist(rec.seq, inrec.seq)

    dist = list()

    for key, val in distance:
         dist.append(val)

    sorted_dist = sorted(dist)
    optimize_dist = _relative_distance(distance, sorted_dist)

    return optimize_dist


def _thetaEK(records, optimize_dist):
    
    """
    Divergence time based optimized blosum scores for amino acid alignment
    
    @records - Alignment record object
    @optimize_dist - divergence time
        """
    
    posScore = dict()
    msaObj = MultipleSeqAlignment(records)
    for i in range(len(msaObj[1])):
        
        posVector = msaObj[:, i:i+1]
        posVectorTest = [str(x.seq) for x in posVector]

        if "*" in posVectorTest:
            raise CoevolError("Too many stop codons found in your data!!")
        elif len(set([x for x in posVectorTest if x != "?" and x != "-"])) == 1:
            i = i + 3
            continue
    
        posVectorStore[i] = (posVectorTest)
    
        store = list()
        posScore[i] = list()
        for a1 in posVector:
            store.append(a1.id)
            inObj = [x for x in posVector if x.id not in store]
            for a2 in inObj:
                posScore[i].append(float(Blossum(str(a1.seq), str(a2.seq)))/optimize_dist)

    for key, val in posScore.items():
        if float(val.count(0))/len(val) >= 0.50:
            posScore.pop(key)

    return posScore, posVectorStore


def _thetaEKcds(records, optimize_dist):
    
    """
        Divergence time based optimized blosum scores for cds alignment
        
        @records - Alignment record object
        @optimize_dist - divergence time
        """
    
    posScore = dict()
    msaObj = MultipleSeqAlignment(records)
    codonList = codonDict.keys()
    posVectorStore = dict()
    
    i = 0
    while i < len(msaObj[1]):
        posVector = msaObj[:, i:i+3]
        posVectorAmino = dict()
        for codons in posVector:
            if str(codons.seq) in codonList:
                posVectorAmino[codons.id] = codonDict[str(codons.seq)]
            elif str(codons.seq).count("-") < 3 and str(codons.seq).count("-") != 0:
                raise CoevolError("Oooops!!!!! Non frame indels found in %s sequence\n" %codons.id)
            elif "TGA" in str(codons.seq) or "TAG" in str(codons.seq) or "TAA" in str(codons.seq):
                raise CoevolError("Too many stop codons found in %s sequence!!" %codons.id)
            else:
                posVectorAmino[codons.id] = '?'

        posVectorTest = [val for key, val in posVectorAmino.items()]
        if float(len([x for x in posVectorTest if x != "?" and x != "-"]))/len(posVectorTest) < 0.3\
            or len(set([x for x in posVectorTest if x != "?" and x != "-"])) == 1:
            i = i + 3
            continue
        
        posVectorStore[i/3] = (posVectorTest)

        store = list()
        posScore[i/3] = list()

        for id1, a1 in posVectorAmino.items():
            store.append(id1)
            inObj = [x for x in posVectorAmino.items() if x[0] not in store]
            for id2, a2 in inObj:
                posScore[i/3].append(float(Blossum(a1, a2))/optimize_dist)
        
        i = i + 3
            
    for key, val in posScore.items():
        if float(val.count(0))/len(val) >= 0.50:
            posScore.pop(key)

    return posScore, posVectorStore



def _meanTheta(posScore):
    meanDist = dict()
    for key, val in posScore.items():
        meanDist[key] = _average(val)

    return meanDist


def _variability(posScore, meanDist):
    varDist = dict()
    
    for key, val in posScore.items():
        varDist[key] = list()
    
    for key, val in posScore.items():
        for inval in val:
            varDist[key].append((inval - meanDist[key])**2)

    return varDist


def _meanVar(varDist):
    meanVarDict = dict()
    for key, val in varDist.items():
        meanVarDict[key] = _average(val)

    return meanVarDict


def _thetaParam(meanVarDict):
    
    """returns alignment positions with mean variability scores greated
        than the parametric value of mean variability
        
        This step is used to optimize the statistical power of the test
        """
    
    lengthAlign = len(meanVarDict)
    thetaParamVal = float(_average(meanVarDict.values()))/len(meanVarDict)
    posCorExec = list()
    for key, val in meanVarDict.items():
        if val > thetaParamVal:
            posCorExec.append(key)

    return posCorExec



def pearson_def(x, y):
    
    """Pearson correlation test"""
    
    assert len(x) == len(y)
    n = len(x)
    assert n > 0
    avg_x = _average(x)
    avg_y = _average(y)
    diffprod = 0
    xdiff2 = 0
    ydiff2 = 0
    for idx in range(n):
        xdiff = x[idx] - avg_x
        ydiff = y[idx] - avg_y
        diffprod += xdiff * ydiff
        xdiff2 += xdiff * xdiff
        ydiff2 += ydiff * ydiff
    
    if xdiff2 == 0 or ydiff2 == 0:
        return 0
    else:
        return float(diffprod)/math.sqrt(xdiff2 * ydiff2)



def _correlation(varDist, varDist2, cory, coreNum, output):
    
    """returns site correlation scores"""
    
    store = list()
    corScoreObj = dict()
    
    if varDist2 == None:
        
        """For intra-molecular coevolution"""
        
        toolbar_width = len(varDist)*(len(varDist)-1)/2
        c = 0

        for key, val in varDist.items():
            store.append((key, val))
            inPos = [x for x in varDist.items() if x not in store]
            d1 = list(); d2 = list()
            for inkey, inval in inPos:
                c = c + 1
                for i, j in zip(val, inval):
                    if i != 0 or j != 0:
                        d1.append(i)
                        d2.append(j)
                    else:
                        continue
            
                if len(d1) <= 5:
                    continue
            
                try:
                    corScoreObj[str(key) + "-" + str(inkey)] = pearson_def(d1, d2)
                except ZeroDivisionError:
                    continue

                if cory != None and coreNum != None:
                    if cory == coreNum-1:
                        p = str((float(c)/toolbar_width)*100)[:5]
                        sys.stdout.write("\r%s%%" %p)
                        sys.stdout.flush()

                else:
                    p = str((float(c)/toolbar_width)*100)[:5]
                    sys.stdout.write("\r%s%%" %p)
                    sys.stdout.flush()

    else:
        
        """For inter-molecular coevolution"""
        
        varList = [varDist, varDist2]
        toolbar_width = len(varDist)*len(varDist2)
        c = 0

        for key, val in varDist.items():
            d1 = list(); d2 = list()
            for key2, val2 in varDist2.items():
                c = c + 1

                for i, j in zip(val, val2):
                    if i != 0 or j != 0:
                        d1.append(i)
                        d2.append(j)
                    else:
                        continue
            
                if len(d1) <= 5:
                    continue

                try:
                    corScoreObj[str(key) + "-" + str(key2)] = pearson_def(d1, d2)
                except ZeroDivisionError:
                    continue

                if cory != None and coreNum != None:
                    if cory == 0:
                        p = str((float(c)/toolbar_width)*100)[:5]
                        sys.stdout.write("\r%s%%" %p)
                        sys.stdout.flush()
                else:
                    p = str((float(c)/toolbar_width)*100)[:5]
                    sys.stdout.write("\r%s%%" %p)
                    sys.stdout.flush()

    if cory != None and coreNum != None:
        return output.put(corScoreObj)
    else:
        return corScoreObj


def gobel_correlation(val1, val2):
    
    """For a given pair of sequences (k, l), each substitution
        at a position (i or j) is associated with a similarity
        score (Xikl and Xjkl, respectively) obtained from the
        McLachlan scoring matrix. The expression used for computing
        the correlation coef®cient (rij) between two sequence positions
        (i, j) in the alignment is
        
        r = 2/N(N-1)(sigma[k=1->N-1]sigma[l=k+1->N](x_ikl - <Xi>)/si)((x_jkl - <Xj>)/si)
        
        where si and sj are the standard deviations of the scores x_ikl
        and x_jkl at positions i and j about their means <Xi> and <Xj>,
        respectively. The use of weights in computing the correlation
        coeficient has been avoided since they not only penalize genuinely
        correlated signals in a group of similar sequences but often cannot
        be quantifed in a universal fashion.
        
        Based on:
        
        Manish C.Saraf, Gregory L.Moore and Costas D.Maranas.
        Using multiple sequence correlation analysis to characterize functionally important protein regions.
        Protein Engineering vol. 16 no. 6 pp. 397-406, 2003
        
        """

    assert len(val1) == len(val2)
    assert len(val1) > 0
    assert len(val2) > 0
    
    std1 = _stDev(_variance(_average(val1), val1))
    std2 = _stDev(_variance(_average(val2), val2))

    rVal = 0
    N = len(val1)
    
    for d1, d2 in zip(val1, val2):
        rVal = rVal + (float(d1 - _average(val1))/std1)*(float(d2 - _average(val2))/std2)
    
    ret_rVal = (float(2)/(N*(N-1)))*rVal

    if ret_rVal > 1:
        ret_rVal = 1
    elif ret_rVal < -1:
        ret_rVal = -1

    return ret_rVal


def _exec_gobel_correlation(posScore1, posScore2, cory, coreNum):
    
    """Executes Gobel's alignment position substitution correlation method"""
    
    store = list()
    corScoreObj = dict()

    if posScore2 == None:
        
        """For intra-molecular coevolution"""
        
        toolbar_width = len(posScore1)*(len(posScore1)-1)
        c = 0
        
        for key, val in posScore1.items():
            store.append((key, val))
            inPos = [x for x in posScore1.items() if x not in store]
            for inkey, inval in inPos:
                c = c + 1
                d1 = list(); d2 = list()
                for i, j in zip(val, inval):
                    if i != 0 and j != 0:
                        d1.append(i)
                        d2.append(j)
                    else:
                        continue
            
                if len(d1) <= 5:
                    continue
            
                try:
                    corScoreObj[str(key) + "-" + str(inkey)] = gobel_correlation(d1, d2)
                except ZeroDivisionError:
                    continue

                if cory != None and coreNum != None:
                    if cory == 0:
                        p = str((float(c)/toolbar_width)*100)[:5]
                        sys.stdout.write("\r%s%%" %p)
                        sys.stdout.flush()
                
                else:
                    p = str((float(c)/toolbar_width)*100)[:5]
                    sys.stdout.write("\r%s%%" %p)
                    sys.stdout.flush()


    else:
        
        """For inter-molecular coevolution"""
        
        toolbar_width = len(posScore1)*len(posScore2)
        c = o
        for key, val in posScore1.items():
            for key2, val2 in posScore2.items():
                c = c + 1
                d1 = list(); d2 = list()
                for i, j in zip(val, val2):
                    if i != 0 and j != 0:
                        d1.append(i)
                        d2.append(j)
                    else:
                        continue
            
                if len(d1) <= 5:
                    continue

                try:
                    corScoreObj[str(key) + "-" + str(key2)] = gobel_correlation(d1, d2)
                except ZeroDivisionError:
                    continue
                
                if cory != None and coreNum != None:
                    if cory == 0:
                        p = str((float(c)/toolbar_width)*100)[:5]
                        sys.stdout.write("\r%s%%" %p)
                        sys.stdout.flush()
                                    
                else:
                    p = str((float(c)/toolbar_width)*100)[:5]
                    sys.stdout.write("\r%s%%" %p)
                    sys.stdout.flush()
                        
                        
    if cory != None and coreNum != None:
        return output.put(corScoreObj)
    else:
        return corScoreObj



def _meanCor(corScore):
    corrs = list()
    for key, val in corScore.items():
        if val < 0:
            corrs.append(-val)
        else:
            corrs.append(val)

    avgCorr = _average(corrs)

    return avgCorr


def _corVar(corScore, avgCorr):
    varCorr = dict()
    for key, val in corScore.items():
        varCorr[key] = (val - avgCorr)**2

    return _average([x[1] for x in varCorr.items()])


#def _zread():
#    zData = dict()
#    print(os.getcwd())
#    zVals = open("Z_scores", 'r').readlines()
#    for lines in zVals:
#        zData[float(lines.split("\t")[0])*10] = (float(lines.split("\t")[1]))
#
#    return zData



def groupy(L):
    """Create group of continuous data"""
    first = last = L[0]
    for n in L[1:]:
        if n - 1 == last:
            last = n
        else:
            yield first, last
            first = last = n
    yield first, last



def _zscore(corScore, avgCorr, varCorr):
    """calculates z score"""
    zscoreVal = dict()
    for key, val in corScore.items():
        zscoreVal[key] = float(val - avgCorr)/math.sqrt(varCorr)

    return zscoreVal


def pvalue(zval, zData):
    """calculates p value from z score"""
    if zval < 0:
        zval = -zval
    if zval == 0:
        return 0.51
    elif zval > 3:
        return 0.0001
    else:
        return(1 - float(zData[int(zval*10)][0] + zData[int(zval*10)][1])/2)


def _step_down_pval(pvals):
    """Retruns Bonferroni stepdown adjusted p-values"""
    pvList = list()
    new_pvals = dict()
    for key, val in pvals.items():
        pvList.append(val)
    
    sorted_pv = sorted(pvList)
    s1 = 0
    counter = len(sorted_pv)
    for key, val in _sort_dict_by_val(pvals):
        new_pvals[key] = max(s1, val*counter)
        counter = counter - 1
        s1 = new_pvals[key]

    return new_pvals



def _resample(corScore):
    newCorScore = [(x, corScore[x]) for x in sorted(corScore, key=corScore.get)]
    P = dict()
    for i, val in enumerate(newCorScore):
        P[val[0]] = float(i + 1)/len(newCorScore)
    
    return P


def _unlist(listSect, pos):
    retList = list()
    for inlist in listSect:
        for obj in inlist:
            if obj != str(pos):
                retList.append(int(obj))

    return retList


def _hydrophobData():
    hydroDict = dict()
    molDict = dict()
    hydroData = hydrophob()
    for key, val in hydroData.items():
        hydroDict[key + "-7"] = (float(val[0]))
        hydroDict[key + "-2"] = (float(val[1])) #[for proteins in acidic envi]
        molDict[key] = (float(val[2]))

    return hydroDict, molDict



def _hydrophobicity(data, posVectorStore2, positions, hydroDict, molDict, pH=7):
    
    """
        Calculates hydrophobicity and molecular weight correlation scores
        @data - First alignment columns
        @posVectorStore2 - Second alignment columns
        @positions - Blosum score based coevolving positions
        @hydrodict - Dictionary element containing hydrophobicity data
        @molDict - Dictionary element containing molecular weight data
        @pH - standard pH value of input dataset in test environment
        """
    
    corrScoreHyd = dict()
    corrScoreMol = dict()
    
    if posVectorStore2 == None:
        posVectorStore2 = data
    
    for pos in positions:
        pos1 = pos.split("-")[0]
        pos2 = pos.split("-")[1]
        
        data1 = list(); data2 = list()
        dataMol1 = list(); dataMol2 = list()
    
        if pH == 2:
            for x, y in zip(data[int(pos1)], posVectorStore2[int(pos2)]):
                if x != "?" and y != "?":
                    data1.append(float(hydroDict[x + "-2"])); data2.append(float(hydroDict[y + "-2"]))
                else:
                    continue

            corrScoreHyd[pos] = pearson_def(data1, data2)
    
        else:
            for x, y in zip(data[int(pos1)], posVectorStore2[int(pos2)]):
                if x != "?" and y != "?":
                    data1.append(float(hydroDict[x + "-7"])); data2.append(float(hydroDict[y + "-7"]))
                else:
                    continue
        
            corrScoreHyd[pos] = pearson_def(data1, data2)
        
        for x, y in zip(data[int(pos1)], posVectorStore2[int(pos2)]):
            if x != "?" and y != "?":
                dataMol1.append(float(molDict[x])); dataMol2.append(float(molDict[y]))
            else:
                continue

        corrScoreMol[pos] = pearson_def(dataMol1, dataMol2)

    return corrScoreHyd, corrScoreMol



def _modulate(number):
    if number < 0:
        return -number
    else:
        return number


def _varTonewVar(varDist, posCorExec):
    newVarDist = dict()
    for key, val in varDist.items():
        if key in posCorExec:
            newVarDist[key] = val

    return newVarDist


def _removeDuplPos(positions):
    for pos in positions:
        if pos.split("-")[1] + "-" + pos.split("-")[0] in positions:
            remove(pos.split("-")[1] + "-" + pos.split("-")[0])

    return positions




def PDBParser(file):
    
    """PDB file parser"""
    
    pdb_data_dict_atom = dict()
    fileData = open(file, "rU").readlines()
    counter = 0
    for lines in fileData:
        if lines[0:4] == "ATOM" or lines[0:4] == "atom":
            pdb_data_dict_atom[counter] = ([x for x in lines.split(" ")  if x != "" and x != "\n"])
            counter = counter + 1

    res_dict = dict()
    res = pdb_data_dict_atom[1][3]
    for key, val in pdb_data_dict_atom.items():
        if len(val[4]) == 1:
            res_dict[val[4]] = dict()
        elif len(val[3]) == 1:
            res_dict[val[3]] = dict()
        else:
            res_dict["A"] = dict()

    for key, val in pdb_data_dict_atom.items():
        if len(val[2]) == 7: #Splits third and fourth column if merged
            obj = val[2]
            val[2] = obj[:4]
            val.insert(3, obj[4:])
        
        pdb_data_dict_atom[key] = (val)
        
        if len(val[4]) == 1:
            res_dict[val[4]][str(key)] = dict()
        elif len(val[3]) == 1:
            res_dict[val[3]][str(key)] = dict()
        else:
            res_dict["A"][str(key)] = dict()

    for key, val in res_dict.items():
        count = 1
        for inkey, inval in val.items():
            val[inkey].clear()
            val[count] = dict()
            count = count + 1
        res_dict[key] = (val)

    counter = 1
    for key, val in pdb_data_dict_atom.items():
        if len(val) == 12:
            if len(val[4]) == 1:
                res_dict[val[4]][int(val[5])][val[2]] = (val[6:9])
            elif len(val[3]) == 1:
                res_dict[val[3]][int(val[4])][val[2]] = (val[6:9])
            else:
                res_dict["A"][int(val[4])][val[2]] = (val[6:9])
        elif len(val) == 11:
            try:
                if len(val[4]) == 1:
                    res_dict[val[4]][int(val[5])][val[2]] = (val[5:8])
                elif len(val[3]) == 1:
                    res_dict[val[3]][int(val[4])][val[2]] = (val[5:8])
                else:
                    res_dict["A"][int(val[4])][val[2]] = (val[5:8])
            except KeyError:
                if len(val[4]) == 1:
                    res_dict[val[4]][int(val[5])][val[2]] = (val[6:9])
                elif len(val[3]) == 1:
                    res_dict[val[3]][int(val[4])][val[2]] = (val[6:9])
                else:
                    res_dict["A"][int(val[4])][val[2]] = (val[6:9])
        elif len(val) == 10:
            try:
                if len(val[3]) == 1:
                    res_dict[val[3]][int(val[4])][val[2]] = (val[4:7])
                else:
                    res_dict["A"][int(val[3])][val[2]] = (val[4:7])
            except KeyError:
                try:
                    if len(val[3]) == 1:
                        res_dict[val[3]][int(val[4])][val[2]] = (val[5:8])
                    else:
                        res_dict["A"][int(val[4])][val[2]] = (val[5:8])
                except KeyError:
                    if len(val[3]) == 1:
                        res_dict[val[3]][int(val[5])][val[2]] = (val[6:9])
                    else:
                        res_dict["A"][int(val[5])][val[2]] = (val[6:9])

    for key, val in res_dict.items():
        for inkey, inval in val.items():
            if inval == {}:
                res_dict[key].pop(inkey)

    return res_dict


def _aa_distance(list1, list2):
    """
        distance = sqrt((x1-x2)^2+(y1-y2)^2+(z1-z2)^2)
        """
    return math.sqrt(sum([(float(x)-float(y))**2 for x,y in zip(list1, list2)]))


def _central_geometry(Obj):
    """Calculate central geometric coordinates of the amnio acid"""
    geo_dist = list()
    xObj = _average([float(val[0]) for key, val in Obj.items()])
    yObj = _average([float(val[1]) for key, val in Obj.items()])
    zObj = _average([float(val[2]) for key, val in Obj.items()])

    return [xObj, yObj, zObj]


def _return_dist(pos1, pos2, pdbData, chain1=None, chain2=None):
    """Calculates amino acid residue distance"""
    
    if chain2 == None:
        return _aa_distance(_central_geometry(pdbData["A"][pos1]), _central_geometry(pdbData["A"][pos2]))
    else:
        return _aa_distance(_central_geometry(pdbData[chain1][pos1]), _central_geometry(pdbData[chain2][pos2]))


def _return_dist_hb(pos1, pos2, pdbData, chain1=None, chain2=None):
    """Calculates H-bond forming atom distance"""
    
    if chain2 == None:
        return _aa_distance(pdbData["A"][pos1], pdbData["A"][pos2])
    else:
        return _aa_distance(pdbData[chain1][pos1], pdbData[chain2][pos2])


def hbond_angle(list1, list2, list3):
    
    """
        If P1 is the vertex of then
        angle = arcos((P12**2 + P13**2 - P23**2) / (2 * P12 * P13))
        where P12 is the length of the segment from P1 to P2, calculated by
        sqrt((P1x - P2x)**2 + (P1y - P2y)**2 + + (P1z - P2z)**2)
        """
    
    P12 = math.sqrt((list1[0] - list3[0])**2 + (list1[1] - list3[1])**2 + (list1[2] - list3[2])**2)
    P13 = math.sqrt((list1[0] - list2[0])**2 + (list1[1] - list2[1])**2 + (list1[2] - list2[2])**2)
    P23 = math.sqrt((list2[0] - list3[0])**2 + (list2[1] - list3[1])**2 + (list2[2] - list3[2])**2)

    return math.acos((P12**2 + P13**2 - P23**2)/(2 * P12 * P13))


def hydrogen_bond(distance):
    if distance <= 3.32:
        return 1
    else:
        return 0


def _matchy_pdb(pdb_keys, rangeVal):
    pre_val = pdb_keys[0] - 1
    ranger = list()
    diff = 0
    for val1, val2 in zip(pdb_keys, rangeVal):
        if val1 - pre_val == 1:
            ranger.append(val2+diff)
            pre_val = val1
        else:
            diff = diff + val1 - pre_val - 1
            pre_val = val1

    return ranger



def _adjust_pdb_pos(pdbData, chain1, chain2, atom_interval1, atom_interval2):
    
    """
        Adjusts amino acid positions in PDB file according to the atom_interval data
        @pdbData - Parsed PDB data dictionary
        @chain1 - First alignment chan ID
        @chain2 - Secind alignment chain ID
        @atom_interval1 - PDB file amino acid position in the first alignment
        @atom_interval2 - PDB file amino acid position in the second alignment
        
        @returns - Adjusted PDB data dictionary
        """
    
    new_pdbData = dict()
    range1 = [x for x in range(int(atom_interval1.split("-")[0]), int(atom_interval1.split("-")[1]) + 1)]
    group_pdb_data = [y for y in groupy([x[0] for x in _sort_dict_by_key(pdbData[chain1])])]
    
    if len(group_pdb_data) != 1:
        range1 = _matchy_pdb([x[0] for x in _sort_dict_by_key(pdbData[chain1])], range1)
    
    if range1[-1] - range1[0] != group_pdb_data[-1][1] - group_pdb_data[0][0]:
        raise CoevolError("First atom interval length does not match with the number of amino acid in chain %s" %chain1)
    
    for key, val in pdbData.items():
        if key == chain1:
            new_pdbData[key] = dict()
            for (inkey, inval), lKey in zip(_sort_dict_by_key(val), range1):
                new_pdbData[key][lKey] = (inval)

    if chain2 != None and atom_interval2 != None:
        range1 = [x for x in range(int(atom_interval2.split("-")[0]), int(atom_interval2.split("-")[1]) + 1)]
        group_pdb_data = [y for y in groupy([x[0] for x in _sort_dict_by_key(pdbData[chain2])])]
        
        if len(group_pdb_data) != 1:
            range2 = _matchy_pdb([x[0] for x in _sort_dict_by_key(pdbData[chain2])], range2)
        
        if len(range2) - 1 != group_pdb_data[-1][1] - group_pdb_data[0][0]:
            raise CoevolError("Second atom interval length does not match with the number of amino acid in chain %s" %chain2)

        for key, val in pdbData.items():
            if key == chain2:
                new_pdbData[key] = dict()
                for (inkey, inval), lKey in zip(_sort_dict_by_key(val), range2):
                    new_pdbData[key][lKey] = (inval)

    return new_pdbData


def PDBTest(pdbfile, newPos, chain1=None, chain2=None, atom_interval1=None, atom_interval2=None):
    pdbData = PDBParser(pdbfile)
    pdbData = _adjust_pdb_pos(pdbData, chain1, chain2, atom_interval1, atom_interval2)
    
    atom_list_hb_N = ['N', 'NH', 'NH1', 'NH2', 'NB', 'NB1', 'NB2', 'NG', 'NG1', 'NG2', 'ND', 'ND1', 'ND2', 'NE', 'NE1', 'NE2', 'NZ', 'NZ1', 'NZ2']
    atom_list_hb_O = ['O', 'OH', 'OH1', 'OH2', 'OB', 'OB1', 'OB2', 'OG', 'OG1', 'OG2', 'OD', 'OD1', 'OD2', 'OE', 'OE1', 'OE2', 'OZ', 'OZ1', 'OZ2']
    
    hbList = list()
    distance_ret = dict()
    for pos in newPos:
        if chain2 == None:
            pos1 = int(pos.split("-")[0])
            pos2 = int(pos.split("-")[1])
            
            try:
                dist = _return_dist(pos1, pos2, pdbData)
            except KeyError:
                continue
            
            for atom in atom_list_hb_N:
                for inatom in atom_list_hb_O:
                    try:
                        hb_flag = hydrogen_bond(_aa_distance(pdbData["A"][pos1][atom], pdbData["A"][pos2][inatom]))
                        if hb_flag == 1:
                            hbList.append(pos)
                            break
                    except KeyError:
                        continue


        elif chain1 != None and chain2 != None:
            dist = _return_dist(pos1, pos2, pdbData, chain1, chain2)

            for atom in atom_list_hb_N:
                for inatom in atom_list_hb_O:
                    try:
                        hb_flag = hydrogen_bond(_aa_distance(pdbData[chain1][pos1][atom], pdbData[chain2][pos2][inatom]))
                        if hb_flag == 1:
                            hbList.append(pos)
                            break
                    except KeyError:
                        continue

        distance_ret[pos] = (dist)

    return hbList, distance_ret



def _datasplit(dataObj, coreNum):
    """
        Program to split data according to the number of available cores
        @Input - List of sequence records
        @Output - List of sublist with sequence data after splitting as per the number of cores.
        """
    
    data = dict()
    totalData = len(dataObj); numGroups = (totalData/coreNum) + 1
    
    for x in range(coreNum):
        data[x] = dict()
    
    i = 0
    counter = 0
    for key, val in dataObj.items():
        if i < numGroups*(counter+1):
            data[counter][key] = (val)
        else:
            counter = counter + 1
            data[counter][key] = (val)

        i = i + 1

    return data



def _sortPos(positions):
    pl = sorted([int(x.split("-")[0]) for x in positions])
    retPos = list()
    for pos in pl:
        for inpos in positions:
            if pos == int(inpos.split("-")[0]) and inpos not in retPos:
                retPos.append(inpos)
    
    return retPos



def collect_record(fileObj):
    form_list = ["fasta", "nexus", "clustal", "stockholm", "phylip", "phylip-relaxed"]
    for i, format in enumerate(form_list):
        handle = open(fileObj, 'rU')
        records_obj = list(SeqIO.parse(handle, format))
        handle.close()
        if records_obj != []:
            break
        else:
            if i == len(form_list) - 1:
                raise CoevolError("Unknown input format. Biopython coevol module supports\
                              fasta, nexus, clustal, stockholm, phylip-interleived and phylip-relaxed file formats")
            else:
                continue

    if "|" in records_obj[1].id:
        for i, rec in enumerate(records_obj):
            records_obj[i].id = rec.id.split("|")[0]

    return records_obj



def posVecStore_update(posVectorStoreObj, selector):
    for key, val in posVectorStoreObj.items():
        posVectorStoreObj[selector[key]] = posVectorStoreObj.pop(key)

    return posVectorStoreObj



def selection(recordObj, selector):
    
    """Creates a record object with alignment positions defined via selector variable"""
    
    seqtor = dict()
    for rec in recordObj:
        seqtor[rec.id] = Seq("", generic_dna)
        for pos in selector:
            seqtor[rec.id] = seqtor[rec.id] + rec.seq[pos:pos+1]

    for i, rec in enumerate(recordObj):
        recordObj[i].seq = seqtor[rec.id]

    return recordObj



def _translate(sequence):
    retlist = list()
    for x in _spliter(sequence, 3):
        if "-" not in x and "?" not in x and "N" not in x and "R" not in x and "Y" not in x and "M" not in x and "S" not in x and "W" not in x:
            retlist.append(codonDict[str(x)])
        else:
            retlist.append('-')
    
    return "".join(retlist)




def mirrorTree(file1, file2, mat=None):
    
    """
        Utilizes pearson correlation coefficinet as a key relationship
        to define correlation between the pairwise distance among the
        taxon pairs between two protein alignments.
        
        Distances are calculated using Kimura's distance method. This is 
        a rough-and-ready distance formula for approximating PAM distance 
        by simply measuring the fraction of amino acids, p, that differs 
        between two sequences and computing the distance as (Kimura, 1983)
        
                    D = - loge ( 1 - p - 0.2 * p^2 )
                    
        This is very quick to do but has some obvious limitations. It does 
        not take into account which amino acids differ or to what amino acids 
        they change, so some information is lost. The units of the distance 
        measure are the fraction of amino acids differing, as also in the case 
        of the PAM distance. If the fraction of amino acids differing gets larger 
        than about 0.8541 the distance becomes infinite.
        
        
        Based on
        
        Hua Zhou, Eric Jakobsson (2013)
        Predicting Protein-Protein Interaction by the Mirrortree Method: Possibilities and Limitations.
        PLoS One 8(12): e81100.
        
        """
    
    records = collect_record(file1)
    records2 = collect_record(file2)
    
    try:
        assert [x.id for x in records] == [x.id for x in records2]
        records = [_translate(x.seq) for x in records]
        records2 = [_translate(x.seq) for x in records2]
    except AssertionError:
        print("\nNon matching taxon Id's found in input dataset\nPicking common taxa from %s and %s\n" %(file1, file2))
        idCommon = [x.id for x in records if x.id in [y.id for y in records2]]
        records = _sort_list_by_key( [(x.id, _translate(x.seq)) for x in records if x.id in idCommon] )
        records2 = _sort_list_by_key( [(x.id, _translate(x.seq)) for x in records2 if x.id in idCommon] )
    
    record_list = [records, records2]
    record_pairwise = dict()
    
    N = len(records)
    c = 0

    if N >= 10:
        print("Number of matching taxa = %s" %N)

        toolbar_width = len(record_list) * N
    
        print("Estimating coevlution results for %s and %s" %(file1, file2))
    
        for i, recordObj in enumerate(record_list):
            record_pairwise[i] = list()
            inRecords = list()
        
            for rec in recordObj:
                c = c + 1
                inRecords.append(rec[0])
                for inrec in recordObj:
                
                    if inrec[0] not in inRecords:
                    
                        p_dist = float(_get_difference(rec[1], inrec[1]))/len(rec[1])
                        #record_pairwise[i].append(pairwise2.align.globaldx(rec, inrec, mat)[0][-1])
                        record_pairwise[i].append( - math.log(1 - p_dist - 0.2*(p_dist**2)))
                
                    p = str((float(c)/toolbar_width)*100)[:5]
                    sys.stdout.write("\r%s%%" %p)
                    sys.stdout.flush()
        
        c_score = pearson_def(record_pairwise[0], record_pairwise[1])
        print("\nCorrelation score = %.2f\n" %c_score)
    
        return c_score

    else:
        return None



def execute_coevol(records=None,
                   records2=None,
                   itype="cds",
                   pH=7,
                   rCut=0.9,
                   fname1=None,
                   fname2=None,
                   pdbfile=None,
                   chain1=None,
                   chain2=None,
                   atom_interval1=None,
                   atom_interval2=None,
                   corr="Pearson",
                   cores=None,
                   select1=None,
                   select2=None
                   ):
    
    """
        Executes coevolution module
        @records - First alignment record
        @records2 - Second alignment record
        @itype - Test type (cds or amino acid[aa])
        @pH - corresponding body pH for the input data
        @rCut - r correlation score cutoff value
        @fname1 - First alignment filename
        @fname2 - Second alignmnet filename
        @pdbfile - PDB file
        @chain1 - PDB chain ID for the first alignment
        @chain2 - PDB chain ID for the second alignment
        @atom_interval1 - Range of PDB amino acid in reference to the first alignment
        @atom_interval2 - Range of PDB amino acid in reference to the second alignment
        @corr - Correlation test type
        @cores - Number of cores to run coevol
        
        """
    
    if records == None and records2 == None:
        raise CoevolError("Empty input alignment dataset")
    
    try:
        assert itype == "cds" or itype == "aa"
    except AssertionError:
        raise CoevolError("itype values 'cds' or 'aa' are allowed")
    
    if records2 != None:
        try:
            assert [x.id for x in records] == [x.id for x in records2]
        except AssertionError:
            try:
                records = sorted(records)
                records2 = sorted(records2)
                assert [x.id for x in records] == [x.id for x in records2]
            except AssertionError:
                print("\nNon matching taxon Id's found in input dataset\nPicking common taxa from %s and %s\n" %(fname1, fname2))
                idCommon = [x.id for x in records if x.id in [y.id for y in records2]]
                if len(idCommon) == 0:
                    raise CoevolError("No matching taxon Id found in the input dataset")
                elif len(idCommon) <= 5:
                    print("Only %s common taxon Id found in the input dataset. This could result in bad inference." %len(idCommon))

                records = [x for x in records if x.id in idCommon]
                records2 = [x for x in records2 if x.id in idCommon]
                #raise CoevolError("Non matching taxon Id's found in the input datasets")


    if select1 != None:
        select1 = sorted(select1)
        records = selection(records, select1)

    if select2 != None:
        select2 = sorted(select2)
        records2 = selection(records2, select2)

    if itype == "cds":
        optimize_dist = li_synonymous(records)
        posScore, posVectorStore = _thetaEKcds(records, optimize_dist)
        
        if records2 != None:
            optimize_dist2 = li_synonymous(records2)
            posScore2, posVectorStore2 = _thetaEKcds(records2, optimize_dist2)
        else:
            posVectorStore2 = None
            posScore2 = None

    else:
        optimize_dist = _optimize(records)
        posScore, posVectorStore = _thetaEK(records, optimize_dist)
    
        if records2 != None:
            optimize_dist2 = _optimize(records2)
            posScore2, posVectorStore2 = _thetaEK(records2, optimize_dist2)
        else:
            posVectorStore2 = None
            posScore2 = None

    if posScore == {} or posScore2 == {}:
        raise CoevolError("Selected position contains all gaps, ? and N's")

    if corr == "Pearson" or corr == "pearson":
        meanDist = _meanTheta(posScore)
        varDist = _variability(posScore, meanDist)
        meanVarDict = _meanVar(varDist)
        posCorExec = _thetaParam(meanVarDict)
        newVarDist = _varTonewVar(varDist, posCorExec)
    
    if records2 != None:
        if  corr == "Pearson" or corr == "pearson":
            meanDist2 = _meanTheta(posScore2)
            varDist2 = _variability(posScore2, meanDist2)
            meanVarDict2 = _meanVar(varDist2)
            posCorExec2 = _thetaParam(meanVarDict2)
            newVarDist2 = _varTonewVar(varDist2, posCorExec2)
    else:
        newVarDist2 = None
        varDist2 = None

    if cores != None:
        output = multiprocessing.Queue()
        coreNum = multiprocessing.cpu_count()

        if cores <= coreNum:
            coreNum = cores
        else:
            raise CoevolError("Only %s cores available on your system" %coreNum)


        print("Using %s cores for correlation test" %coreNum)
        varDistList = _datasplit(varDist, coreNum)

    if corr == "Pearson" or corr == "pearson":
        if cores != None:
            processes = [multiprocessing.Process(target=_correlation, args=(varDistList[x], varDist2, x, coreNum, output)) for x in range(coreNum)]
            for p in processes:
                p.start()
            for p in processes:
                p.join()
        
            corScore = dict()
            corScoreList = [output.get() for p in processes]
            for Obj in corScoreList:
                for key, val in Obj.items():
                    corScore[key] = (val)

        else:
            corScore = _correlation(newVarDist, newVarDist2, cory=None, coreNum=None, output=None)
                
                
    elif corr == "Gobel" or corr == "gobel":
        if cores != None:
            posScoreList = _datasplit(posScore)
            processes = [multiprocessing.Process(target=_correlation, args=(posScoreList[x], posScore2, x, coreNum, output)) for x in range(coreNum)]
            for p in processes:
                p.start()
            for p in processes:
                p.join()

            corScore = dict()
            corScoreList = [output.get() for p in processes]
            for Obj in corScoreList:
                for key, val in Obj.items():
                    corScore[key] = (val)

        else:
            corScore = _exec_gobel_correlation(posScore, posScore2, cory=None, coreNum=None, output=None)

    if select1 != None or select2 != None:
        for key, val in corScore.items():
            if  select1 != None and select2 != None:
                newkey = str(select1[int(key.split("-")[0])]) + "-" + str(select2[int(key.split("-")[1])])
            elif  select1 != None and select2 == None:
                newkey = str(select1[int(key.split("-")[0])]) + "-" + key.split("-")[1]
            elif  select2 != None and select1 == None:
                newkey =  key.split("-")[0] + "-" + str(select1[int(key.split("-")[1])])

        corScore[newkey] = corScore.pop(key)

    if select1 != None:
        posVectorStore = posVecStore_update(posVectorStore, select1)

    if select2 != None and posVectorStore2 != None:
        posVectorStore2 = posVecStore_update(posVectorStore2, select2)


    avgCorr = _meanCor(corScore)
    varCorr = _corVar(corScore, avgCorr)
    zData = z_table()
    zscoreVal = _zscore(corScore, avgCorr, varCorr)

    Pvals = _resample(corScore)
    positions = [key for key, vals in Pvals.items() if _modulate(corScore[key]) >= rCut]


    hydroDict, molDict = _hydrophobData()
    if pH == 2:
        corrScoreHyd, corrScoreMol = _hydrophobicity(posVectorStore, posVectorStore2, positions, hydroDict, molDict, pH=2)
    else:
        corrScoreHyd, corrScoreMol = _hydrophobicity(posVectorStore, posVectorStore2, positions, hydroDict, molDict, pH=7)

    #positions = [x for x in positions if _modulate(corrScoreHyd[x]) >= rCut and _modulate(corrScoreMol[x]) >= rCut]

    newPos = list()

    if records2 == None:
        positions = _removeDuplPos(positions)

    print("\n\n")
    if records2 != None:
        print("\n\n%s-%s\n"%(fname1, fname2))

    if records2 == None:
        positions = _sortPos(positions)
    
    retDict = dict()



    for pos in positions:
        if pvalue(zscoreVal[pos], zData) <= 0.05:
            print("%s: Blosum Substitution Correlation Score = %s, Hydrophobicity Correlation = %s, Molecular Weight Correlation = %s, p value = %s"\
                  %(pos, corScore[pos], corrScoreHyd[pos], corrScoreMol[pos], pvalue(zscoreVal[pos], zData)))
                  
            newPos.append(pos)
            retDict[pos] = ([corScore[pos], corrScoreHyd[pos], corrScoreMol[pos], pvalue(zscoreVal[pos], zData)])
        
        else:
            continue


    if newPos != []:
        if records2 == None:
            store = list()
            for pos in newPos:
                store.append(pos.split("-")[0])
                #store.append(pos.split("-")[1])

            store = set(store)

        else:
            store = list(); store2 = list()
            for pos in newPos:
                store.append(pos.split("-")[0])
                store2.append(pos.split("-")[1])

            store = set(store)
            store2 = set(store2)

        print("\n\n")
        if records2 == None:
            #print("%s coevolving sites found\n" %len(store))
            for pos in sorted([int(x) for x in store]):
                print("%s is coevolving with: %s" %(pos, _unlist([str(x).split("-") for x in newPos\
                                                                  if str(pos) in str(x).split("-")], pos)))
    
            if pdbfile != None:
                print("\nInitiating PDB data analysis.....\n\n")
                hbList, distance_ret = PDBTest(pdbfile, newPos, chain1, chain2, atom_interval1, atom_interval2)
                flag = False
                for pos, val in distance_ret.items():
                    print("%s - %s distance" %(pos, distance_ret[pos]))
    
                if hbList != []:
                    print("\n\nHydrogen bond interactions were found in \n\n%s"%hbList)
                else:
                    print("\n\nNo hydrogen bond interactions were found between coevolving sites\n")


        else:
            print("%s sites of %s are coevolving with %s sites of %s\n" %(len(store), fname1, len(store2), fname2))
            for pos in sorted([int(x) for x in store]):
                print("%s is coevolving with: %s" %(pos, [str(x).split("-")[1] for x in newPos if str(pos) == str(x).split("-")[0]]))
            
            if pdbfile != None:
                print("\nInitiating PDB data analysis.....\n\n")
                hbList, distance_ret = PDBTest(pdbfile, newPos, chain1, chain2, atom_interval1, atom_interval2)
                flag = False
                for pos, val in pdb_int_data.items():
                    print("%s - %s distance" %(pos, distance_ret[pos]))

                if hbList != []:
                    print("\n\nHydrogen bond interactions were found in \n\n%s"%hbList)
                else:
                    print("\n\nNo hydrogen bond interactions were found between coevolving sites\n")
                        

        with open("Coevol_log.txt", 'w') as fp:
            for key, val in retDict.items():
                fp.write("%s\t%.2f\t%.2f\t%.2f\t%.2f\n" %(key, val[0], val[1], val[2], val[3]))

                    
    else:
        retDict[pos] = (None)
        print("0 coevolving sites found\n")

    return retDict



def coevol(files=None,
           folder=None,
           itype="cds",
           pH=7,
           rCut=0.9,
           method="Intra",
           pdbfile=None,
           chain1=None,
           chain2=None,
           atom_interval1=None,
           atom_interval2=None,
           select1=None,
           select2=None,
           corr="Pearson",
           cores=None
           ):
    
    
    """
        Performs intra-alignment coevolution test
        
        @files - Input file list
        @itype - Test type (cds or amino acid[aa])
        @pH - corresponding body pH for the input data
        @rCut - r correlation score cutoff value
        @method - Type of coevol execution (Intra, Inter or Mega)
        @pdbfile - PDB file
        @chain1 - PDB chain ID for the first alignment
        @chain2 - PDB chain ID for the second alignment
        @atom_interval1 - Range of PDB amino acid in reference to the first alignment
        @atom_interval2 - Range of PDB amino acid in reference to the second alignment
        @select1 - Alignment1 selected positions to execute coevolution analysis
        @select2 - Alignment2 selected positions to execute coevolution analysis
        @corr - Correlation test type
        @cores - Number of cores to run coevol
        
        
        -> Multiprocessing mode is currently inactive
        
        
        Note - If you execute Coevol in cds mode and supply selected positions then 
        make sure you enter all three positions of codon sites that codes for the 
        specific amino acid position in the alignment.
        
        
        Execution command
        
        >>> from Coevol import coevol
        
        Mega mode
        
        >>> fileList = ["file1.nex", "file2.nex", "file3.nex"]
        >>> output = coevol(fileList, method = "Mega")
        
        or
        
        >>> output = coevol(folder="folder_name", method="Mega")
        
        
        Inter mode
        >>> coevol('file1.nex', 'file2.nex', method = "Inter", pdbfile="1T5C.pdb", chain1="A", chain2="A", 
            atom_interval1="4-339", atom_interval2="75-224", corr="Gobel")
        
        Intra mode
        >>> coevol('file.nex', method = "Intra", pdbfile="1T5C.pdb", chain1="A", atom_interval1="4-339", corr="Pearson")
        
            
        """




    if method != "Mega" and folder != None:
        raise CoevolError("Folder option can only be used in Mega method")

    if folder != None and files != None:
        raise CoevolError("Any one argument among files and folder can be used")


    if method == "Intra" or method == "intra" and len(files) > 1:
        raise CoevolError("Coevol takes one input dataset in Intra mode")


    if pdbfile != None and atom_interval1 == None or pdbfile == None and atom_interval1 != None:
        raise CoevolError("Coevol module requires atom_interval information with pdbfile")


    if pdbfile != None and method == "Inter" and atom_interval2 == None:
        raise CoevolError("Second atom interval information required to parse and adjust the position of\
                          amino acid coordinates in PDB file corresponding to the second alignment")


    if select2 != None and method == "Intra":
        raise CoevolError("Cannot select second input alignment positions in Intra mode")


    if select2 != None and file2 == None:
        raise CoevolError("Requires second input alignmnet file to select second alignment coordinates")



    if files != None:
        
        if type(files) == list:
            file1 = files[0]
        elif type(files) == str:
            file1 = files

        if file1 == None:
            raise IOError("Funtion requires input alignment dataset")

        if len(files) == 2 and type(files) == list:
            file2 = files[1]
        elif len(files) > 2 and type(files) == list:
            pass
        else:
            file2 = None


    if method == "Inter" or method == "inter" or method == "Intra" or  method == "intra":
        records = collect_record(file1)
        if records == []:
            raise CoevolError("Empty input alignment dataset")


    if method == "Inter" or method == "inter":
        if file2 == None:
            raise IOError("Coevol requires two alignment dataset as input in 'Inter' mode")

        records2 = collect_record(file2)

    elif method == "Mega" or method == "mega":
        pass

    else:
        records2 = None


    if method == "Intra" or method == "intra":
        correl = execute_coevol(records,
                       records2,
                       itype,
                       pH,
                       rCut,
                       pdbfile=pdbfile,
                       chain1=chain1,
                       chain2=chain2,
                       atom_interval1=atom_interval1,
                       atom_interval2=atom_interval2,
                       corr=corr,
                       cores=cores,
                       select1=select1,
                       select2=select2
                       )

    elif method == "Inter" or method == "inter":
        correl = execute_coevol(records,
                       records2,
                       itype,
                       pH,
                       rCut,
                       file1,
                       file2,
                       pdbfile=pdbfile,
                       chain1=chain1,
                       chain2=chain2,
                       atom_interval1=atom_interval1,
                       atom_interval2=atom_interval2,
                       corr=corr,
                       cores=cores,
                       select1=select1,
                       select2=select2
                       )


    elif method == "Mega" or method == "mega":
        if folder != None:
            files = glob.glob(folder + "/*.*")
        
        correl = dict()
        #mat = matlist.blosum62
        inFiles = list()
        for file in files:
            inFiles.append(file)
            for infile in files:
                if infile not in inFiles:
                    try:
                        retMirTree = mirrorTree(file, infile)
                        if retMirTree != None:
                            correl[file.split("/")[1] + "-" + infile.split("/")[1]] = (retMirTree)
                    except:
                        pass
                    
                    print("\n\n")



        if len(files) >= 10:
            avgCorr = _meanCor(correl)
            varCorr = _corVar(correl, avgCorr)
            zData = z_table()
            zscoreVal = _zscore(correl, avgCorr, varCorr)
        
            correl = _sort_dict_by_val(correl)

            retCorrel = dict()
            for (key, val) in correl:
                retCorrel[key] = (val, pvalue((zscoreVal[key]), zData))
            
            correl = _sort_dict_by_val_doub(retCorrel)

            with open("Coevol_log.txt", 'w') as fp:
                fp.write("Filename\tCorrelation Score\tP-Value\n")
                for key, val in correl:
                    fp.write("%s\t%.2f\t%.3f\n" %(key, val[0], val[1]))
                        
        else:
            correl = _sort_dict_by_val(correl)
            with open("Coevol_log.txt", 'w') as fp:
                fp.write("Filename\tCorrelation Score\n")
                for key, val in correl:
                    fp.write("%s\t%.2f\n" %(key, val))
                    

    return correl




"""
    In Mega mode correlation coefficients are tested for significance under
    a normal distribution when the number of input file exceeds 10. 
    
    Z = float(corr_ij - corr_mean)/square_root(variability(corr_scores))
     
    
    """





