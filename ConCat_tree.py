from Bio.Phylo.Applications import RaxmlCommandline
from src.Functions import Convert
import subprocess
import os

Convert('nexus', 'phylip-relaxed', 'Combined.nex')

with open('RaxML/Partition.txt', 'w') as fp:
    fdata = open('Partition.txt', 'r').readlines()
    for data in fdata:
        fp.write(data)

print "Enter parsimony seed value \n"
bootInp = input('\n')
print "Enter bootstrap value \n"
nrep = input('\n')
raxml_cline = RaxmlCommandline(sequences="Combined.phy", model="GTRCAT", name="Test", partition_filename='Partition.txt', parsimony_seed=bootInp, num_replicates=nrep)
os.chdir('RaXML')
subprocess.call(raxml_cline(), shell=True)
print "All Done \n"