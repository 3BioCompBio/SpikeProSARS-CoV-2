## SpikeProSARS-CoV-2



#List of files in this directory

1) AbListgit.dat: The list of 31 Antibody-Spike protein PDB complex used in the paper with their characteristics


## Usage of SpikePro

To compile the c++ program type this command:

c++ SpikePro.cpp edlib/src/edlib.cpp CSVparser.cpp -o SpikePro -I edlib/include/ -std=c++11

To run the code and predict the fitness of a viral strain 

./SpikePro TEST.fasta go

where TEST.fasta is the sequence of the viral spike protein in fasta format.  


