## SpikeProSARS-CoV-2

## Usage of SpikePro

SpikePro algorithm predicts the fitness of a SARS-CoV-2 strain from the sequence of its spike protein. Given the target sequence in fasta format, the algorithm aligns it to the reference SARS-CoV-2 spike protein (Uniprot P0DTC2), list all mutations with respect to the reference and compute the fitness for each mutations as well as for the overal viral strain. You can find more details on our preprint (Pucci and Rooman, [Prediction and evolution of the molecular fitness of SARS-CoV-2 variants: Introducing SpikePro](https://www.biorxiv.org/content/10.1101/2021.04.11.439322v1), submitted).   


To compile the c++ program type this command:

c++ SpikePro.cpp edlib/src/edlib.cpp CSVparser.cpp -o SpikePro -I edlib/include/ -std=c++11

To run the code and predict the fitness of a viral strain 

./SpikePro TEST.fasta go

where TEST.fasta is the sequence of the considered variant of the SARS-CoV-2 spike protein in fasta format.  


#List of files in this directory

1) Structures: All the PDB structures used in the main paper (two PDB models for the full spike protein, ACE2-Spike protein complex PDB code 6M0J, and 31 Antibody-Spike protein complexes)
2) LICENSE
3) README
4) SpikePro.cpp: Main .cpp file
5) Edlib, CVParser.cpp, CVParser.hpp, P0DTC2.fasta and PIO_6.csv: Dependencies
6) TEST.fasta: Fasta file to use as example input 



