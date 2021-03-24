## SpikeProSARS-CoV-2



#List of files in this directory

1) Structures: All the PDB structures used in the main paper (two PDB models for the full spike protein, ACE2-Spike protein complex PDB code 6M0J, and 31 Antibody-Spike protein complexes)
2) LICENSE
3) README
4) SpikePro.cpp: Main .cpp file
5) Edlib, CVParser.cpp, CVParser.hpp, P0DTC2.fasta and PIO_6.csv: Dependencies


## Usage of SpikePro

To compile the c++ program type this command:

c++ SpikePro.cpp edlib/src/edlib.cpp CSVparser.cpp -o SpikePro -I edlib/include/ -std=c++11

To run the code and predict the fitness of a viral strain 

./SpikePro TEST.fasta go

where TEST.fasta is the sequence of the viral spike protein in fasta format.  


