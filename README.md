## SpikePro. New Webserver version (July 2022)

The SpikePro algorithm has now been updated. We have also made it available via a webserver (http://babylone.3bio.ulb.ac.be/SpikePro) in which we additionaly provide experimental fitness data collected from the literature and the 3D visualization of the Spike proiten structure with the localization of all mutations. 

If you use SpikePro please cite the following papers:

[1] Pucci F, Rooman M. Prediction and Evolution of the Molecular Fitness of SARS-CoV-2 Variants: Introducing SpikePro. Viruses. 2021; 13(5):935. 

[2] Cia G, Kwasigroch J, Rooman M, Pucci F. SpikePro: a webserver to predict the fitness of SARS-CoV-2 variants. Bioinformatics. 2022; btca517.

For technical details on SpikePro we refere to these two references.  

## SpikePro. Usage

SpikePro algorithm predicts the fitness of a SARS-CoV-2 strain from the sequence of its spike protein. Given the target sequence in fasta format, the algorithm aligns it to the reference SARS-CoV-2 spike protein (Uniprot P0DTC2), list all mutations with respect to the reference and compute the fitness for each mutations as well as for the overal viral strain. You can find more details on our preprint (Pucci and Rooman, [Prediction and evolution of the molecular fitness of SARS-CoV-2 variants: Introducing SpikePro](https://www.biorxiv.org/content/10.1101/2021.04.11.439322v1), submitted).   


To compile the c++ program type this command:

c++ SpikePro.cpp edlib/src/edlib.cpp CSVparser.cpp -o SpikePro -I edlib/include/ -std=c++11

To run the code and predict the fitness of a viral strain 

./SpikePro TEST.fasta go

where TEST.fasta is the sequence of the considered variant of the SARS-CoV-2 spike protein in fasta format.  


## List of files in this directory

1) Structures: All PDB structures used in the main paper (two PDB models for the full spike protein, ACE2-Spike protein complex PDB code 6M0J, 31 Antibody-Spike protein complexes) and list of all RBD epitopes. 
2) LICENSE
3) README
4) SpikePro.cpp: Main .cpp file
5) Edlib, CVParser.cpp, CVParser.hpp, P0DTC2.fasta, PIO_6.csv, PIO_8.csv: Dependencies
6) TEST.fasta: Fasta file to use as example input 



