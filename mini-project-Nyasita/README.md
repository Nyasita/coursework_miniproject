# Welcome to **GeneFinder** 

Gene Finder contains:
* A ```MENU.py```notebook that has a GUI :smile:
* A```gene_finder.py``` python script with all code required in this program
* A```load.py``` to enable you load your data with ease
* An ```amino_acids.py``` script so that you don't have to write code for translation :wink:
* A ```data``` folder that has sample data you can use to test the functionality of the Gene Finder program 

#### Requirements
- Python version 3.7.8
- Biolinux/MacOS command line

#### About this program and its functions
The _Menu_ notebook is helpful if you would like a friendly interface, however if you are proficient, import functions from genefinder as required. This program requres that you have a fasta file on which you would like to perform biological manipulations. You can download some fasta files [here](https://www.ncbi.nlm.nih.gov/), if you need additional data.

**USAGE**

- [x] Get the complementary nucleotide: enter a nucleotide and get its complement
- [x] Get the reverse complementary sequence of DNA: enter a dna string and get its reverse complementary DNA sequence represented as a string
- [x] Find the rest_of_ORF: enter a DNA sequence that is assumed to begin with a start
        codon and returns the sequence up to but not including the
        first in frame stop codon.  If there is no in frame stop codon,
        returns the whole string.
- [x] Find all non-nested open reading frames: finds all non-nested open reading frames in the given DNA
        sequence and returns them as a list.
- [x] Find all non-nested open reading frames in all 3 possible frames: finds all non-nested open reading frames in the given DNA sequence in
        all 3 possible frames and returns them as a list.
- [x] Find all non-nested open reading frames on both strands: finds all non-nested open reading frames in the given DNA sequence on both
        strands and returns a list
- [x] Find the longest ORF on both strands: finds the longest ORF on both strands of the specified DNA and returns it
        as a string
- [x] Compute the maximum length of the longest non-coding ORF over num_trials shuffles
        of the specfied DNA sequence
- [x] Compute the Protein encoded by a sequence of DNA
    + This function
        does not check for start and stop codons (it assumes that the input
        DNA sequence represents an protein coding region).
- [x] Find the amino acid sequences that are likely coded by the specified dna
    + Returns the amino acid sequences that are likely coded by the specified dna as a list

#### Contribution
Pull requests are welcome to contribute to this project because of course team work makes the dream work!



Remember,
> Once "information" has passed into protein it cannot get out again. 
>> In more detail, the transfer of information from nucleic acid to nucleic acid, or from nucleic acid to protein may be possible, but transfer from protein to protein, or from protein to nucleic acid is impossible. Information means here the precise determination of sequence, either of bases in the nucleic acid or of amino acid residues in the protein.

Enjoy!
