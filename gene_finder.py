# -*- coding: utf-8 -*-
"""
THE GENE FINDER PROGRAM

@author: LAURAH ONDARI

"""

import random #This is important for the first function shuffle_strings
from amino_acids import aa, codons, aa_table  # compiles a dictionary for translation to protein. Important for coding strand to AA function
from load import load_seq #for loading your fasta file whose data will be used for manipulation


def shuffle_string(s):
    """Shuffles the characters in the input string
        NOTE: this is a helper function, you do not
        have to modify this in any way """
    return ''.join(random.sample(s, len(s)))


def get_complement(nucleotide):
    """ Returns the complementary nucleotide

        nucleotide: a nucleotide (A, C, G, or T) represented as a string
        returns: the complementary nucleotide
    >>> get_complement('A')
    'T'
    >>> get_complement('C')
    'G'
    """
    nucleotide=nucleotide.upper() #Use upper() method to ensure even if the input is lowercase its converted to uppercase
    complimentary = {"A": "T", "C": "G", "T": "A", "G": "C"} #Dictionary of nucleotide keys and corresponding complimentary nucleotide as values
    if nucleotide not in complimentary:
        print('You entered an invalid nucleotide')
        print('Please try again')
    else:    
        return complimentary[nucleotide]
    pass


def get_reverse_complement(dna):
    """ Computes the reverse complementary sequence of DNA for the specfied DNA
        sequence

        dna: a DNA sequence represented as a string
        returns: the reverse complementary DNA sequence represented as a string
    >>> get_reverse_complement("ATGCCCGCTTT")
    'AAAGCGGGCAT'
    >>> get_reverse_complement("CCGCGTTCA")
    'TGAACGCGG'
    """
    
    reverse_compliment=[] #open an empty list to append your complimented nucleotides after every iteration
    compliment = {"A": "T", "C": "G", "T": "A", "G": "C"} 
    for i in range(len(dna)): # Using the len function to ensure all nucleotides are covered
        if dna[-i - 1] in compliment.keys(): #starting from the nucleotide at the end of the dna string whose index is -1
            reverse_compliment.append(compliment[dna[-i - 1]]) 
        else: 
            reverse_compliment.append(dna[-i - 1]) #to cater for N nucleotides that denote any nuclotide
    reverse_compliment = ''.join(reverse_compliment) #Use join to return a string
    return reverse_compliment
    pass


def rest_of_ORF(dna):
    """ Takes a DNA sequence that is assumed to begin with a start
        codon and returns the sequence up to but not including the
        first in frame stop codon.  If there is no in frame stop codon,
        returns the whole string.

        dna: a DNA sequence
        returns: the open reading frame represented as a string
    >>> rest_of_ORF("ATGTGAA")
    'ATG'
    >>> rest_of_ORF("ATGAGATAGG")
    'ATGAGA'
    """
    
    Orf = [] #Start an empty list to store the Orf
    for codon in range(0, len(dna), 3 ): #split the entire dna in blocks of 3 which is a codon
        Orf.append( dna[0+int(codon):3+int(codon)] ) #append the codons in the orf list

    StopCodons = [] #start an empty list to store stop codons found in the Orf
    stops = ['TAG', 'TAA', 'TGA']
    for StopCodon in stops:
        if StopCodon in Orf:
                StopCodons.append(StopCodon)

    indexs_stops = [] #Start an empty list to store the index of each stop codon found
    for i in StopCodons:
        indexs_stops.append(int(Orf.index(i)))
    sorted_indexs_stops = sorted(indexs_stops) #Sorts the stop codons so that once the first one is found the iteration stops there
    
    start = (Orf.index("ATG")) #Return first index of ATG found
    
    if len(sorted_indexs_stops) >= 1:
        stop = sorted_indexs_stops[0] #the only considered stop codon is the first one that is found i.e at index 0
        return "".join(Orf[start:stop]) #returns a string that starts from the start to the stop codon
    else:
        return "".join(Orf[start:]) #returns a string that starts from the start to the end of dna
    pass


def find_all_ORFs_oneframe(dna):
    """ Finds all non-nested open reading frames in the given DNA
        sequence and returns them as a list.  This function should
        only find ORFs that are in the default frame of the sequence
        (i.e. they start on indices that are multiples of 3).
        By non-nested we mean that if an ORF occurs entirely within
        another ORF, it should not be included in the returned list of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_oneframe("ATGCATGAATGTAGATAGATGTGCCC")
    ['ATGCATGAATGTAGA', 'ATGTGCCC']
    """
    
    listOfOrf = [] #Start an empty list to store all Orfs
    frame = [] 
    
    frame.append([dna[i:i + 3] for i in range(0, len(dna), 3)]) #use list comprehension to iterate over every codon made up of three characters
    
    for i in range(0,len(frame),1):
        start=0 
        while start <len(frame[i]): 
            if frame[i][start]=="ATG": #start codon is ATG
                for stop in range(start+1,len(frame[i]),1):
                    if frame[i][stop]=="TAA" or  frame[i][stop]=="TAG" or  frame[i][stop]=="TGA" : #if i contains a stop codon as stated
                        listOfOrf.append(' '.join(frame[i][start:stop])) #return orf from start codon to stop codon
                        break
                else:
                     listOfOrf.append(' '.join(frame[i][start:])) #return the whole frame from start to last codon if there is no stop codon

            start+=1 #enables state change for the while loop
            
    one_frame_orf=(",".join(listOfOrf).replace(" ","")).split(",") #use .join() to return a string and replace all spaces between characters to ensure a continuous string split the new strings with a comma
    return one_frame_orf

    pass


def find_all_ORFs(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence in
        all 3 possible frames and returns them as a list.  By non-nested we
        mean that if an ORF occurs entirely within another ORF and they are
        both in the same frame, it should not be included in the returned list
        of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs

    >>> find_all_ORFs("ATGCATGAATGTAG")
    ['ATGCATGAATGTAG', 'ATGAATGTAG', 'ATG']
    """
    
    listOfOrf = []
    frames = [] 
    
    frames.append([dna[i:i + 3] for i in range(0, len(dna), 3)]) #Starts from index 0 i.e first nucleotide
    frames.append([dna[i:i + 3] for i in range(1, len(dna), 3)]) #Starts from index 1 i.e second nucleotide
    frames.append([dna[i:i + 3] for i in range(2, len(dna), 3)]) #Starts from index 2 i.e third nucleotide
    
    for i in range(0,len(frames),1): #iterates over the lists in frames, 1 list at a time
        start=0 #Start counter
        while start <len(frames[i]):
            if frames[i][start]=="ATG":
                for stop in range(start+1,len(frames[i]),1):
                    if frames[i][stop]=="TAA" or  frames[i][stop]=="TAG" or  frames[i][stop]=="TGA" :
                        listOfOrf.append(' '.join(frames[i][start:stop])) 
                        break
                else:
                     listOfOrf.append(' '.join(frames[i][start:]))

            start+=1 #increment to ensure the while condition has new input after every iteration and will eventually stop i.e. state change
            
    all_orf=(",".join(listOfOrf).replace(" ","")).split(",") #ensures that the string that is returned does not have empty spaces and each of the strings are split using commas
    return all_orf
    
    pass


def find_all_ORFs_both_strands(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence on both
        strands.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_both_strands("ATGCGAATGTAGCATCAAA")
    ['ATGCGAATG', 'ATGCTACATTCGCAT']
    """
    
    Orf_forward =  find_all_ORFs(dna) #Calls function to find all Orfs in the dna provided
    reverseCdna =  get_reverse_complement(dna) #Calls function to reverse compliment the dna provided
    Orf_reverse = find_all_ORFs(reverseCdna) # finds all orfs in the reverse strand
    Orf_Both_strands = Orf_forward + Orf_reverse #add both lists to return one list of all orfs
    return Orf_Both_strands

    pass


def longest_ORF(dna):
    """ Finds the longest ORF on both strands of the specified DNA and returns it
        as a string
    >>> longest_ORF("ATGCGAATGTAGCATCAAA")
    'ATGCTACATTCGCAT'
    """
    
    Orf2 = find_all_ORFs_both_strands(dna) #Calls the function to find Orfs in both strands
    thelength = [] #Start an empty list to store the length of each orf
    for i in Orf2: #i equates orf in this case
        thelength.append(len(i))
    seqlen=dict((j,i) for j,i in zip(thelength,Orf2)) #Create a dictionary that has the length as key and orf as value
    thelength=sorted(thelength) #sort the lengths in ascending order
    return seqlen[thelength[-1]] #retun the value at index -1 which is the largest key

    pass


def longest_ORF_noncoding(dna, num_trials):
    """ Computes the maximum length of the longest ORF over num_trials shuffles
        of the specfied DNA sequence

        dna: a DNA sequence
        num_trials: the number of random shuffles
        returns: the maximum length longest ORF """
    
    number_of_shuffles=0 #Start counter at 0
    listofshuffled=[] #Start an empty list to store shuffled dna
    Orfs=[] #Start an empty list to store orfs from whence the longest will be computed
    while number_of_shuffles < num_trials:
        listofshuffled.append(shuffle_string(dna)) # calls the function for shuffling DNA strings
        number_of_shuffles+=1 #increments number of shuffles
        
    for shuffledDna in listofshuffled:
        Orfs.append(longest_ORF(shuffledDna)) #calls the function for finding the longest orf
    return max(Orfs,key=len)
    pass


def coding_strand_to_AA(dna):
    """ Computes the Protein encoded by a sequence of DNA.  This function
        does not check for start and stop codons (it assumes that the input
        DNA sequence represents an protein coding region).

        dna: a DNA sequence represented as a string
        returns: a string containing the sequence of amino acids encoded by the
                 the input DNA fragment

        >>> coding_strand_to_AA("ATGCGA")
        'MR'
        >>> coding_strand_to_AA("ATGCCCGCTTT")
        'MPA'
    """

    return ''.join([aa_table[dna[pos:pos + 3]] for pos in range (0, len(dna) -2, 3)])

    pass


def gene_finder(dna):
    """ Returns the amino acid sequences that are likely coded by the specified dna

        dna: a DNA sequence
        returns: a list of all amino acid sequences coded by the sequence dna.
    """
    
    mygenes=[] #Start an empty list to store genes found
    for i in find_all_ORFs_both_strands(dna):
        mygenes.append(''.join(coding_strand_to_AA(i)).replace(",","")) #use join to make it a string and replace the comas with 'nothing' to return a conutinuous string
    return mygenes

    pass

if __name__ == "__main__":
    import doctest
    doctest.testmod()
