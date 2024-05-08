import sys


def read_file(filename):
    
    """
    Reads sequences from a specified file and returns them as a list of strings.

    Keyword arguments:
    filename -- the path to the file from which sequences are to be read
    """

    # Initializes an empty list to store the sequences in
    sequences = []
    # Opens the file in read mode
    with open(filename, 'r') as f:
        # Reads the file line by line
        for line in f:
            # Removes any leading or trailing whitespace
            line = line.strip()
            # Adds the sequence to the list
            sequences.append(line)
    # Returns the list of sequences.
    return sequences

def find_kmers(sequence, k):
    
    """
    Identify all substrings (k-mers) of size k in the given sequence and map each to a set of its immediate subsequent substrings.

    Keyword arguments:
    sequence -- the DNA sequence to analyze
    k -- the length of the k-mers to find
    """

    # Creates a dictionary to store Kmer subsequent substrings, the key is the original kmer
    kmers = {}
    # Goes all the way to the end of the sequence - the length of the substrings.
    for i in range(len(sequence) - k):
        # Selects the current Kmer from the sequence.
        current_kmer = sequence[i:i+k]
        # Identifies the Kmer after the sequence
        next_kmer = sequence[i+1:i+1+k]
        if current_kmer not in kmers:
            # If the Current_Kmer has not been seen before, add a set to the dict
            kmers[current_kmer] = set()
        # Fill the set with the subsequent substring.
        kmers[current_kmer].add(next_kmer)  
    # Returns the dictionary of kmers and their substrings.
    return kmers

def process_sequences(sequences, k):
    
    """
    Aggregate all k-mers and their subsequent substrings across multiple sequences and store them in a dictionary.

    Keyword arguments:
    sequences -- a list of DNA sequences
    k -- the length of the k-mers to find
    """

    # Initializes an empty dictionary to store ALL Kmers
    all_kmers = {}
    # Iterates through the list of sequences
    for sequence in sequences:
        # Finds the Kmers and their subsequent substrings
        kmers = find_kmers(sequence, k)
        # Selects the Kmer and the next Kmer from the dictionary for processing.
        for kmer, next_kmers in kmers.items():
            # If the Kmer does not yet exist in the list, add an empty set to store the substrings.
            if kmer not in all_kmers:
                all_kmers[kmer] = set()
            all_kmers[kmer].update(next_kmers)
    return all_kmers

def find_minimum_k(sequences):

    """
    Determine the smallest k such that every k-mer in the sequences has exactly one unique subsequent k-mer.

    Keyword arguments:
    sequences -- a list of DNA sequences
    """

    # Begins with K at the minimum of 1
    k = 1
    while True:
        # Creates a dictionary for all Kmers
        all_kmers = {}
        # ASsume the Kmers are unique to begin with
        unique_follow = True
        
        # Iterating through all of the sequences
        for sequence in sequences:
            # Finds the Kmers for the specific sequence
            kmers = find_kmers(sequence, k)
            for kmer, next_kmers in kmers.items():
                if kmer not in all_kmers:
                    all_kmers[kmer] = set()
                all_kmers[kmer].update(next_kmers)
                # If at any point a k-mer has more than one subsequent k-mer, set unique_follow to False
                if len(all_kmers[kmer]) > 1:
                    unique_follow = False
                    break
            if not unique_follow:
                break
        # If every k-mer has exactly one subsequent k-mer, return k
        if unique_follow:
            return k
        k += 1


if __name__ == "__main__":
    # Ensures that the argument is provided
    if len(sys.argv) < 2:
        print("Usage: python script_name.py <filename>")
        sys.exit(1)
    # ASsigns the filename to the first system arg given
    filename = sys.argv[1]
    # Reads in the sequences
    sequences = read_file(filename)
    # FInds the minimum K where each Kmer has exactly one subsequent kmer
    minimum_k = find_minimum_k(sequences)
    print(f"The smallest k where each k-mer has exactly one subsequent k-mer is: {minimum_k}")



