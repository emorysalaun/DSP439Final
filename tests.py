from Exam import read_file

def test_read_file():
    expected_output = ["ATCG", "GCTA", "TACG"]
    assert read_file("test_sequences.txt") == expected_output



from Exam import find_kmers

def test_find_kmers():
    sequence = "ATCGAT"
    k = 3
    expected_output = {'ATC': {'TCG'}, 'TCG': {'CGA'}, 'CGA': {'GAT'}}
    assert find_kmers(sequence, k) == expected_output


from Exam import process_sequences

def test_process_sequences():
    sequences = ["ATCG", "GCTA"]
    k = 2
    expected_output = {'AT': {'TC'}, 'TC': {'CG'}, 'GC': {'CT'}, 'CT': {'TA'}}
    assert process_sequences(sequences, k) == expected_output

from Exam import find_minimum_k

def test_find_minimum_k():
    sequences = ["ATCGAT", "ATCGAT"]
    # Since the sequences are identical, minimum k where each k-mer leads to a unique subsequent k-mer should be 1
    expected_k = 1
    assert find_minimum_k(sequences) == expected_k

def test_find_kmers_empty():
    sequence = ""
    k = 3
    assert find_kmers(sequence, k) == {}

def test_process_sequences_empty():
    sequences = []
    k = 2
    assert process_sequences(sequences, k) == {}

def test_find_kmers_large_k():
    sequence = "ATCG"
    k = 5
    assert find_kmers(sequence, k) == {}

def test_process_sequences_large_k():
    sequences = ["ATCG", "GCTA"]
    k = 10
    assert process_sequences(sequences, k) == {}