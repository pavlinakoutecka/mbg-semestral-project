from Bio import SeqIO
from Bio.Seq import Seq

import csv
  
def count_mutations_array(input_file, num_sequences, seq_length):
    """
    This function counts the number of mutations (changes of nucleotide)
    at each position of aligned sequences where the number of mutations > 0
    and returns the counts as a 2D array.

    :param input_file: FASTA file with aligned sequences
    :return: 2D array with mutations counts at positions with at least 1 mutation
    :save names of sequences into .csv file
    :save mutations counts into .csv file
    """
   
    fasta_sequences = SeqIO.parse(open(input_file), 'fasta')
    
    seq_names = ['-']*num_sequences
    # columns = seq_length, rows = num_sequences
    sequences = [['-' for i in range(seq_length)] for j in range(num_sequences)]

    i = 0
    for fasta in fasta_sequences:
        seq_names[i], sequences[i] = fasta.id, str(fasta.seq)
        i += 1
        
    nucleotides = []
    mutations_counts = [[] for i in range(2)]

    for i in range(seq_length):  # cols
        nucleotides.clear()
        
        for j in range(num_sequences):  # rows 
            if sequences[j][i] not in nucleotides:
                nucleotides.append(sequences[j][i])

        n = (len(nucleotides) - 1)
        if n > 0:
            mutations_counts[0].append(i)
            mutations_counts[1].append(n)

    print("Mutations on all positions counted.")

    with open(f"pca_input_data/{num_sequences}s_names.csv", 'w', newline = '') as csvfile:
        my_writer = csv.writer(csvfile, delimiter = ' ')
        my_writer.writerow(seq_names)
    csvfile.close()

    with open(f"pca_input_data/{num_sequences}s_mutations_counts.csv", 'w', newline = '') as csvfile:
        my_writer = csv.writer(csvfile, delimiter = ' ')
        my_writer.writerows(mutations_counts)
    csvfile.close()
    
    return mutations_counts

def preprocess_pca_data(input_file, output_file, num_sequences):
    """
    This function preprocesses data into .csv file used further for PCA.

    All sequences are processed into a binary matrix, where 1 at position (r,c) indicates
    that the sequence r differs from a referential sequence at position c.
    The reference is stored in the first row, so the first row equals 0 in all columns.

    :param input_file: FASTA file with aligned sequences
    ::param output_file: output .csv file containing binary matrix representing sequences
    """

    print(f"Preprocessing data for PCA from file {input_file}.")

    fasta_sequences = SeqIO.parse(open(input_file), 'fasta')
           
    for fasta in fasta_sequences:
        ref_sequence = str(fasta.seq)
        seq_length = len(ref_sequence)
        break
   
    # count mutations on positions
    mutations_counts = count_mutations_array(input_file, num_sequences, seq_length)

    fasta_sequences = SeqIO.parse(open(input_file), 'fasta')

    n_positions = len(mutations_counts[0])
    
    seq_diffs = [[0 for c in range(n_positions)] for r in range(num_sequences)]

    j = 0
    for fasta in fasta_sequences:  # rows
        sequence = str(fasta.seq)
                  
        for i in range(n_positions): # columns
            idx = mutations_counts[0][i]
            if sequence[idx] != ref_sequence[idx]:
                seq_diffs[j][i] = 1
        j += 1
                
    with open(output_file, 'w', newline = '') as csvfile:
        my_writer = csv.writer(csvfile, delimiter = ' ')
        my_writer.writerow(mutations_counts[0])
        my_writer.writerows(seq_diffs)                   

    csvfile.close()

    print(f"Preprocessed PCA data saved into file {output_file}.")


if __name__ == "__main__":

    # ----- PREPARE INPUT AND OUTPUT FILE NAME AND VARIABLES
    number_of_sequences = 1000
    input_file_name = f'sequences/{number_of_sequences}_aligned_sequences.fasta'
    output_file_name = f'pca_input_data/{number_of_sequences}s_pca_input.csv'

    # ----- CREATE FILE WITH DESIRED INFORMATIONS ABOUT SEQUENCES FOR FURTHER ANALYSIS USING PCA
    preprocess_pca_data(input_file_name, output_file_name, number_of_sequences)


