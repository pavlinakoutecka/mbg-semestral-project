from Bio import SeqIO
from Bio.Seq import Seq

import csv

def count_mutations(input_file, output_file, num_sequences):
    """
    This function counts the number of mutations (changes of nucleotide)
    at each position of aligned sequences where the number of mutations > 0.

    :param input_file: FASTA file with aligned (and cut) sequences
    :save mutations counts into .csv file
    """

    print(f"Counting mutations in file {input_file}")
   
    fasta_sequences = SeqIO.parse(open(input_file), 'fasta')

    # get the length of (aligned and cut) sequences
    for fasta in fasta_sequences:
        seq_length = len(fasta.seq)
        break

    # reset the fasta generator
    fasta_sequences = SeqIO.parse(open(input_file), 'fasta')

    # store fasta sequences as 2D array of characters (significantly faster analysis)
    # columns = seq_length, rows = num_sequences
    sequences = [['-' for i in range(seq_length)] for j in range(num_sequences)]

    i = 0
    for fasta in fasta_sequences:
        sequences[i] = str(fasta.seq)
        i += 1
        
    nucleotides = []
    mutations_counts = [[] for j in range(2)]

    # count different characters at each position (column)
    # and store the numbers with indices into 2D array
    for i in range(seq_length):  # cols
        nucleotides.clear()
        
        for j in range(num_sequences):  # rows 
            if sequences[j][i] not in nucleotides:
                nucleotides.append(sequences[j][i])

        n = (len(nucleotides) - 1)
        if n > 0:
            mutations_counts[0].append(i)
            mutations_counts[1].append(n)

    print(f"Mutations counted and stored in {output_file}")

    # store the array into .csv file
    with open(output_file, 'w', newline = '') as csvfile:
        my_writer = csv.writer(csvfile, delimiter = ' ')
        my_writer.writerows(mutations_counts)
    csvfile.close()


def count_diff_sequences(input_file, output_file, num_sequences):
    """
    This function counts aligned sequences which differ from the reference at given position,
    for that positions where the number of mutations > 0.

    :param input_file: FASTA file with aligned (and cut) sequences
    ::param output_file: output file with different sequences counts
    """

    print(f"Counting different sequences in file {input_file}")

    fasta_sequences = SeqIO.parse(open(input_file), 'fasta')

    # get the referential sequence (in our data always first sequence in the file)
    # and the length of (aligned and cut) sequences
    for fasta in fasta_sequences:
        reference_seq = str(fasta.seq)
        seq_length = len(reference_seq)
        break

    # reset the fasta generator
    fasta_sequences = SeqIO.parse(open(input_file), 'fasta')
         
    # store fasta sequences as 2D array of characters (significantly faster analysis)
    # columns = seq_length, rows = num_sequences
    sequences = [['-' for i in range(seq_length)] for j in range(num_sequences)]

    i = 0
    for fasta in fasta_sequences:
        sequences[i] = str(fasta.seq)
        i += 1

    # for all position where the number of mutations > 0,
    # count sequences which, on given position,
    # have different character then reference
        
    diff_seq_counts = [[] for j in range(2)]

    for i in range(seq_length): #cols
        n = 0
            
        for j in range(num_sequences): #rows
            seq = sequences[j]
            if seq[i] != reference_seq[i]:
                n += 1

        if n > 0:
            diff_seq_counts[0].append(i)
            diff_seq_counts[1].append(n)

    # store the array into .csv file
    with open(output_file, 'w', newline = '') as csvfile:
        my_writer = csv.writer(csvfile, delimiter = ' ')
        my_writer.writerows(diff_seq_counts) 
    csvfile.close()

    print(f"Different sequences counted and stored in {output_file}")        


def compute_nucleotides_matrix(input_file, output_file, num_sequences):
    """
    This function counts occurrences of each nucleotide at each position.
    
    Deletion/isertion represented by '-' is considered as one of the nucleotides,
    hence the number of different nucleotides is 6 (a, c, g, t, -, and other).

    :param input_file: FASTA file with aligned (and cut) sequences
    ::param output_file: output file containing matrix of nucleotides' occurrences

    Counted occurencies are stored into a matrix where different nucleotides
    are represented as rows and positions are represented as columns.
    Order of the nucleotides is as follows: a, c, g, t, -, other.
    """

    print(f"Counting nucleotides' occurencies in file {input_file}...")

    fasta_sequences = SeqIO.parse(open(input_file), 'fasta')

    # get the length of (aligned and cut) sequences       
    for fasta in fasta_sequences:
        seq_length = len(fasta.seq)
        break

    # reset the fasta generator
    fasta_sequences = SeqIO.parse(open(input_file), 'fasta')
         
    # store fasta sequences as 2D array of characters (significantly faster analysis)
    # columns = seq_length, rows = num_sequences
    sequences = [['-' for i in range(seq_length)] for j in range(num_sequences)]

    i = 0
    for fasta in fasta_sequences:
        sequences[i] = str(fasta.seq)
        i += 1

    # compute nucleotides' occurences and store the numbers in a matrix    
    nucleotides_matrix = [[0 for i in range(seq_length)] for j in range(6)]

    for i in range(seq_length): #cols
        a, c, g, t, dash, other = 0, 0, 0, 0, 0, 0
      
        for j in range(num_sequences): #rows
            nucleotide = sequences[j][i]

            if nucleotide == "a": a += 1
            elif nucleotide == "c": c += 1
            elif nucleotide == "g": g += 1
            elif nucleotide == "t": t += 1
            elif nucleotide == "-": dash += 1
            else: other += 1

        nucleotides_matrix[0][i] = a
        nucleotides_matrix[1][i] = c
        nucleotides_matrix[2][i] = g
        nucleotides_matrix[3][i] = t
        nucleotides_matrix[4][i] = dash
        nucleotides_matrix[5][i] = other        

    # store the array into .csv file
    with open(output_file, 'w', newline = '') as csvfile:
        my_writer = csv.writer(csvfile, delimiter = ' ')
        my_writer.writerows(nucleotides_matrix) 
    csvfile.close()
   
    print(f"Nucleotides' occurences counted and stored in {output_file}.")

    
if __name__ == "__main__":
    
    # ----- PREPARE INPUT FILE NAME AND VARIABLES
    number_of_sequences = 100
    input_file_name = f'sequences/{number_of_sequences}_cut_sequences.fasta'
    
    # ----- COUNT MUTATIONS IN SEQUENCES
    #output_file_name = f"mutations_counts/{number_of_sequences}s_mutations_counts.csv"
    #count_mutations(input_file_name, output_file_name, number_of_sequences)
    
    # ----- COUNT DIFFERENT SEQUENCES
    #output_file_name = f'mutations_counts/{number_of_sequences}s_diff_seq_counts.csv'
    #count_diff_sequences(input_file_name, output_file_name, number_of_sequences)
        
    # ----- COUNT MATRIX OF NUCLEOTIDES' OCCURRENCES
    output_file_name = f'mutations_counts/{number_of_sequences}s_nucleotides_occurrences.csv'
    compute_nucleotides_matrix(input_file_name, output_file_name, number_of_sequences)
    



