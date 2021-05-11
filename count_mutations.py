from Bio import SeqIO
from Bio.Seq import Seq

def count_mutations(input_file, output_file, th):
    """
    This function counts the number of different mutations at each position of aligned sequences.

    :param input_file: FASTA file with aligned sequences
    ::param output_file: output file with mutations counts
    """

    print(f"Starting the process of counting different mutations in file {input_file}...")

    fasta_sequences = SeqIO.parse(open(input_file), 'fasta')
           
    for fasta in fasta_sequences:
        sequence = str(fasta.seq)
        seq_length = len(sequence)
        break

    # print(f"Length of sequences is {seq_length}")
   
    nucleotides = []

    file = open(output_file,"a")
    file.write("Position mutations_count\n")
    
    for idx in range(seq_length):
        nucleotides.clear()
        fasta_sequences = SeqIO.parse(open(input_file), 'fasta')    # reset the generator!
        
        for fasta in fasta_sequences:
            sequence = str(fasta.seq)
            #print(f"Position {idx}, nucleotide {sequence[idx]}")
            
            if sequence[idx] not in nucleotides:
                nucleotides.append(sequence[idx])
                
        mutations_count = len(nucleotides) - 1

        if mutations_count >= th:
            file.write(f"{idx} {mutations_count}\n")

        if idx % 1000 == 0:
            print(f"Progress: position {idx} out of {seq_length - 1} counted.")

    file.close()

    print(f"Mutations succesfully counted and stored in {output_file}.")


def count_diff_sequences(input_file, output_file, th):
    """
    This function counts aligned sequences which differ from the reference at given position,
    for all positions in aligned sequences.

    :param input_file: FASTA file with aligned sequences
    ::param output_file: output file with different sequences counts
    """

    print(f"Starting the process of counting different sequences in file {input_file}...")

    fasta_sequences = SeqIO.parse(open(input_file), 'fasta')
           
    for fasta in fasta_sequences:
        reference_seq = str(fasta.seq)
        seq_length = len(reference_seq)
        break

    seq_length = 10001
    # print(f"Length of sequences is {seq_length}")

    file = open(output_file,"a")
    file.write("Position diff_sequences_count\n")

    for idx in range(seq_length):
        diff_seq_count = 0
        fasta_sequences = SeqIO.parse(open(input_file), 'fasta')    # reset the generator!
        
        for fasta in fasta_sequences:
            sequence = str(fasta.seq)
            #print(f"Position {idx}, nucleotide {sequence[idx]}")
            
            if sequence[idx] != reference_seq[idx]:
                diff_seq_count += 1

        if diff_seq_count >= th:
            file.write(f"{idx} {diff_seq_count}\n")

        if idx % 1000 == 0:
            print(f"Progress: position {idx} out of {seq_length - 1} counted.")

    file.close()

    print(f"Different sequences succesfully counted and stored in {output_file}.")        


def compute_mutations_matrix(input_file, output_file):
    """
    This function counts occurrences of each nukleotide at each position in aligned sequences.

    :param input_file: FASTA file with aligned sequences
    ::param output_file: output file containing matrix of mutations occurrences
    """

    print(f"Starting the process of computing mutations matrix for file {input_file}...")

    fasta_sequences = SeqIO.parse(open(input_file), 'fasta')
           
    for fasta in fasta_sequences:
        sequence = str(fasta.seq)
        seq_length = len(sequence)
        break

    # print(f"Length of sequences is {seq_length}")

    file = open(output_file,"a")
    file.write("Position - a c g t other\n")

    for idx in range(seq_length):
        dash_cnt, a_cnt, c_cnt, g_cnt, t_cnt, other_cnt = 0, 0, 0, 0, 0, 0
        fasta_sequences = SeqIO.parse(open(input_file), 'fasta')    # reset the generator!
        
        for fasta in fasta_sequences:
            sequence = str(fasta.seq)
            #print(f"Position {idx}, nucleotide {sequence[idx]}")

            nucleotide = sequence[idx]    # "-", "a", "c", "g", "t", "other"

            if nucleotide == "-": dash_cnt += 1
            elif nucleotide == "a": a_cnt += 1
            elif nucleotide == "c": c_cnt += 1
            elif nucleotide == "g": g_cnt += 1
            elif nucleotide == "t": t_cnt += 1
            else: other_cnt += 1
        
        file.write(f"{idx} {dash_cnt} {a_cnt} {c_cnt} {g_cnt} {t_cnt} {other_cnt}\n")

        if idx % 1000 == 0:
            print(f"Progress: position {idx} out of {seq_length - 1} counted.")

    file.close()

    print(f"Mutation occurences succesfully counted and stored in {output_file}.")
    

if __name__ == "__main__":

    c = 2   # what to count: 1 = mutations, 2 = different sequence, 3 = mutations occurences

    th_mutations = 0  # threshold for mutations count
    th_diff_seq = 5  # threshold for diff_sequences count
    
    # ----- PREPARE INPUT FILE NAME AND VARIABLES
    number_of_sequences = 100
    input_file_name = f'sequences/{number_of_sequences}_aligned_sequences.fasta'

    if c == 1:
        # ----- COUNT DIFFERENT MUTATIONS IN SEQUENCES
        output_file_name = f'mutations_counts/{number_of_sequences}s_mutations_counts_th{th_mutations}.txt'
        count_mutations(input_file_name, output_file_name, th_mutations)
        
    elif c == 2:
        # ----- COUNT POSITION-DIFFERENT SEQUENCES
        output_file_name = f'mutations_counts/{number_of_sequences}s_diff_seq_counts_th{th_diff_seq}.txt'
        count_diff_sequences(input_file_name, output_file_name, th_diff_seq)

    else:
        # ----- COUNT MATRIX OF MUTATIONS OCCURRENCES
        output_file_name = f'mutations_counts/{number_of_sequences}s_mutations_occurrences.txt'
        compute_mutations_matrix(input_file_name, output_file_name)




