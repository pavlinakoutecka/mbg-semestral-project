from Bio import SeqIO
from Bio.Seq import Seq


def count_sequence_edges(input):
    """
    This function counts the number of gaps at the start and end for each sequence and collects the info into
    the start and end dictionary. It also stores the maximum of the gap at the beginning and end over all sequences.

    :param input: FASTA file with sequences or sequences themselves
    :return: dictionary with start and end gaps, maximum gap for the start and end
    """

    print("Starting the process of counting gaps in the sequences' ending...")

    fasta_sequences = SeqIO.parse(open(input), 'fasta') if isinstance(input, str) else input
    start_gaps, end_gaps, maximum_start_gap, maximum_end_gap = {}, {}, 0, 0
    for fasta in fasta_sequences:
        name, sequence = fasta.id, str(fasta.seq)
        start_idx, end_idx, start_gap, end_gap = 0, len(sequence)-1, 0, 0

        while sequence[start_idx] == '-':
            start_gap += 1
            start_idx += 1

        while sequence[end_idx] == '-':
            end_gap += 1
            end_idx -= 1

        start_gaps[name], end_gaps[name] = start_gap, end_gap
        if start_gap > maximum_start_gap:
            maximum_start_gap = start_gap
        if end_gap > maximum_end_gap:
            maximum_end_gap = end_gap

    print(f"Maximum start gap is {maximum_start_gap} and maximum end gap is {maximum_end_gap}.")

    return start_gaps, end_gaps, maximum_start_gap, maximum_end_gap


def cut_sequence_edges(input, output, maximum_start_gap, maximum_end_gap):
    """
    This function takes sequences and cuts out specified number of characters from the beginning (maximum_start_gap)
    and end (maximum_end_gap) of them. After that it saves cleaned sequences to output file.

    :param input: FASTA file with sequences or sequences themselves
    :param output: output FASTA file with sequences
    :param maximum_start_gap: number of characters to be cut from the beginning of the sequence
    :param maximum_end_gap: number of characters to be cut from the end of the sequence
    :return: None
    """

    print("Starting the process of cutting sequences' ending using maximum gap values...")

    fasta_sequences = SeqIO.parse(open(input), 'fasta') if isinstance(input, str) else input
    with open(output, 'w') as file:
        for fasta in fasta_sequences:
            sequence = str(fasta.seq)
            fasta.seq = Seq(sequence[maximum_start_gap:len(sequence) - maximum_end_gap])
            SeqIO.write(fasta, file, 'fasta')

    print(f"Cutting edges ended successfully, sequences were saved into file {output}.")


def cut_sequence_edges_selective(input, output, start_gaps, end_gaps, deviation_tolerance):
    """
    This function takes sequences, filter out sequences that have more edge gaps than the mean value + deviation
    tolerance, and cuts out specified number of characters from the beginning (maximum_start_gap) and
    end (maximum_end_gap) of the filtered sequences. After that it saves cleaned sequences to output file.

    :param input: FASTA file with sequences
    :param output: output FASTA file with sequences
    :param start_gaps: dictionary with sequences' name and number of gap characters at the beginning of sequence
    :param end_gaps: dictionary with sequences' name and number of gap characters at the end of sequence
    :param deviation_tolerance: relative value expressing the tolerance of deviation from the mean value of gaps
    :return: None
    """

    print("Starting the process of eliminating short sequences and cutting sequences' ending using maximum gap values...")

    # obtain gaps for all sequences
    fasta_sequences = SeqIO.parse(open(input), 'fasta')
    number_of_sequences, sum_start_gap, sum_end_gap = 0, 0, 0
    for fasta in fasta_sequences:
        name = fasta.id
        number_of_sequences += 1
        sum_start_gap += start_gaps[name]
        sum_end_gap += end_gaps[name]

    # compute means
    mean_start_gap = round(sum_start_gap/number_of_sequences)
    mean_end_gap = round(sum_end_gap/number_of_sequences)
    print(f"Mean start gap is {mean_start_gap} and mean end gap is {mean_end_gap}.")

    # eliminate sequences with gaps bigger than threshold value (mean + deviation tolerance)
    fasta_sequences = SeqIO.parse(open(input), 'fasta')
    number_all, number_delete = 0, 0
    valid_sequences = []
    for idx, fasta in enumerate(fasta_sequences):
        name, sequence = fasta.id, str(fasta.seq)
        number_all += 1

        if start_gaps[name] > mean_start_gap*(1+deviation_tolerance) or end_gaps[name] > mean_end_gap*(1+deviation_tolerance):
            number_delete += 1
            continue
        valid_sequences.append(fasta)

    print(f'Eliminated {number_delete} out of {number_all} sequences.')

    # use filtered sequences as input for edge counting and cutting function
    _, _, maximum_start_gap, maximum_end_gap = count_sequence_edges(valid_sequences)
    cut_sequence_edges(valid_sequences, output, maximum_start_gap, maximum_end_gap)


def cut_number_of_sequences(input_file, output_file, number):
    """
    This function takes sequence from input FASTA file, cut first N sequences and saves
    them to output FASTA file.

    :param input_file: input FASTA file with sequences
    :param output_file: output FASTA file with sequences
    :param number: number of sequences to be cut
    :return: None
    """
    print(f"Starting the process of cutting first {number} sequences from file {input_file}...")

    fasta_sequences = SeqIO.parse(open(input_file), 'fasta')
    with open(output_file, 'w') as file:
        for index, fasta in enumerate(fasta_sequences):
            if index >= number:
                break
            SeqIO.write(fasta, file, 'fasta')

    print(f"Cutting ended successfully, sequences were saved into file {output_file}!")


if __name__ == "__main__":

    # ----- PREPARE I/O FILE NAMES AND VARIABLES
    number_of_sequences = 1000
    input_file_name = f'sequences/{number_of_sequences}_aligned_sequences.fasta'
    output_file_name = f'sequences/{number_of_sequences}_aligned_cut_sequences.fasta'

    # # ----- CUT GIVEN NUMBER OF SEQUENCES
    # cut_number_of_sequences(input_file_name, output_file_name, number_of_sequences)

    # ----- CUT SEQUENCES' EDGES
    start_gaps, end_gaps, maximum_start_gap, maximum_end_gap = count_sequence_edges(input_file_name)
    # cut_sequence_edges(input_file_name, output_file_name, maximum_start_gap, maximum_end_gap)
    cut_sequence_edges_selective(input_file_name, output_file_name, start_gaps, end_gaps, 0.8)
