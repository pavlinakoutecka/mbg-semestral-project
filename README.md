# Molecular biology and genetics course - semestral project

This project was developed by students as a part of Molecular biology and genetics course (B4M36MBG) that took place in 
the spring of 2021. It aims on the analysis of SARS-CoV-2 sequence data and mutation / protein relationships in 
individual variants.

## Source code

### Sequence preprocessing
This part is located in `sequence_preprocessing.py` file.

It aims to properly preprocess downloaded sequences - mainly
to align them and cut their edges to obtain sequences of same length. This process can be divided into following steps:
1. create random bunches of input sequences - all sequences cannot be aligned in a reasonable time
2. align the bunches of sequences using [MAFFT tool](https://mafft.cbrc.jp/alignment/software/) - for this step, 
computing power of the [RCI Cluster](https://login.rci.cvut.cz/wiki/start) was used
3. cut the gaps created by sequence aligning at both edges - to obtain sequences of same length and without gaps at
the edges, two methods for cutting the edges were developed

    - simple edge cutting algorithm - counts the number of gaps at the start and end for each sequence, finds the
     maximum of the gap at the beginning and end over all sequences, takes sequences and cuts out the maximum at the beginning
    and end
    - selective edge cutting algorithm - counts the number of gaps at the start and end for each sequence, takes 
    sequences and filters out those, which have more edge gaps than the mean value of gaps + deviation tolerance; then
    it follows the simple edge cutting algorithm process

### Mutations counting
This part is located in `mutations_counting.py` file.

The file contains three functions:
1. count_mutations(input_file, output_file, num_sequences)
    It aims to count mutations (i.e. different nucleotides) on all positions and store the results into a .csv file, omitting positions which have no mutations.

2. count_diff_sequences(input_file, output_file, num_sequences)
3. compute_nucleotides_matrix(input_file, output_file, num_sequences)

### PCA data preprocessing
This part is located in `pca_data_preprocessing.py` file.

It aims to properly preprocess aligned sequences for further PCA analysis. This process can be divided into following steps:
1. store sequences as a 2D array of characters (for significantly faster analysis)
2. count mutations on all positions using a function which is very similar to `count_mutations` described above
3. create a .csv file of sequences' names for their further identification
4. create a binary matrix (2D array) where 1 at position (r,c) indicates that the sequence r differs from a referential sequence at position c
5. store this matrix into a .csv file which is easily useable for PCA

### PCA analysis
<mark> TODO: Pavel </mark>

### Statistical analysis
<mark> TODO: Marki </mark>

### Indicating the most significant mutations
<mark> TODO: Matěj </mark>

## Data
SARS-CoV-2 Sequences were obtained from the [NCBI library](https://www.ncbi.nlm.nih.gov/sars-cov-2/). At the time of
preparing and downloading data, around 96 000 whole genome sequences of SARS-CoV-2 were available.

Code of the reference sequence (complete genome isolated ) is **NC_045512.2**.


### Data storage
As the used data are big, they cannot be store in this repository. For that reason, all the used sequences can be
found [on the ownCloud](https://owncloud.cesnet.cz/index.php/s/jXG08slIJbDptIo).

**Hierarchy of the sequences' data storage:**
- `sequences` - original sequences and their subsets/bunches
- `sequences_aligned` - aligned version of sequences stored in `sequences` folder
- `sequences_cut` - sequences from `sequences_aligned` folder with cut gaps at both edges of the sequence
- `mutations_counts` - <mark> TODO: Lukáš </mark>
- `pca_input_data` - <mark> TODO: Lukáš </mark>


## Links
[Variants and Genomic Surveillance for SARS-CoV-2](https://www.cdc.gov/coronavirus/2019-ncov/variants/index.html)

[Nextstrain SARS-CoV-2](https://nextstrain.org/sars-cov-2)
