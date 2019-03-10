# Nd1_PacBio
These scripts are associated with the Nd-1 genome assembly based on PacBio data (https://doi.org/10.1101/407627).


1) assembly statistics

python contig_stats.py \
--input <FILENAME> \
optional: \
--min_contig_len <INTEGER> [default=500] \
--out <FULL_PATH_TO_OUTPUT_DIRECTORY>


2) sort contigs on ref

python sort_contigs_on_ref.py \
--contig_file <FULL_PATH_TO_FILE> \
--ref_file <FULL_PATH_TO_FILE> \
--output_dir <FULL_PATH_TO_DIR> \
--species <SOME_SUFFIX_FOR_PSEUDOCHROMOSOME_NAMES>


3) transfer TE annotation

annotate_TEs_in_Nd1.py is customized for the annotation of the Nd-1 genome assembly based on Araport11. Adjustments in the script are required for application to other data sets.

4) RBH identification

python identify_RBHs.py \
--prefix <FULL_PATH_TO_DIRECTORY_FOR_TMP_DATA_AND_RESULTS> \
--input1 <FULL_PATH_TO_INPUT_FILE1> \
--input2 <FULL_PATH_TO_INPUT_FILE2> \
--seq_type <'nucl'|'prot'>


5) best protein match

python match_proteins.py \
--prefix <FULL_PATH_TO_DIRECTORY_FOR_TMP_DATA_AND_RESULTS> \
--input1 <FULL_PATH_TO_INPUT_FILE1> \
--input2 <FULL_PATH_TO_INPUT_FILE2>


6) dotplot heatmap script

python dot_plot_heatmap.py \
--in1 <FULL_PATH_TO_FASTA_FILE1> \
--in2 <FULL_PATH_TO_FASTA_FILE2> \
--out <FULL_PATH_TO_OUTPUT_DIRECTORY>[.] \
OPTIONAL: \
--show	dot plot heatmap will be displayed as interactive figure \
--cite	will not run the script, but display the reference to cite for it

7) map annotation




8) check TE overlap
9) construct coverage file
10) check P/A genes
11) check gaps
12) read length distribution
