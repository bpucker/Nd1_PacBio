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

python map_annotation.py \
--input_file <FULL_PATH_TO_INPUT_FILE>

WARNING: path to annotation file must be added in this script!


8) check TE overlap

check_TEs_for_overlap_with_genes.py is customized for the analysis of Nd-1. Modification of paths within this script are required to apply it to other data sets.


9) construct coverage file

python construct_cov_file.py \
--in <BAM_FILE> \
--out <OUTPUT_FILE> \
--bam_is_sorted <PREVENTS_EXTRA_SORTING_OF_BAM_FILE_By_POSITION>


10) check P/A genes
coverage_assessment.py is customized for the investigation of presence/absence variations in Arabidopsis thaliana accessions based on a read mapping against the Nd-1 assembly. Editing is required to apply this script to other data sets.



11) check gaps

check_gaps.py and analyze_sequences.py are customized for the investigation of gaps in the Col-0 reference genome sequence. Editing of paths in both scripts is necessary to apply them to other data sets.




12) read length distribution

python construct_read_lenght distribution_figure.py \
--in <FULL_PATH_TO_INPUT_FILE (FASTA)> \
--out <FULL_PATH_TO_OUTPUT_DIRECTORY>



13) tRNA gene annotation analysis
python compare_INFERNAL_to_tRNAscan.py \
--infernal <FULL_PATH_TO_INFERNAL_RESULT_FILE> \
--trnascan <FULL_PATH_TO_tRNAscan_RESULT_FILE>


14) extraction of sequence blocks from the assembly

python seqex.py \
--in <FULL_PATH_TO_INPUT_FILE> \
--out <FULL_PATH_TO_OUTPUT_FILE> \
--contig <STRING, name of contig> \
--start <INT, start of region to extract> \
--end <INT, end of region to extract> \
