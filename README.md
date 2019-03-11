[![DOI](https://zenodo.org/badge/168144094.svg)](https://zenodo.org/badge/latestdoi/168144094)


# Nd1_PacBio
These scripts are associated with the Nd-1 genome assembly based on PacBio data (https://doi.org/10.1101/407627).


1) assembly statistics: this script generates statistics about a given FASTA file and allows trimming of the assembled sequences therein.

python contig_stats.py \
--input <FILENAME> \
OPTIONAL: \
--min_contig_len <INTEGER> [default=500] \
--out <FULL_PATH_TO_OUTPUT_DIRECTORY>


2) sort contigs on ref: this script runs a BLASTn search against a reference genome sequence and sorts the contigs based on the position of the best BLAST hit per contig.

python sort_contigs_on_ref.py \
--contig_file <FULL_PATH_TO_FILE> \
--ref_file <FULL_PATH_TO_FILE> \
--output_dir <FULL_PATH_TO_DIRECTORY> \
--species <SOME_SUFFIX_FOR_PSEUDOCHROMOSOME_NAMES>


3) transfer TE annotation: this script uses the Araport11 TE annotation to annotate elements in a new assembly. BLASTn is used for the mapping of sequences. Overlapping hits are filtered to allow only one hit per position.

annotate_TEs_in_Nd1.py is customized for the annotation of the Nd-1 genome assembly based on Araport11. Adjustments in the script are required for application to other data sets.


4) RBH identification: this script identifies reciprocal best BLAST hits (RBHs) between two sequence sets. The comparison can be compared on the DNA or peptide level.

python identify_RBHs.py \
--prefix <FULL_PATH_TO_DIRECTORY_FOR_TMP_DATA_AND_RESULTS> \
--input1 <FULL_PATH_TO_INPUT_FILE1> \
--input2 <FULL_PATH_TO_INPUT_FILE2> \
--seq_type <'nucl'|'prot'>


5) best protein match: this script identifies the best match in another set of peptide sequences. First, RBHs are identified between both sets. Next, remaining query sequences are assigned to the unidirectional best BLASTp hit.

python match_proteins.py \
--prefix <FULL_PATH_TO_DIRECTORY_FOR_TMP_DATA_AND_RESULTS> \
--input1 <FULL_PATH_TO_INPUT_FILE1> \
--input2 <FULL_PATH_TO_INPUT_FILE2>


6) dotplot heatmap script: this script compares two FASTA files by running BLAST to identify similarity. Best matches for small chunks of the one sequence against the other one are used to generate a dot plot. The similarity of matches is represented by the colour of dots.

python dot_plot_heatmap.py \
--in1 <FULL_PATH_TO_FASTA_FILE1> \
--in2 <FULL_PATH_TO_FASTA_FILE2> \
--out <FULL_PATH_TO_OUTPUT_DIRECTORY>[.] \
OPTIONAL: \
--show	(dot plot heatmap will be displayed as interactive figure) \
--cite	(will not run the script, but display the reference to cite for it)


7) map annotation: this script maps the Araport11 annotation to a file with Arabidopsis gene identifiers (AGIs).

python map_annotation.py \
--input_file <FULL_PATH_TO_INPUT_FILE>

WARNING: path to annotation file must be added in this script!


8) check TE overlap: this script checks annotated TEs in a genome sequence for overlap with annotated protein coding genes. Results are visualized as histogram.

check_TEs_for_overlap_with_genes.py is customized for the analysis of Nd-1. Modification of paths within this script are required to apply it to other data sets.


9) construct coverage file: this script generate a simple text file which contains the coverage of each position in a genome sequence.

python construct_cov_file.py \
--in <BAM_FILE> \
--out <OUTPUT_FILE> \
--bam_is_sorted (PREVENTS_EXTRA_SORTING_OF_BAM_FILE_BY_POSITION)


10) check P/A genes: this script allows pan-genomic analysis in Arabidopsis thaliana by assessing the presence/absence of genes based on read mapping results.

coverage_assessment.py is customized for the investigation of presence/absence variations in Arabidopsis thaliana accessions based on a read mapping against the Nd-1 assembly. Editing is required to apply this script to other data sets.



11) check gaps: this script identifies gaps in the Col-0 reference sequence, which are spanned by the Nd-1 assembly.

check_gaps.py and analyze_sequences.py are customized for the investigation of gaps in the Col-0 reference genome sequence. Editing of paths in both scripts is necessary to apply them to other data sets.




12) read length distribution: this script assesses the length distribution of reads.

python construct_read_lenght distribution_figure.py \
--in <FULL_PATH_TO_INPUT_FILE (FASTA)> \
--out <FULL_PATH_TO_OUTPUT_DIRECTORY>



13) tRNA gene annotation analysis: this script calculates the overlap between two different tRNA annotation methods. 

python compare_INFERNAL_to_tRNAscan.py \
--infernal <FULL_PATH_TO_INFERNAL_RESULT_FILE> \
--trnascan <FULL_PATH_TO_tRNAscan_RESULT_FILE>


14) extraction of sequence blocks from the assembly: this scripts allows the extraction of sequence blocks from FASTA file based on the sequence name and positional information.

python seqex.py \
--in <FULL_PATH_TO_INPUT_FILE> \
--out <FULL_PATH_TO_OUTPUT_FILE> \
--contig <STRING, name of contig> \
--start <INT, start of region to extract> \
--end <INT, end of region to extract>

