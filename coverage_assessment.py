### Boas Pucker ###
### bpucker@cebitec.uni-bielefeld.de ###

import re, sys
import numpy as np

# --- end of imports --- #


def get_all_gene_positions( gff3_file ):
	"""! @brief get all gene positions from GFF3 file """
	
	gene_positions = {}
	with open( gff3_file, "r" ) as f:
		line = f.readline()
		while line:
			if line[0] != '#':
				parts = line.strip().split('\t')
				if parts[2] == "gene":
					ID = re.findall( "NdCChr\d\.g\d+", parts[-1] )[0]
					gene_positions.update( { ID: { 'id': ID, 'chr': parts[0], 'start': int( parts[3] ), 'end': int( parts[4] ) }} )
			line = f.readline()
	return gene_positions


def load_cov( cov_file ):
	"""! @brief load coverage values """
	
	coverage = {}
	with open( cov_file, "r" ) as f:
		line = f.readline()
		chromosome = line.split('\t')[0]
		cov = []
		while line:
			parts = line.strip().split('\t')
			if parts[0] != chromosome:
				coverage.update( { chromosome: cov } )
				chromosome = parts[0]
				cov = []
			cov.append( float( parts[-1] ) )
			line  = f.readline()
		coverage.update( { chromosome: cov } )
	return coverage


def load_genes( gene_of_interest_file ):
	"""! @brief load genes of interest """
	
	genes = []
	with open( gene_of_interest_file, "r" ) as f:
		f.readline()
		line = f.readline()
		while line:
			genes.append( line.strip() )
			line = f.readline()
	return genes


def main( arguments ):
	"""! @brief run all parts of this script """
	
	gff_file = "AthNd1_v2c.gff3"
	cov_file = arguments[ arguments.index('--cov')+1 ]
	#report_file = cov_file.replace( ".cov", "_report.txt" )
	report_file = cov_file.replace( ".cov", "_report.txt" )
	
	gene_positions = get_all_gene_positions( gff_file )
	coverage = load_cov( cov_file )

	genes = gene_positions.keys() #load_genes( gene_of_interest_file )

	avg_cov = np.median( [ each for sublist in coverage.values() for each in sublist ] )
	print avg_cov
	
	with open( report_file, "w" ) as out:
		for gene in sorted( genes ):
			info = gene_positions[ gene ]
			cov = np.median( coverage[ info['chr'] ][ info['start']:info['end'] ] )
			out.write( "\t".join( map( str, [ gene, avg_cov, cov ] ) ) + '\n' )


if '--cov' in sys.argv:
	main( sys.argv )
