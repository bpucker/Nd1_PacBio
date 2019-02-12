### Boas Pucker ###
### bpucker@cebitec.uni-bielefeld.de ###
### v0.1 ###

__usage__ = """
					python compare_INFERNAL_to_tRNAscan.py
					--infernal <FULL_PATH_TO_INFERNAL_RESULT_FILE>
					--trnascan <FULL_PATH_TO_tRNAscan_RESULT_FILE>
					"""


import sys

# --- end of imports --- #

def load_tRNAscan_results( filename ):
	"""! @brief load predicted tRNAs from tRNAscan-SE result file """
	
	tRNAs = []
	with open( filename, "r" ) as f:
		f.readline()	#header1
		f.readline()	#header2
		f.readline()	#header3
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			tRNAs.append( { 'chr': parts[0].strip(), 'start': int( parts[2].strip() ), 'end': int( parts[3].strip() ) } )
			line = f.readline()
	return tRNAs


def load_infernal_results( infernal_file ):
	"""! @brief load predicted RNA genes from given input file """
	
	infernal = []
	with open( infernal_file, "r" ) as f:
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			if "tRNA" in parts[-1]:
				infernal.append( { 'chr': parts[0], 'start': int( parts[3] ), 'end': int( parts[4] ) } )
			line = f.readline()
	return infernal


def main( arguments ):
	"""! @brief runs all parts of this script """
	
	infernal_file = arguments[ arguments.index( '--infernal' )+1 ]
	tRNAscan_file = arguments[ arguments.index( '--trnascan' )+1 ]
	
	tRNAs = load_tRNAscan_results( tRNAscan_file )
	print "number of tRNAs in tRNAscan-SE: " + str( len( tRNAs ) )
	infernal = load_infernal_results( infernal_file )
	print "number of tRNAs in INFERNAL: " + str( len( infernal ) )
	
	counter = 0
	for gene1 in tRNAs:
		for gene2 in infernal:
			if gene1['chr'] == gene2['chr']:
				if gene1['start'] <= gene2['end']:
					if gene1['end'] >= gene2['start']:
						counter += 1
	print "number of overlapping annotations: " + str( counter )


if '--infernal' in sys.argv and '--trnascan' in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )
