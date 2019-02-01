### Boas Pucker ###
### bpucker@cebitec.uni-bielefeld.de ###

import re
import matplotlib.pyplot as plt

# --- end of imports --- #

def load_all_gene_positions( gene_gff3_file ):
	
	gene_positions = []
	
	with open( gene_gff3_file, "r" ) as f:
		line = f.readline()
		while line:
			if line[0] != '#':
				parts = line.strip().split('\t')
				if parts[2] == "gene":
					ID = re.findall( "Nd[CFMY]*Chr\d\.g\d+", line )[0]
					gene_positions.append( { 'id': ID, 'chr': parts[0], 'start': int( parts[3] ), 'end': int( parts[4] ) } )
			line = f.readline()
	return gene_positions


def load_all_TE_positions( TE_gff3_file ):
	
	TE_positions = []
	
	with open( TE_gff3_file, "r" ) as f:
		line = f.readline()
		while line:
			if line[0] != '#':
				parts = line.strip().split('\t')
				try:
					ID = re.findall( "AT\dG\d+", line )[0]
				except:
					ID = re.findall( "AT\dTE\d+", line )[0]
				TE_positions.append( { 'id': ID, 'chr': parts[0], 'start': int( parts[3] ), 'end': int( parts[4] ) } )
			line = f.readline()
	return TE_positions


def get_TE_genes( te_positions, gene_positions ):
	
	TE_effected_genes = []
	TE_overlap_lens = []
	TE_overlap_fraction = []
	
	for idx, gene in enumerate( gene_positions ):
		overlap_len = 0
		for TE in te_positions:
			if gene['chr'] == TE['chr']:
				if gene['start'] < TE['end']:
					if gene['end'] > TE['start']:
						m1, m2 = sorted( [ gene['start'], TE['end'], gene['end'], TE['start'] ] )[1:3]
						overlap_len += m2-m1
		if overlap_len > 0:
			TE_overlap_lens.append( overlap_len )
			TE_effected_genes.append( gene['id'] )
			TE_overlap_fraction.append( float( overlap_len ) / ( gene['end']-gene['start'] ) )
		if idx % 500 == 0:
			print str( idx ) + "/" + str( len( gene_positions ) )
	return TE_effected_genes, TE_overlap_lens, TE_overlap_fraction


if __name__ == '__main__':
	
	nd1_te_gff3_file = "Nd1_TEs.gff3"
	nd1_genes_gff3_file = "AthNd1_v2m.gff3"
	
	nd1_te_positions = load_all_TE_positions( nd1_te_gff3_file )
	print "number of loaded TE positions: " + str( len( nd1_te_positions ) )
	nd1_gene_positions = load_all_gene_positions( nd1_genes_gff3_file )
	print "number of loaded gene positions: " + str( len( nd1_gene_positions ) )
	
	TE_effected_genes, TE_overlap_lens, TE_overlap_fraction = get_TE_genes( nd1_te_positions, nd1_gene_positions )
	print "number of TE effected genes: " + str( len( TE_effected_genes ) )
	
	plt.hist( TE_overlap_fraction, bins=100 )
	plt.title( "Overlapping fraction of TEs and genes" )
	plt.xlabel( "fraction of overlap" )
	plt.ylabel( "number of genes" )
	plt.xlim( [0, 1] )
	plt.savefig( "Nd1_TE_genes.png", dpi=600 )
	#plt.show()
		
	cutoff = 0.8
	
	final_TE_genes = []
	for idx, gene in enumerate( TE_effected_genes ):
		if TE_overlap_fraction[ idx ] > cutoff:
			final_TE_genes.append( gene )
		
	TE_gene_file = "Nd1_TE_genes.txt"
	with open( TE_gene_file, "w" ) as out:
		out.write( "\n".join( final_TE_genes ) )
		
	print "number of final TE genes: " + str( len( final_TE_genes ) )

