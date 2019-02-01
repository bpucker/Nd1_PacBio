### Boas Pucker ###
### bpucker@cebitec.uni-bielefeld.de ###

import matplotlib.pyplot as plt
from scipy import stats

# --- end of imports --- #


def load_sequences( multiple_fasta_file ):
	"""! @brief load candidate gene IDs from file """
	
	sequences = {}
	
	with open( multiple_fasta_file ) as f:
		header = f.readline()[1:].strip().split(' ')[0]
		seq = ""
		line = f.readline()
		while line:
			if line[0] == '>':
					sequences.update( { header: seq } )
					header = line.strip()[1:].split(' ')[0]
					seq = ""
			else:
				seq += line.strip()
			line = f.readline()
		sequences.update( { header: seq } )
	return sequences


def chunks( seq, n):
    """! @brief yield successive n-sized chunks from seq """
    
    for i in range( 0, len( seq ), n ):
		if i+n < len( seq ):
			yield seq[ i : i + n ]
		else:
			yield seq[ i: ]


def construct_plot_per_seq( seq, outputfile ):
	"""! @brief construct plot for sequence composition """
	
	chunk_size = 100
	seq = seq.upper()
	blocks = chunks( seq, chunk_size)
	x_pos = []
	a_y = []
	c_y = []
	g_y = []
	t_y = []
	n_y = []
	
	for idx, block in enumerate( blocks ):
		x_pos.append( idx*chunk_size + 0.5*chunk_size )
		a_y.append( 100*block.count('A') / float( chunk_size ) )
		c_y.append( 100*block.count('C')  / float( chunk_size ))
		g_y.append( 100*block.count('G')  / float( chunk_size ))
		t_y.append( 100*block.count('T') / float( chunk_size ) )
		n_y.append( 100*block.count('N') / float( chunk_size ))
	
	fig, ax = plt.subplots()
	ax.scatter( x_pos, a_y, color="green", label="A" )
	ax.scatter( x_pos, c_y, color="blue", label="C" )
	ax.scatter( x_pos, g_y, color="black", label="G" )
	ax.scatter( x_pos, t_y, color="red", label="T" )
	ax.scatter( x_pos, n_y, color="purple", label="N" )
	
	ax.set_ylim( 0, 100 )
	ax.set_xlim( 0, max( x_pos ) )
	
	#plt.show()
	fig.savefig( output_file, dpi=300 )
	plt.close( "all" )


if __name__ == '__main__':
	
	input_file1 = "gaps.fasta"
	
	prefix1 = "OUTPUT_DIR1"
	
	input_file2 = "control_seqs.fasta"
	prefix2 = "OUTPUT_DIR2"
	
	final_figure = "summary.png"
	
	seqs1 = load_sequences( input_file1 )
	
	lengths1 = []
	frequencies1 = []
	
	for key in sorted( seqs1.keys() ):
		seq = seqs1[ key ]
		print "length: " + str( len( seq ) )
		print "AAAA: " + str( seq.count( "AAAA" ) / ( len( seq )*0.001 ) )
		print "CCCC: " + str( seq.count( "CCCC" ) / ( len( seq )*0.001 ) )
		print "GGGG: " + str( seq.count( "GGGG" ) / ( len( seq )*0.001 ))
		print "TTTT: " + str( seq.count( "TTTT" ) / ( len( seq )*0.001 ))
		value = (seq.count( "AAAA" )+seq.count( "CCCC" )+seq.count( "GGGG" )+seq.count( "TTTT" ) ) / ( len( seq )*0.001 )
		frequencies1.append( value )
		lengths1.append( len( seq )*0.001 )
		print "tetra homopolymer frequency: " + str( value )
		print ""
		output_file = prefix1 + key + ".png"
		construct_plot_per_seq( seq, output_file )
	
	seqs2 = load_sequences( input_file2 )
	
	lengths2 = []
	frequencies2 = []
	for key in sorted( seqs2.keys() ):
		seq = seqs2[ key ]
		print "length: " + str( len( seq ) )
		value = (seq.count( "AAAA" )+seq.count( "CCCC" )+seq.count( "GGGG" )+seq.count( "TTTT" ) ) / ( len( seq )*0.001 )
		frequencies2.append( value )
		lengths2.append( len( seq )*0.001 )
		print "tetra homopolymer frequency: " + str( value )
		print ""
		output_file = prefix2 + key + ".png"
		construct_plot_per_seq( seq, output_file )
	
	
	fig, ax = plt.subplots()
	ax.scatter( lengths1, frequencies1, color="red", label="gaps" )
	ax.scatter( lengths2, frequencies2, color="green", label="control" )
	
	ax.legend()
	
	ax.set_xlabel( "gap length [kb]" )
	ax.set_ylabel( "homotranucleotides per kb" )
	ax.set_xlim( 0, max( lengths1+lengths2 ) )
	ax.set_ylim( 0, max( frequencies1+frequencies2 ) )
	
	fig.savefig( final_figure, dpi=900 )
	#plt.show()
	
	print stats.mannwhitneyu( frequencies1, frequencies2 )
	
	
	print "all done!"
