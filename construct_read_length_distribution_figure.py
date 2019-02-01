### Boas Pucker ###
### bpucker@cebitec.uni-bielefeld.de ###

__usage__ = """
					python construct_read_lenght distribution_figure.py
					--in <FULL_PATH_TO_INPUT_FILE (FASTA)>
					--out <FULL_PATH_TO_OUTPUT_DIRECTORY>
					"""

import matplotlib.pyplot as plt
import sys

# ---- end of imports --- #

def load_readlen( input_file ):
	"""! @brief get lengthes of all sequences in given file """
	
	seq_lens = []
	
	with open( input_file, "r" ) as f:
		line = f.readline()
		seq_len = 0
		while line:
			if line[0] == ">":
				seq_lens.append( seq_len )
				seq_len = 0
			else:
				seq_len += len( line.strip() )
			line = f.readline()
	return seq_lens


def main( arguments ):
	"""! @brief runs everything """
	
	input_file = arguments[ arguments.index( '--in' )+1 ]
	seq_lens = load_readlen( input_file )
	
	prefix = arguments[ arguments.index( '--out' )+1 ]
	
	if prefix[-1] != "/":
		prefix += "/"
	if not os.path.exists( prefix ):
		os.makedir( prefix )
	
	tmp_file = prefix + "tmp.txt"
	with open( tmp_file, "w" ) as out:
		out.write( "\n".join( map( str, seq_lens ) ) )
	
	with open( tmp_file, "r" ) as f:
		content = f.read()
	
	seq_lens = map( int, content.strip().split('\n') )
	
	print "max seq length: " + str( max( seq_lens ) )
	
	fig_file = prefix + "read_len_distribution.png"
	fig, ax = plt.subplots()
	ax.hist( seq_lens, bins=1000 )
	ax.set_xlim( [ 0, 70000 ] )
	ax.set_xlabel( "read length [nt]" )
	ax.set_ylabel( "number of reads" )
	fig.savefig( fig_file, dpi=600 )
	
	print "all done!"


if __name__ == '__main__':
	if '--in' in sys.argv and '--out' in sys.argv:
		main( sys.argv )
	else:
		sys.exit( __usage__ )
	
