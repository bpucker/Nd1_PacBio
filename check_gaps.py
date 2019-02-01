### Boas Pucker ###
### bpucker@cebitec.uni-bielefeld.de ###

import os

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


def load_gap_positions( gap_file ):
	"""! @brief load all gap infos """
	
	gaps = []
	
	with open( gap_file, "r" ) as f:
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			start, end = sorted( map( int, parts[3:] ) )
			gaps.append( { 'chr': parts[0], 'start': start, 'end': end } )
			line = f.readline()
	return gaps


def construct_query_file( gaps, col_seqs, query_file, flank_len ):
	"""! @brief construct file with gap flanking sequences """
	
	counter = 0
	with open(  query_file, "w" ) as out:
		for idx, gap in enumerate( gaps ):
			if gap['start'] > flank_len and len( col_seqs[ gap['chr'] ] )-gap['end'] > flank_len:
				out.write( '>gap'  + str( idx ).zfill(2) + '_1\n' + col_seqs[ gap['chr'] ][ gap['start']-flank_len:gap['start'] ] + '\n' )
				out.write( '>gap'  + str( idx ).zfill(2) + '_2\n' + col_seqs[ gap['chr'] ][ gap['end']:gap['end']+flank_len ] + '\n' )
				counter += 1
	print "number of query sequences: " + str( counter )


def load_best_blast_hit_per_query( blast_result_file ):
	"""! @brief load best BLAST hit """
	
	best_hits = {}
	with open( blast_result_file, "r" ) as f:
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			try:
				hit = best_hits[ parts[0] ]
				if float( parts[-1] ) > hit['score']:
					del best_hits[ parts[0] ]
					best_hits.update( { parts[0]: { 'chr': parts[1], 'start': int( parts[8] ), 'end': int( parts[9] ), 'score': float( parts[-1] ) } } )
			except KeyError:
				best_hits.update( { parts[0]: { 'chr': parts[1], 'start': int( parts[8] ), 'end': int( parts[9] ), 'score': float( parts[-1] ) } } )
			line = f.readline()
	return best_hits


def identify_matching_seqs( best_hits, gaps, new_genome_seqs ):
	"""! @brief get matching sequences for all given gaps """
	
	gap_results = []
	
	for idx, gap in enumerate( gaps ):
		try:
			hit1 = best_hits[ 'gap' + str( idx ).zfill(2) + "_1" ]
			hit2 = best_hits[ 'gap' + str( idx ).zfill(2) + "_2" ]
			
			if hit1['chr'] == hit2['chr']:
				if hit1['end'] < hit2['start']:
					seq_of_interest = new_genome_seqs[ hit1['chr'] ][ hit1['end']:hit2['start'] ]
					start = hit1['end']
					end = hit2['start']
				else:
					end, start = sorted( [ hit2['start'], hit1['end'], hit2['end'], hit1['start'] ] )[1:3]
					seq_of_interest = new_genome_seqs[ hit1['chr'] ][ start:end ]
				
				if len( seq_of_interest ) > 1 and len( seq_of_interest )<800000:
					gap_results.append( { 'id': idx, 'seq': seq_of_interest, 'chr': hit1['chr'], 'start': start, 'end': end  } )
			else:
				print "hits on different sequences! - stopping analysis."
		except KeyError:
			print gap
	return gap_results


if __name__ == '__main__':
	
	gap_file = "TAIR10.gaps.txt"
	TAIR10_file = "TAIR10.fa"
	
	new_genome_seq_file = "AthNd1_v2c.fasta"
	
	prefix = "OUTPUT_DIRECTORY"
	
	
	# --- load gap position data --- #
	gaps = load_gap_positions( gap_file )
	col_seqs = load_sequences( TAIR10_file )
	new_genome_seqs = load_sequences( new_genome_seq_file )
	
	# --- construct_query_file --- #
	query_file = prefix + "query.fasta"
	construct_query_file( gaps, col_seqs, query_file, flank_len=30000 )
	
	
	# --- run BLAST --- #
	blast_result_file = prefix + "blast_result_file.txt"
	blast_db = prefix + "blastdb"
	os.popen( "makeblastdb -in " + new_genome_seq_file + " -out " +  blast_db + " -dbtype nucl" )
	os.popen( "blastn -query " + query_file + " -db " + blast_db + " -out " + blast_result_file + " -outfmt 6 -evalue 0.000000000000000000000000000000000000000000000001 -num_threads 16" )
	
	# --- get gap matching sequences --- #
	best_hits = load_best_blast_hit_per_query( blast_result_file )
	gap_data = identify_matching_seqs( best_hits, gaps, new_genome_seqs )
	
	# --- save gap matching sequences --- #
	print "number of closed gaps: " + str( len( gap_data ) ) + "(total gaps: " + str( len( gaps ) ) + ")"
	result_file = prefix + "results.fasta"
	with open( result_file, "w" ) as out:
		closed_gaps = []
		for gap in gap_data:
			if not 'N' in gap['seq']:
				status = True
				for each in closed_gaps:
					if each['chr'] == gap['chr']:
						if each['start'] < gap['end']:
							if each['end'] > gap['start']:
								status = False
				if status:
					closed_gaps.append( gap )
					out.write( '>gap' + str( gap['id']  ).zfill( 2 ) + '_' + gap['chr'] + '_' + str( gap['start'] ) + "_" + str( gap['end'] ) + '\n' + gap['seq'] + '\n' )
	
	print "all done!"
