### Boas Pucker ###
### bpucker@cebitec.uni-bielefeld.de ###


import re, os
from operator import itemgetter

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
					sequences.update( { header: seq.upper() } )
					header = line.strip()[1:].split(' ')[0]
					seq = ""
			else:
				seq += line.strip()
			line = f.readline()
		sequences.update( { header: seq.upper() } )
	return sequences


def load_TE_positions( col_TE_file ):
	"""! @brief load TE positions from GFF file to construct query file """
	
	TE_fragments = []
	
	with open( col_TE_file, "r" ) as f:
		line = f.readline()
		while line:
			if line[0] != '#':
				parts = line.strip().split('\t')
				if parts[2] == "transposon_fragment":
					ID = re.findall( "AT\dTE\d+:transposon_fragment:\d+", line )[0]
					parent = re.findall( "AT\dTE\d+", line )[0]
					TE_fragments.append( { 'id': ID, 'parent': parent, 'chr': parts[0], 'start': int( parts[3] ), 'end': int( parts[4] ), 'orientation': parts[6] } )
				elif parts[2] == "exon":
					ID = re.findall( "AT\dG\d+:exon:\d+", line )[0]
					parent = re.findall( "AT\dG\d+", line )[0]
					TE_fragments.append( { 'id': ID, 'parent': parent, 'chr': parts[0], 'start': int( parts[3] ), 'end': int( parts[4] ), 'orientation': parts[6] } )
			line = f.readline()
	
	return sorted( TE_fragments, key=itemgetter( 'chr', 'start', 'end' ) )


def construct_query_file( TEs, ref_seq, query_file ):
	"""! @brief construction of query file """
	
	with open( query_file, "w" ) as out:
		for TE in TEs:
			if TE['orientation'] == '+':
				seq = ( ref_seq[ TE['chr'] ][ TE['start']-1:TE['end'] ] ).upper()
			else:
				seq = revcomp( ref_seq[ TE['chr'] ][ TE['start']-1:TE['end'] ] )
			
			out.write( '>' + "_%_".join( map( str, [ TE['id'], TE['parent'], TE['chr'], TE['start'], TE['end'] ] ) ) + '\n' + seq + '\n' )


def revcomp( seq ):
	"""! @brief construct reverse complement of sequence """
	
	new_seq = []
	
	bases = { 'A':'T', 'T':'A', 'C':'G', 'G':'C' }
	for nt in seq:
		try:
			new_seq.append( bases[nt] )
		except:
			new_seq.append( 'N' )
	return ''.join( new_seq[::-1] )


def load_best_blast_hit( blast_result_file ):
	"""! @brief load best blast hit per query """
	
	best_hits = {}
	
	with open( blast_result_file, "r" ) as f:
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			try:
				data = best_hits[ parts[0] ]
				if float( parts[-1] ) > data['score']:
					del best_hits[ parts[0] ]
					best_hits.update( { parts[0]: float( parts[-1] ) } )
			except:
				best_hits.update( { parts[0]: float( parts[-1] ) } )
			line = f.readline()
	return best_hits


def load_target_blast_results( target_blast_result_file, best_self_hit_per_query, score_ratio_cutoff ):
	"""! @brief load all target hits above the score cutoff """
	
	results_per_query_id = {}
	
	with open( target_blast_result_file, "r" ) as f:
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			try:
				if ( float( parts[-1] ) / best_self_hit_per_query[ parts[0] ] ) > score_ratio_cutoff:	#remove low quality hits due to repetetive sequences
					try:
						hits = results_per_query_id[ parts[0] ]
						del results_per_query_id[ parts[0] ]
						if int( parts[8] ) < int( parts[9] ):
							hits.append( { 'id': parts[0], 'chr': parts[1], 'start': int( parts[8] ), 'end': int( parts[9] ), 'orientation': '+', 'score': float( parts[-1] ) } )
						else:
							hits.append( { 'id': parts[0], 'chr': parts[1], 'start': int( parts[9] ), 'end': int( parts[8] ), 'orientation': '-', 'score': float( parts[-1] ) } )
						results_per_query_id.update( { parts[0]: hits } )
					except KeyError:
						if int( parts[8] ) < int( parts[9] ):
							results_per_query_id.update( { parts[0]: [ { 'id': parts[0], 'chr': parts[1], 'start': int( parts[8] ), 'end': int( parts[9] ), 'orientation': '+', 'score': float( parts[-1] ) } ] } )
						else:
							results_per_query_id.update( { parts[0]: [ { 'id': parts[0], 'chr': parts[1], 'start': int( parts[9] ), 'end': int( parts[8] ), 'orientation': '-', 'score': float( parts[-1] ) } ] } )
			except KeyError:
				print parts[0] + " - hit is missing in Col-0 but present in Nd-1 !!"
			line = f.readline()
	print "number of recovered fragments: " + str( len( results_per_query_id.keys() ) )
	return results_per_query_id


def get_all_non_overlapping_hits( results_per_query_id ):
	"""! @brief walk through all elements and identify all non-overlapping hits """
	
	final_hits = [ [], [], [], [], [] ]
	raw_hits = [ item for sublist in results_per_query_id.values() for item in sublist ]
	print "number of raw hits: " + str( len( raw_hits ) )
	t1 = datetime.now()
	t2 = datetime.now()
	for idx, hit in enumerate( raw_hits ):
		index = int( hit['chr'][-1] ) -1
		to_be_deleted = []
		status = True
		for k, each in enumerate( final_hits[ index ] ):
			if hit['start'] < each['end'] and hit['end'] > each['start']:
				if hit['score'] > each['score']:
					to_be_deleted.append( k )
				else:
					status = False
					break
		if status:
			for k in to_be_deleted[::-1]:	#start removing elements from the back to avoid index changes in front of this element
				del final_hits[ index ][ k ]
			final_hits[ index ].append( hit )
		
		if idx % 50000 == 0:
			print str( idx ) + "/" + str( len( raw_hits ) ) + "\t consumed time for last 10k elements: " + str( datetime.now()-t2 ) + "\ttotal time: " + str( datetime.now()-t1 ) 
			t2 = datetime.now()
	
	return final_hits


def construct_target_TE_GFF3_file( target_TE_file, non_overlapping_hits ):
	"""! @brief construct target GFF3 file containing mapped TE features """
	
	counter = 0
	with open( target_TE_file, "w" ) as out:
		out.write( "#gff-version\t3\n#Chr\tSource\tType\tStart\tEnd\tXXX\tStrand\tScore\tComment\n" )	#header lines
		for idx, chromosome in enumerate( non_overlapping_hits ):
			sorted_TEs = sorted( chromosome, key=itemgetter('start') )
			for run_num, TE in enumerate( sorted_TEs ):
				new_line = [ TE['chr'], "mapped_from_Araport11", "TE_fragment", TE['start'], TE['end'], ".", TE['orientation'], ".", "NdAT"+str(idx+1)+"TE"+(str(run_num).zfill(4))+"0;" + TE['id'] ]
				out.write( "\t".join( map( str, new_line ) ) + '\n' )
				counter += 1
	print "number of produced TE entries: " + str( counter )



if __name__ == '__main__':
	
	col_TE_file = "TAIR10_GFF3_transposons.gff3"
	col_ref_seq_file = "TAIR10.fa"
	
	nd1_ref_seq_file = "pseudochromosomes.fasta"
	
	prefix = "OUTPUT_DIRECTORY"
	
	TEs = load_TE_positions( col_TE_file )	#sorted by position on reference chromosome
	col_ref_seq = load_sequences( col_ref_seq_file )
	
	query_file = prefix + "Col0_TE_queryfile.fasta"
	construct_query_file( TEs, col_ref_seq, query_file )
	
	
	self_blast_db = prefix + "self_blast_db"
	target_blast_db = prefix + "target_blast_db"
	
	print "BLASTn ..."
	os.popen( "makeblastdb -in " + col_ref_seq_file + " -out " + self_blast_db + " -dbtype nucl" )
	os.popen( "makeblastdb -in " + nd1_ref_seq_file + " -out " + target_blast_db + " -dbtype nucl" )
	
	
	blastn_script = "/vol/cluster-data/bpucker/bin/scripts/run_blastn_on_cluster.py"
	tmp_result_file = query_file + "_final_BLAST_results.txt"
	
	self_blast_result_file = prefix + "self_blast_result_file.txt"
	target_blast_result_file = prefix + "target_blast_result_file.txt"
	
	#BLASTn vs. Col-0:
	print "... col"
	os.popen( "python blastn -query " + query_file + " -db " + self_blast_db + " -out " + self_blast_result_file+ " -outfmt 6 -evalue 0.0001" )
	
	
	#BLASTn vs. Nd-1:
	print "... nd1"
	os.popen( "python blastn -query_file " + query_file + " -db " + target_blast_db + " -out " + target_blast_result_file + " -outfmt 6 -evalue 0.0001" )
	
	
	print "processing BLASTn results ..."
	score_ratio_cutoff = 0.9
		
	best_self_hit_per_query = load_best_blast_hit( self_blast_result_file )
	target_blast_results = load_target_blast_results( target_blast_result_file, best_self_hit_per_query, score_ratio_cutoff )
	
	non_overlapping_hits = get_all_non_overlapping_hits( target_blast_results )
	
	print "producing GFF3 result file ... "
	target_TE_file = prefix + "Nd1_TEs.gff3"
	construct_target_TE_GFF3_file( target_TE_file, non_overlapping_hits )
	
	print "all done!"
