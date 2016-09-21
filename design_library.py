import argparse, tempfile, os, itertools, subprocess
import twobitreader
import sqlite3
import numpy
import time
import math
from operator import itemgetter
import csv

#guide design parameters
GUIDE_LENGTH = 20
PAM_LIST = ['AGG', 'TGG', 'GGG', 'CGG']
PAM_LENGTH = len(PAM_LIST[0])
CLEAVAGE_SITE = 17 #distance to 5' end of guide

#seqmap parameters
N_PROBES = 10000
MAX_PROCESSES = 1
MAX_MISMATCHES = 3
tf_counter = 0

#weights for off-target score calculations
weights = numpy.array([0,0,0.014,0,0,0.395,0.317,0,0.389,0.079,0.445,0.508,0.613,0.851,0.732,0.828,0.615,0.804,0.685,0.583])

#flanking sequences around spacer for gecko and sam libraries
gecko_flank = ['TTTCTTGGCTTTATATATCTTGTGGAAAGGACGAAACACCG', 'GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGT']
sam_flank = ['TTTCTTGGCTTTATATATCTTGTGGAAAGGACGAAACACCG', 'GTTTTAGAGCTAGGCCAACATGAGGATCACC']

def revcomp(sequence):
	"""
	returns the reverse complement of sequence
	"""
	basecomplement = {'A':'T', 'C':'G', 'T':'A', 'G':'C', 'N':'N'} 
	letters = list(sequence) 
	letters.reverse() 
	dna = ''
	for base in letters:
		dna += basecomplement[base] 
	return dna 


def indexList(s, item, i=0):
	"""
	make a list of indexes of item in s
	"""
	i_list = []
	while True:
		try:
			i = s.index(item, i)
			i_list.append(i)
			i += 1
		except:
			break
	return i_list

def Target_stretch(guide):
	"""
	returns true if guide does not contain any of the following homopolymer target stretches
	"""
	st1 = ('AAAA')
	st2 = ('TTTT')
	st3 = ('GGGG')
	st4 = ('CCCC')

	if not ((st1 in guide) or (st2 in guide) or (st3 in guide) or (st4 in guide)):
		return 'true'


def GC_content(GC_cutoff, guide):
	"""
	takes guide sequence as input, returns true if GC content above threshold defined above
	"""
	N = guide.count("G")
	N += guide.count("C")
	percent = float(N)/len(guide)*100
	if percent > GC_cutoff:
		return 'true'
	else: 
		return 'false'

def get_b_guides(region, GC_cutoff, index_list):
	"""
	takes a sequence s and and a list of indices in sequence that indicate the start of the
	reverse complement of the PAM sequence and returns a list of 20bp bottom guide sequences
	that have been filtered for GC content and target stretch
	"""
	guides = []
	for i in index_list:
		if len(region) > i + GUIDE_LENGTH + PAM_LENGTH:
			guide = (revcomp(region[i+PAM_LENGTH:i+PAM_LENGTH+GUIDE_LENGTH]))
			if 'N' not in guide:
				if  GC_content(GC_cutoff, guide) and Target_stretch(guide):
					guides.append([guide, i])
	return guides 

def get_t_guides(region, GC_cutoff, index_list):
	"""
	takes a sequence s and and a list of indices in sequence that indicate the start of the
	PAM sequence and returns a list of 20bp top guide sequences that have been filtered
	for GC content and target stretch
	"""
	guides = []
	for i in index_list:
		if i > GUIDE_LENGTH:
			guide = (region[i-GUIDE_LENGTH:i])
			if 'N' not in guide:
				if  GC_content(GC_cutoff, guide) and Target_stretch(guide):
					guides.append([guide, i])
	return guides 


def get_location(gene, b_guides, t_guides):
	"""
	returns lists of bottom guide and top guide cleavage site distances from the start of the
	target genomic region
	"""
	b_guide_loc = [(long(gene["start"]) + GUIDE_LENGTH + PAM_LENGTH - CLEAVAGE_SITE + x[1]) for x in b_guides]
	t_guide_loc = [(long(gene["start"]) - GUIDE_LENGTH + CLEAVAGE_SITE + x[1]) for x in t_guides]
	return b_guide_loc, t_guide_loc

def get_guides(region, GC_cutoff):
	"""
	finds all top and bottom guides in region and returns them with indices
	"""
	#find all the indices in region with pam sequence
	i_list_b = []
	i_list_t = []
	for pam in PAM_LIST:
		i_list_b.extend(indexList(region, revcomp(pam)))
		i_list_t.extend(indexList(region, pam))

	#find all the guides that correspond with pam indices
	b_guides = get_b_guides(region, GC_cutoff, i_list_b)  
	t_guides = get_t_guides(region, GC_cutoff, i_list_t)

	return b_guides, t_guides

def generate_fa_files(input_prefix, genes_file, GC_cutoff):
	"""
	generate a genome fa file from a genome 2bit file
	generate a guide fa file that contains filtered unique guides that target the genome
	at the regions specified in genes file
	"""
	genome_2bit_file = input_prefix + '.2bit'
	genome_fa_file = input_prefix + '.fa'
	guide_fa_file = input_prefix + '_all_guides.fa'
	tbf = twobitreader.TwoBitFile(genome_2bit_file)
	all_guides = set([])

	#iterate through the chromosomes in the genome 2bit file and write to genome fa file
	with open(genome_fa_file, 'wb') as genome_fa:
		for chrom in tbf:
			if '_' not in chrom:
				region = tbf[chrom][0:].upper()
				genome_fa.write('>{0}\n'.format(chrom))
				genome_fa.write(region+'\n')

	#for each genomic region specified in genes file identify filtered unique guides
	#and write to guide fa file
	with open(genes_file, 'rb') as gf:
		f = [row for row in csv.reader(gf.read().splitlines())]
		for i,l in enumerate(f): # i is index, l is entry
			if i == 0: 
				columns = l
				continue

			#fetch the current gene and region
			gene = dict([(columns[i],e) for i,e in enumerate(l)])
			region_bounds = [long(gene["start"]), long(gene["end"]) + 1]
			region = tbf[gene["chrom"]][region_bounds[0]:region_bounds[1]]
			region = region.upper()
			if "N" in region: 
				print "found N in target region of", gene["name"]
				continue

			#identify and filter guides that target region
			(b_guides, t_guides) = get_guides(region, GC_cutoff)
			current_guides = set([g[0] for g in (b_guides + t_guides)])
			all_guides = all_guides | current_guides

	#write filtered unique guides to an output guide_fa_file
	guide_count = 0
	with open(guide_fa_file, 'wb') as guide_fa:
		for guide in all_guides:
			guide_count += 1
			guide_fa.write('>{0}\n'.format(guide_count))
			guide_fa.write(guide+'\n')


def find_offtargets(input_prefix):
	"""
	calls seqmap to find all close matches to a given sgRNA listed in guide fa file in genome fa file
	and prints results to an offtargets file
	"""
	global tf_counter
	genome_fa_file = input_prefix + '.fa'
	guide_fa_file = input_prefix + '_all_guides.fa'
	offtargets_file = input_prefix + '_offtargets.tsv'

	#break up sgrnas by n_probes
	tempfiles_in = []
	with open(guide_fa_file) as guide_file_pointer:
		guide_file_pointer.seek(0)
		for k,g in itertools.groupby(enumerate(guide_file_pointer), 
									 key = lambda x:int(x[0]/(N_PROBES * 2))): 
			lines = list(e for i,e in g)
			if len(lines) == 0 : continue
			f_in = tempfile.NamedTemporaryFile(mode='w', suffix='.{0}.probes.input.fa'.format(tf_counter), prefix='temp.', delete=False)
			tf_counter +=1
			for line in lines: f_in.writelines(line)
			tempfiles_in.append(f_in)
			f_in.close()
	tempfiles_out = []

	#submit sgnra / probe scans {MAX_PROCESSES} at a time and wait for completion in groups
	for k,tfs_group in itertools.groupby(enumerate(tempfiles_in), key= lambda x:int(x[0] / MAX_PROCESSES)):
		f_out = tempfile.NamedTemporaryFile(mode='w', suffix='.{0}.seqmap.output'.format(tf_counter), prefix='temp.', delete=False)
		tf_counter +=1
		tempfiles_out.append(f_out)
		f_out.close()
		processes = []
	
		for k,f_in in tfs_group:
			print "PROCESSING IN PARALLEL!!!"
			cmd = "seqmap-1.0.13-src/seqmap {0} {1} {2} {3} /output_all_matches /do_not_output_probe_without_match".format(
				MAX_MISMATCHES, f_in.name, genome_fa_file, f_out.name)
			processes.append(subprocess.Popen(cmd, shell=True))
		for p in processes:
			p.communicate()

	#merge output and clean up
	for tf in tempfiles_in:
		os.remove(tf.name)
	with open(offtargets_file, "w") as offtargets_file_pointer:
		for tf in tempfiles_out:
			with open(tf.name) as f_out:
				lines = []
				for i,l in enumerate(f_out):
					if i == 0: continue
					lines.append(l)
					if l[-1] != '\n':
						print l
						raise Exception()
				print "sorting {0} offtarget hits".format(len(lines))
				lines_sorted = sorted(lines, key = lambda x:x.split('\t')[4])
				os.remove(tf.name)
				offtargets_file_pointer.writelines(lines_sorted)        
	return

def make_db(input_prefix):
	"""
	Creates a sqlite database of guides and offtarget scores based on the offtargets file
	and outputs to database file
	"""
	offtargets_file = input_prefix + '_offtargets.tsv'
	offtarget_scores_file = input_prefix + '_offtarget_scores.csv'
	database_file = input_prefix + '_database.sqlite'
	to_db =[]
	ot_number = 0

	with open(offtargets_file) as otf:
		for ontarget_sequence,g in itertools.groupby(otf, key = lambda x:x.split("\t")[4]):
			rows = list(g)
			offtarget_sequences = [e.split('\t')[2] for e in rows]

			if len(ontarget_sequence) != GUIDE_LENGTH:
				print rows[0].split('\t')
				raise Exception("improperly formatted SGRNA")
			for i,s in enumerate(offtarget_sequences):
				if len(s) != GUIDE_LENGTH: 
					raise Exception("improperly formatted OT")

			scores = [score_one_offtarget(ontarget_sequence, e) for e in offtarget_sequences]
			ot_number += 1
			if ot_number%10000 == 0:
				print ot_number, "guides scored"
				if ot_number == 10000:
					tic = time.clock()
				elif ot_number ==20000:
					toc = time.clock()
					print "time to process 10000 guides:", toc - tic
			if len(scores) == 0: raise Exception('ERROR: no matches for target.')
			total_score = score_sgrna(scores, has_ontarget = True)
			to_db.append((ontarget_sequence, total_score))

		# make a csv file with unique isoforms  
		os.system('touch ' + offtarget_scores_file)
		with open(offtarget_scores_file,'wb') as csvfile:
			mywriter = csv.writer(csvfile)
			for guide in to_db:
				mywriter.writerow(guide)  

	print "opening sqlite connection to {0}".format(database_file)
	dbpath = os.path.abspath(database_file)
	if os.path.isfile(database_file):
		os.remove(database_file)
	dbaddress = "//{0}".format(dbpath)
	print dbaddress
	con = sqlite3.connect(dbaddress)
	cur = con.cursor()
	cur.execute("CREATE TABLE sgrnas (seq TEXT PRIMARY KEY, score NUM);")          
	cur.executemany("INSERT INTO sgrnas (seq, score) VALUES (?, ?);", to_db)
	print "committing"
	con.commit()
	print "done creating databse with {0} entries".format(len(to_db))
	
def score_one_offtarget(sgrna_sequence, offtarget_sequence):
	"""
	scores a single offtarget match
	"""
	mismatches = numpy.array([i for i in range(GUIDE_LENGTH) 
				  if sgrna_sequence[i] != offtarget_sequence[i]])
	if len(mismatches) == 0:
		score = 100
	else:
		score = 100 * (1 - weights[mismatches]).prod()
		if len(mismatches) > 1:
			mean_pairwise =float(sum(mismatches[1:] - mismatches[:-1])) / (len(mismatches)-1)
			mpw_factor = ((float((19-mean_pairwise))/19)*4 + 1)
			scl_factor = pow(len(mismatches),2)
			score  = score / ( mpw_factor * scl_factor )
			score = max([score,0])
	return score

def score_sgrna(scores, has_ontarget = True):
	"""
	computes a total score for an sgRNA guide sequence from all offtargets
	"""
	sum_scores = float(sum(scores))
	norm_score = 100 / sum_scores
	return norm_score

def remove_overlap(all_guides_sorted, spacing):
	"""
	removes guides that have cleavage site distances that are less than specified by spacing
	and returns a filtered list of guides
	"""
	all_guides_filtered = [all_guides_sorted[0]]
	prev_guide = all_guides_sorted[0]
	prev_loc = prev_guide[-1]

	for guide in all_guides_sorted[1:]:
		loc = guide[-1]
		if abs(loc - prev_loc) > spacing:
			all_guides_filtered.append(guide)
			prev_guide = guide
			prev_loc = loc

	return all_guides_filtered

def get_sorted_guides(region, gene, GC_cutoff, spacing, input_prefix):
	"""
	returns a list of filtered guides in region sorted by distance to the start of the
	targeted region in the form:
	[name, spacer sequence, strand (b/t), chromosome, and cleavage site location]
	"""
	all_b_guides = []
	all_t_guides = []
	(b_guides, t_guides) = get_guides(region, GC_cutoff)
	#add location of cleavage site to each guide
	(b_guide_loc, t_guide_loc) = get_location(gene, b_guides, t_guides)
	for i, loc in enumerate(b_guide_loc): 
		all_b_guides.append([gene["name"], b_guides[i][0], "b", gene["chrom"], loc]) 
	for i, loc in enumerate(t_guide_loc):
		all_t_guides.append([gene["name"], t_guides[i][0], "t", gene["chrom"], loc])
	all_guides = all_b_guides + all_t_guides #makes one nested list of bottom and top guides
	all_guides_sorted = sorted(all_guides, key=itemgetter(-1))  #sorts by location
	all_guides_filtered = remove_overlap(all_guides_sorted, spacing)

	return all_guides_filtered

def get_ot_guides(guide_list, input_prefix):
	"""
	connects to the offtarget database to fetch ot scores for guides in guide_list
	returns a list of guides with respective off-target scores
	"""
	database_file = input_prefix + '_database.sqlite'
	guides_scored = []
	sequences = []

	for guide in guide_list:
		sequences.append(guide[1])

	dbpath = os.path.abspath(database_file)
	dbaddress = "//{0}".format(dbpath)
	conn = sqlite3.connect(dbaddress)
	c = conn.cursor()
	query = "select * FROM sgrnas WHERE seq IN ({0}) ORDER BY score".format(','.join(['?']*len(sequences)))
	try:
		c.execute(query, sequences)
		guides_scored = c.fetchall() 
	   
	except:
		print "did not find guides"

	return guides_scored

def list_sgrnas(genes_file, input_prefix, GC_cutoff, spacing, guides_per_gene, gecko, sam):
	"""
	Returns a list of (ontarget) sgrna sequences using a genome file and list of
	transcription start sites form a .csv file.
	"""
	genome_2bit_file = input_prefix + '.2bit'
	tbf = twobitreader.TwoBitFile(genome_2bit_file)
	final_guides = []

	with open(genes_file, 'rb') as gf:
		f = [row for row in csv.reader(gf.read().splitlines())]
		for i,l in enumerate(f):
			if i == 0: 
				columns = l
				continue

			#fetch the current gene and region
			gene = dict([(columns[i],e) for i,e in enumerate(l)])
			region_bounds = [long(gene["start"]), long(gene["end"]) + 1]
			region = tbf[gene["chrom"]][region_bounds[0]:region_bounds[1]]
			region = region.upper()
			if "N" in region: 
				print "found N in target region of", gene["name"]
				continue

			#identify and filter guides that target region
			guides = get_sorted_guides(region, gene, GC_cutoff, spacing, input_prefix)

			#add offtarget scores to filtered guides and select guides with higher offtarget scores
			ot_guides_sql = get_ot_guides(guides, input_prefix) 
			ot_guides_dict = dict(ot_guides_sql)
			for g in guides:
				spacer = g[1]
				g.append(ot_guides_dict[spacer])

			#sort and add guides with the highest offtarget scores to final guides
			guides = sorted(guides, key=itemgetter(-1), reverse=True)

			if len(guides) <= guides_per_gene:
				final_guides.extend(guides)
			else:
				final_guides.extend(guides[:guides_per_gene])

	# add gecko or sam flanking sequences to the spacer for the oligo library
	if sam or gecko:
		for guide in final_guides:
			spacer = guide[1]
			if gecko:
				oligo = gecko_flank[0] + spacer + gecko_flank[1]
			if sam:
				oligo = sam_flank[0] + spacer + sam_flank[1]
			guide.append(oligo)

	return final_guides

def writecsv(data, filename):
	"""
	write data to a csv file
	"""
	with open(filename, 'wb') as csvfile:
		csvwriter = csv.writer(csvfile)
		for row in data:
			csvwriter.writerow(row)      

def __main__():
	parser = argparse.ArgumentParser(
		description='Design oligo library sequences for custom library cloning')
	parser.add_argument('-o', '--output', type=str, dest='guide_file',
						help='output file name', default='final_guides.csv')
	parser.add_argument('-i', '--input', type=str, dest='input_prefix',
						help='input genome prefix', default='hg19')
	parser.add_argument('-g', '--genes', type=str, dest='genes_file',
						help='input gene file name', default='genes.csv')
	parser.add_argument('-gc', '--gc', type=int, dest='GC_cutoff',
						help='gc content cutoff', default=25)
	parser.add_argument('-s', '--spacing', type=int, dest='spacing',
						help='minimum spacing between cleavage sites', default=20)
	parser.add_argument('-n', '--guides-per-gene', type=int, dest='guides_per_gene',
						help='maximum number of guides per gene', default=3)
	parser.add_argument('-db', dest='db', help='use existing off-target database', action='store_false')
	parser.set_defaults(db=True)
	parser.add_argument('-gecko', dest='gecko', help='add gecko flanking sequences', action='store_true')
	parser.set_defaults(gecko=False)
	parser.add_argument('-sam', dest='sam', help='add sam flanking sequences', action='store_true')
	parser.set_defaults(sam=False)
	args = parser.parse_args()

	if args.db:
		generate_fa_files(args.input_prefix, args.genes_file, args.GC_cutoff)
		find_offtargets(args.input_prefix)
		make_db(args.input_prefix)  
	final_guides = list_sgrnas(args.genes_file, args.input_prefix, args.GC_cutoff, args.spacing, args.guides_per_gene, args.gecko, args.sam)
	writecsv(final_guides, args.guide_file)

if __name__ == "__main__":
	__main__()

