import csv
import sys
import argparse

# flanking sequences around spacer for gecko and sam libraries
gecko_flank = ['TTTCTTGGCTTTATATATCTTGTGGAAAGGACGAAACACCG', 'GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGT']
sam_flank = ['TTTCTTGGCTTTATATATCTTGTGGAAAGGACGAAACACCG', 'GTTTTAGAGCTAGGCCAACATGAGGATCACC']

def design_oligos(output_file, library_file, genes_file, gecko, sam):
	"""
	creates a list of spacers corresponding to the target genes, writes to output_file
	output_file: subset of spacers corresponding to the target genes
	library_file: all RefSeq genes with corresponding spacers
	genes_file: target genes
	"""
	# list of spacer sequences for the library
	spacer_list = []

	# open gene file for list of target genes
	with open(genes_file, 'rb') as infile:
		target_genes = [row[0] for row in csv.reader(infile.read().splitlines())]

	# open library file for all RefSeq genes
	with open(library_file, 'rb') as infile:
		for row in csv.reader(infile.read().splitlines()):
			gene = row[0]
			spacer = row[1]

			# check if each gene is in the list of target genes
			if gene in target_genes:
				spacer_list.append(row)

				# add gecko or sam flanking sequences to the spacer for the oligo library
				if gecko:
					oligo = gecko_flank[0] + spacer + gecko_flank[1]
					spacer_list[-1].append(oligo)
				if sam:
					oligo = sam_flank[0] + spacer + sam_flank[1]
					spacer_list[-1].append(oligo)
		# sort spacer sequences in gene order			
		spacer_list = sorted(spacer_list, key=lambda t: t[0])

	# write spacer list to output file
	with open(output_file, 'wb') as outfile:
		csvwriter = csv.writer(outfile)
		for s in spacer_list:
			csvwriter.writerow(s)

if __name__ == '__main__':
	parser = argparse.ArgumentParser(
		description='Design oligo library sequences for targeted library cloning')
	parser.add_argument('-o', '--output', type=str, dest='output_file',
						help='output file name', default='oligos.csv')
	parser.add_argument('-l', '--library', type=str, dest='library_file',
						help='input file name', default='annotated_library.csv')
	parser.add_argument('-g', '--genes', type=str, dest='genes_file',
					help='input file name', default='target_genes.csv')
	parser.add_argument('-gecko', dest='gecko', help='add gecko flanking sequences', action='store_true')
	parser.set_defaults(gecko=False)
	parser.add_argument('-sam', dest='sam', help='add sam flanking sequences', action='store_true')
	parser.set_defaults(sam=False)
	args = parser.parse_args()

	design_oligos(args.output_file, args.library_file, args.genes_file, args.gecko, args.sam)
