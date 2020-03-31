
import time
from docopt import docopt
from collections import defaultdict
from Bio import SeqIO
import os
import sys


def blast_search(toxinNT, toxinPro, antitoxinNT, dataPath, outPath, threads):

	for index, file in enumerate([toxinNT, toxinPro, antitoxinNT]):
		if not os.path.isfile(file):
			# if index == 2: # if no at
			# 	sys.exit()
			continue
		else:
			if index == 0:
				blast_command('blastn', file, f'{dataPath}Homo_sapiens.GRCh38.dna.primary_assembly.fa', f'{outPath}t_nt_blast.txt', threads)
			elif index == 2:
				blast_command('blastn', file, f'{dataPath}Homo_sapiens.GRCh38.dna.primary_assembly.fa', f'{outPath}at_nt_blast.txt', threads)
			elif index == 1:
				blast_command('blastn', file, f'{dataPath}Homo_sapiens.GRCh38.dna.primary_assembly.fa', f'{outPath}t_pro_blast.txt', threads)

def blast_command(blastType, tasInput, databaseFile, out, threads):
	if os.path.isfile(databaseFile):
		os.system(f'nice {blastType} -num_threads {threads} -query {tasInput} -db {databaseFile} -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen gaps stitle sstrand" -out {out}') #add qlen to outputformat, to get query length

#docopt for command line arguments 	

args = docopt("""blast_search.

Usage:
    box_1_blast_search.py [-h] [--tnt FILE] [--tpro FILE] [--atnt FILE] [--t INT] -o DIR -d FILE 

Options:
    -h --help    show this message and exit.
    -o DIR       Dir for output
    -d FILE      genome database
    --tnt FILE   FASTA file with one or more toxin sequences (nucleotide)
    --tpro FILE  FASTA file with one or more toxin sequences (protein)
    --atnt FILE  FASTA file with one ot more antitoxin sequences (nucleotide)
    --t INT      number of reminders to give. [default: 1]

""")

if __name__ == '__main__':
	print(args)
	bacDatabase = args["-d"]
	savePath = args["-o"]

	tNtInfile = args['--tnt']
	tProInfile = args['--tpro']
	atInfile = args['--atnt']
	threads = args['--t']

	if not os.path.exists(savePath):
		os.system(f'mkdir -p {savePath}')

	start = time.time()

	blast_search(toxinNT=tNtInfile, toxinPro=tProInfile, antitoxinNT=atInfile, dataPath=bacDatabase, outPath=savePath, threads=threads)

	end = time.time()
	t = end-start
	print('time to blast ', t)
