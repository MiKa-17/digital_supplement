# /home/co68mol/programs/Python3.6/bin/python3

import time
from docopt import docopt
from collections import defaultdict
from Bio import SeqIO
import os
import sys

#projects/master/ta_system/TASYSTEM/blast/blast_filter

def seq_similarity(fastaFile):
	simSeq = False
	simDict = defaultdict(list) # Use the sequence as the key and then have a list of id's as the value
	for seq in SeqIO.parse(fastaFile, 'fasta'):
		simDict[str(seq.seq)].append(seq.id)
		if len(simDict[str(seq.seq)]) > 1:
			simSeq = True

	if simSeq: # if there are similar sequences: renew Fastafile
		base = fastaFile.split('.')[0]
		newfile = f'{base}_unique.fa'
		with open(newfile, 'w') as simOut:
			for seq, ids in simDict.items():
				simOut.write(">{}\n".format(''.join(ids))) #join header ids with same seq
				simOut.write(seq + "\n")

		return newfile
	else:
		return fastaFile #if there are no similar seq: take old fasta file

##blastn: 2.7.1+ Package: blast 2.7.1, build Oct 18 2017 19:57:24
##blastp: 2.7.1+ Package: blast 2.7.1, build Oct 18 2017 19:57:24


def blast_search(taName, toxinNT, toxinPro, antitoxinNT, dataPath, outPath, threads):

	for index, file in enumerate([toxinNT, toxinPro, antitoxinNT]):
		if not os.path.isfile(file):
			# if index == 2: # if no at
			# 	sys.exit()
			continue
		else:
			######################### generate file_name out of INPUT !!!!!!!!!!!!!!!!!!!!!
			fileName = file.split('/')[-1].split('.')[0]
			if index == 0:
				fileName=f't_nt_{taName}'
				blast_command('blastn', file, f'{dataPath}chromosome/full_chromosome_bacteria.fna', f'{outPath}{fileName}_chromosome_blastn.txt', threads)
				blast_command('blastn', file, f'{dataPath}genome/full_genome_bacteria.fna', f'{outPath}{fileName}_genome_blastn.txt', threads)
			elif index == 2:
				fileName=f'at_nt_{taName}'
				blast_command('blastn', file, f'{dataPath}chromosome/full_chromosome_bacteria.fna', f'{outPath}{fileName}_chromosome_blastn.txt', threads)
				blast_command('blastn', file, f'{dataPath}genome/full_genome_bacteria.fna', f'{outPath}{fileName}_genome_blastn.txt', threads)
			elif index == 1:
				fileName=f't_pro_{taName}'
				blast_command('blastp', file, f'{dataPath}proteome/full_proteome_bacteria_chr_un.faa', f'{outPath}{fileName}_chromosome_blastp.txt', threads)
				blast_command('blastp', file, f'{dataPath}proteome/full_proteome_bacteria_un.faa', f'{outPath}{fileName}_genome_blastp.txt', threads)
				blast_command('tblastn', file, f'{dataPath}chromosome/full_chromosome_bacteria.fna', f'{outPath}{fileName}_chromosome_tblastn.txt', threads)
				blast_command('tblastn', file, f'{dataPath}genome/full_genome_bacteria.fna', f'{outPath}{fileName}_genome_tblastn.txt', threads)

def blast_command(blastType, tasInput, databaseFile, out, threads):
	if os.path.isfile(databaseFile):
		os.system(f'nice {blastType} -num_threads {threads} -query {tasInput} -db {databaseFile} -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen gaps stitle sstrand" -out {out}') #add qlen to outputformat, to get query length

#docopt for command line arguments 	

args = docopt("""blast_search.

Usage:
    box_1_blast_search.py [-h] [--tnt FILE] [--tpro FILE] [--atnt FILE] [--t INT] -o DIR -d FILE -n STR

Options:
    -h --help    show this message and exit.
    -o DIR       Dir for output
    -d FILE      genome database
    -n STR       Name of the TAS
    --tnt FILE   FASTA file with one or more toxin sequences (nucleotide)
    --tpro FILE  FASTA file with one or more toxin sequences (protein)
    --atnt FILE  FASTA file with one ot more antitoxin sequences (nucleotide)
    --t INT      number of reminders to give. [default: 1]

""")

if __name__ == '__main__':
	print(args)
	bacDatabase = args["-d"]
	savePath = args["-o"]
	taName = args["-n"]


	tNtInfile = args['--tnt']
	tProInfile = args['--tpro']
	atInfile = args['--atnt']
	threads = args['--t']

	noSeq=[]

	if os.path.isfile(tNtInfile):
		noSeq.append(False)
		tNtInfile = seq_similarity(tNtInfile)#check Fasta-Infiles, if there are similar sequences with different header
	else:
		noSeq.append(True)
		print('no nucleotide toxin given.')
		
	if os.path.isfile(tProInfile):	
		noSeq.append(False)
		tProInfile = seq_similarity(tProInfile)#check Fasta-Infiles, if there are similar sequences with different header
	else:
		noSeq.append(True)
		print('no protein toxin given.')

	if os.path.isfile(atInfile):
		noSeq.append(False)	
		atInfile = seq_similarity(atInfile)#check Fasta-Infiles, if there are similar sequences with different header
	else:
		noSeq.append(True)
		print('no nucleotide antitoxin given.') 

	if all(noSeq):
		print('ERROR: please given at least one Toxin-Antitoxin-Sequence')
		exit()


	if not os.path.exists(savePath):
		os.system(f'mkdir -p {savePath}')

	start = time.time()

	blast_search(taName=taName, toxinNT=tNtInfile, toxinPro=tProInfile, antitoxinNT=atInfile, dataPath=bacDatabase, outPath=savePath, threads=threads)

	end = time.time()
	t = end-start
	print('time to blast ', t)
