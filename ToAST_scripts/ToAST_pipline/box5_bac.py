#!/usr/bin/env python3
import csv
import time
from docopt import docopt
import os
import numpy as np
import itertools as it
import subprocess
import pandas as pd
import re




class TasInBacteria:

	def read_comb_summary(inFolder, outFolder, dbFile):

		# for dirs in os.listdir(inFolder):
			# folder = os.listdir(f'{inFolder}{dirs}/combine_bi/')
		# combdir = f'{inFolder}{dirs}/3_combined_analysis/'
		combdir = inFolder

		if os.path.isdir(combdir):

			ctasDir=f'{combdir}cTAS/'

			# dir=$(ls -d "$combDir"cTAS/*/ | tail -1)

			folder = [ name for name in os.listdir(ctasDir) if os.path.isdir(os.path.join(ctasDir, name)) ]
			
			maxFolderNumb = max([int(f.split('v')[-1]) for f in folder])
			maxFolder = f'cTASv{maxFolderNumb}'

			summary=f'{combdir}cTAS/{maxFolder}/summary_cTAS.csv'

			data = pd.read_csv(summary, header=None, names=['tasID', 'bac', 'genome_chrom', 'strand', 'genomeID', 'TAsys', 'TAheader', 'blastType', 'eValue', 'genomeStart', 'genomeStop', 'alignLen', 'startTA', 'stopTA', 'lengthTA', 'genomeLength'], sep='\t', dtype={'eValue': object, 'genomeStart': int, 'genomeStop': int, 'startTA': int, 'stopTA': int, 'lengthTA': int, 'genomeLength': int})
			data = data.sort_values(by=['bac'])

			data['bac'] = data.bac.str.replace('UNVERIFIED_ORG:_', '')
			data['bac'] = data.bac.str.replace('UNVERIFIED_ORG_', '')

			new = data['bac'].str.split('_', n=2, expand=True)
			data.insert(0, 'bac_species', new[0] + '_' + new[1])
			data.insert(1, 'bac_strain', new[2])
			# data.drop(columns=['bac'], inplace=True)

			groupData = data.groupby(['bac_species'])

			for key, group in groupData:
				bac = group['bac_species'].iloc[0]
				if not os.path.exists(f'{outFolder}/{bac}/'):
					os.system(f'mkdir -p {outFolder}/{bac}/')
				group.is_copy = None
				group.drop(columns=['bac_species', 'bac_strain'], inplace=True) # inplce=overwrite group, no coy
				group.to_csv(f'{outFolder}/{bac}/{bac}_summary.csv', sep='\t', header=False, index=False, mode='a')

				for index, line in group.iterrows():
					OutputFastaTas.fasta(line, f'{outFolder}/{bac}/', dbFile, bac)

	# def write_summery():
		# bac _ taSys1 _ numer of cTAs on one genome _ taSys2 _ number of cTAs on one genome ...

class OutputCsv():

	def write_csv(out, data):
		with open(out, 'ab') as f:
			np.savetxt(f, data, delimiter="\t", fmt='%s')


class OutputFastaTas():

	def fasta(group, outFolder, databasePath, bacSpecies):

		tasID = group['tasID']
		bac = group['bac']
		blastType = group['blastType']
		genomeType = group['genome_chrom']
		genomeID = group['genomeID']
		taSys = group['TAsys']
		taType = '_'.join(taSys.split('_')[:2])

		strand = group['strand']
		genomeLength = group['genomeLength']
		start = group['genomeStart']
		stop = group['genomeStop']

		dbfile = f'{databasePath}{OutputFasta.decide_file(blastType, genomeType)}'

		OutputFasta.fastacmd(dbfile, genomeID, start, stop, 1, outFolder) if strand == '+' else OutputFasta.fastacmd(dbfile, genomeID, stop, start, 2, outFolder) if strand == '-' else print('Error: incorrect strand label') & exit()

		newHeader = f'>{tasID}_bac:{bac}_tasys:{taSys}_genome:{genomeID}_genomeType:{genomeType}_blastType:{blastType}_position:{start}to{stop}_length:{genomeLength}\n'
		saveFasta = f'{outFolder}/{taType}_{bacSpecies}.fa'
		OutputFasta.write_fasta(saveFasta, newHeader, outFolder)

class OutputFasta():

	def decide_file(blastType, genomeType):

		if blastType == "blastp":
			dbPath = 'proteome/full_proteome_bacteria_chr_un.faa' if genomeType=='chromosome' else 'proteome/full_proteome_bacteria_un.faa' if genomeType=='genome' else print('Error: incorrect genomeType label')  & exit()
			return dbPath
		elif blastType == "blastn" or blastType == "tblastn" or blastType == 'infernal':
			dbPath = 'chromosome/full_chromosome_bacteria.fna' if genomeType=='chromosome' else 'genome/full_genome_bacteria.fna' if genomeType=='genome' else print('Error: incorrect genomeType label')  & exit()
			return dbPath
		else:
			print('Error: incorrect blastType label')
			exit()


	def fastacmd(db, subjectID, start, stop, S, save):
		processID = os.getpid()
		os.system(f'fastacmd -d {db} -s {subjectID} -S {S} -L {start},{stop} -o {save}{processID}_fastacmd_temp.txt')

	def write_fasta(saveResults, newHeader, save):
		processID = os.getpid()
		with open(f'{save}{processID}_fastacmd_temp.txt') as cmdTemp:
			with open(saveResults, 'a') as saveFasta:
				for index, line in enumerate(cmdTemp.readlines()):
					if line.startswith('>'):
						saveFasta.write(newHeader)
					else:
						saveFasta.write(line)
		os.remove(f'{save}{processID}_fastacmd_temp.txt')


### inputfolder: projects/master/ta_systems
### destination: ../ta_systems/*/combine_bi/blast_infernal_summary.csv


args = docopt("""box5_bac.

Usage:
    box5_bac.py [-h] [--t INT] -i DIR -o DIR -d DIR

Options:
    -h --help   show this message and exit.
    -o DIR      Folder for output
    -i DIR     File with msa fasta
    -d DIR     genome database
    --t INT     number of reminders to give. [default: 1]

""")


if __name__ == '__main__':
	print(args)

	inFolder = args['-i']
	outFolder = args['-o']
	bacDatabase = args['-d']
	threads = int(args['--t'])

	start = time.time()

	if not os.path.exists(outFolder):
		os.system(f'mkdir -p {outFolder}')
	
		
	TasInBacteria.read_comb_summary(inFolder, outFolder, bacDatabase)

	end = time.time()
	print(f'time to: {end-start}')
