#!/usr/bin/env python3
import csv
import time
from docopt import docopt
import os
import numpy as np
import itertools as it
import subprocess

class Extend:
	def preprocess_extend(inFile, outDir, databasePath, extraLength):

		linetasID = int((subprocess.check_output(f'wc -l {inFile}', shell=True)).decode("utf-8").split(' ')[0]) - 1
		with open(inFile) as sumFile:

			seq = ''
			header = ''
			for index, line in enumerate(sumFile):
				if line.startswith('>'):
					if seq:
						Extend.extend_seq(header, seq, inFile, outDir, databasePath, extraLength)
				
					header = line.replace('\n', '')
					seq = ''
					continue

				elif index == linetasID:
					seq += line.replace('\n', '')
					Extend.extend_seq(header, seq, inFile, outDir, databasePath, extraLength)

			
				else:
					seq += line.replace('\n', '')


	def extend_seq(header, seq, inFile, outDir, databasePath, extraLength):

		startGaps = Extend.count_gaps(seq)
		endGaps = Extend.count_gaps(reversed(seq))

		tasID = header.split('_bac')[0].split('>')[1]
		bac = Extend.split_header(header, '_bac:', '_tasys:')
		TAsys = Extend.split_header(header, '_tasys:', '_genome:')
		genomeID = Extend.split_header(header, '_genome:', '_genomeType:')
		genomeType = Extend.split_header(header, '_genomeType:', '_blastType')
		blastType = Extend.split_header(header, '_blastType:', '_position:')
		startGenome, stopGenome = Extend.split_header(header, '_position:', '_genomeLength:').split('to')
		genomeLength = int(header.split('_genomeLength:')[1])

		dbfile = f'{databasePath}{Output.decide_file(blastType, genomeType)}'
		
		if int(startGenome) < int(stopGenome): # plus strand
			strand = '+'
			newStart, newStop = Extend.extend_coordinates(int(startGenome), int(stopGenome), startGaps, endGaps, genomeLength)
			Output.fastacmd(dbfile, genomeID, newStart, newStop, 1, outDir)
			alignLen = newStop - newStart
		
		elif int(startGenome) > int(stopGenome): # minus strand
			strand = '-'
			newStop, newStart = Extend.extend_coordinates(int(stopGenome), int(startGenome), endGaps, startGaps, genomeLength)
			alignLen = newStart - newStop
			Output.fastacmd(dbfile, genomeID, newStop, newStart, 2, outDir)
		else:
			print('ERROR:')
			exit()

		## Save new Summary
		# fileName = inFile.split('/')[-1].replace('.aln','_extend.csv')
		fileName = 'extend_summary.csv'
		lineArray = [tasID, bac, genomeType, strand, genomeID, TAsys, 'unknown', blastType, 0, newStart, newStop, alignLen, 0, 0, 0, genomeLength]

		Output.write_table(lineArray, f'{outDir}{fileName}')
		## Save Fasta
		saveFasta = f'{outDir}{TAsys}_{blastType}_extend.fa'
		newHeader = f'>{tasID}_bac:{bac}_tasys:{TAsys}_genome:{genomeID}_genomeType:{genomeType}_blastType:{blastType}_position:{newStart}to{newStop}_genomeLength:{genomeLength}\n'
		Output.write_fasta(saveFasta, newHeader)

	def split_header(string, key1, key2):
		value = string.split(key1)[1].split(key2)[0]
		return value

	def count_gaps(string):
		gaps = 0
		for i in string:
			if i == '-':
				gaps += 1
			else:
				break
		return gaps

	def extend_coordinates(a, b, minus, plus, genomeLength):
		aNew=1 if a<minus else (a-minus)
		bNew=genomeLength if (b+plus)>genomeLength else (b + plus)

		return aNew, bNew
		

class Output:
	def write_table(row, save):
		with open(save, mode='a') as tableFile:
			tableWriter = csv.writer(tableFile, delimiter='\t')
			tableWriter.writerow(row)

	def write_fasta(saveResults, newHeader):
		processID = os.getpid()
		save = os.path.dirname(saveResults)
		with open(f'{save}/{processID}_fastacmd_temp.txt') as cmdTemp:
			with open(saveResults, 'a') as saveFasta:
				for index, line in enumerate(cmdTemp.readlines()):
					if line.startswith('>'):
						saveFasta.write(newHeader)	
					else:
						saveFasta.write(line)

		os.remove(f'{save}/{processID}_fastacmd_temp.txt')

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


	def fastacmd(db, subjectID, start, stop, S, out):
		processID = os.getpid()
		os.system(f'fastacmd -d {db} -s {subjectID} -S {S} -L {start},{stop} -o {out}{processID}_fastacmd_temp.txt')

args = docopt("""extend.

Usage:
    box1_extend.py [-h] [--add INT] [--t INT] -i FILE -o DIR -n STR -d FILE

Options:
    -h --help   show this message and exit.
    -o DIR      Dir for output
    -i FILE     File with msa fasta
    -n STR      name of TA System
    -d FILE     genome database
    --add INT   additional length, must be divisible by 3 [default: 0]
    --t INT     number of reminders to give. [default: 1]

""")


if __name__ == '__main__':
	print(args)

	msa = args['-i']
	out = args['-o']
	bacDatabase = args['-d']
	extraLength = int(args['--add'])
	threads = int(args['--t'])

	start = time.time()

	if not os.path.exists(out):
		os.system(f'mkdir -p {out}')
	
	if (extraLength % 3) != 0:
		print('Error: parameter --add must be divisible by 3, to extend also the Proteinsequence.')
		exit()
	
		
	Extend.preprocess_extend(msa, out, bacDatabase, extraLength)

	end = time.time()
	print(f'time to extend: {end-start}')
