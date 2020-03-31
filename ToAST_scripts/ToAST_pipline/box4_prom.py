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
from Bio import SeqIO
# from fuzzywuzzy import fuzz
import re

#  /home/co78wid/projects/master/ta_systems/aapa_isoa/combine_bi/cTASv0/full_summary_merged_complete_tas.csv


class PromoterAnalysis:

	def is_fasta(inFile):


	    with open(inFile, "r") as handle:
	        fasta = SeqIO.parse(handle, "fasta")
	        return any(fasta)

	def get_prom_seq(inFile, outDir, dbDir):

		data = pd.read_csv(inFile, header=None, names=['tasID', 'bac', 'genome_chrom', 'strand', 'genomeID', 'TAsys', 'TAheader', 'blastType', 'eValue', 'genomeStart', 'genomeStop', 'alignLen', 'startTA', 'stopTA', 'lengthTA', 'genomeLength'], sep='\t', dtype={'eValue': object, 'genomeStart': int, 'genomeStop': int, 'startTA': int, 'stopTA': int, 'lengthTA': int, 'genomeLength': int})

		for index, line in data.iterrows():

			# strand = line['strand']
			startGenome = line['genomeStart']
			stopGenome = line['genomeStop']

			if line['blastType'] == 'infernal':
				if stopGenome > startGenome: # plus strand
					line['strand'] = '+'
					strand = '+'
					for i in range(50, 0, -1): # count 50 to 40
						try:
							newStart, newStop = PromoterAnalysis.PromoterAnalysis_coordinates(int(startGenome), int(stopGenome), i, 0, None)
							OutputFastaTas.fasta(line, f'{outDir}', dbDir, newStart, newStop, startGenome, stopGenome)
							break
						except Exception as e:
							print('exept')
							raise e
				elif stopGenome < startGenome: # minus strand
					line['strand'] = '-'
					strand = '-'
					for i in range(50, 0, -1):
						try:
							newStop, newStart = PromoterAnalysis.PromoterAnalysis_coordinates(int(stopGenome), int(startGenome), 0, i, None)
							OutputFastaTas.fasta(line, f'{outDir}', dbDir, newStart, newStop, startGenome, stopGenome)
							break
						except Exception as e:
							print('exept')
							raise e
				else:
					print('ERROR:')
					exit()

			else:
				genomeLength = line['genomeLength']

				if stopGenome > startGenome: # plus strand
					line['strand'] = '+'
					strand = '+'
					newStart, newStop = PromoterAnalysis.PromoterAnalysis_coordinates(int(startGenome), int(stopGenome), 50, 0, genomeLength)
					OutputFastaTas.fasta(line, f'{outDir}', dbDir, newStart, newStop, startGenome, stopGenome)
				elif stopGenome < startGenome: # minus strand
					line['strand'] = '-'
					strand = '-'
					newStop, newStart = PromoterAnalysis.PromoterAnalysis_coordinates(int(stopGenome), int(startGenome), 10, 0, genomeLength)
					OutputFastaTas.fasta(line, f'{outDir}', dbDir, newStart, newStop, startGenome, stopGenome)
				else:
					print('ERROR:')
					exit()



	def PromoterAnalysis_coordinates(a, b, minus, plus, genomeLength):

		aNew=1 if a<minus else (a-minus)
		if genomeLength:
			bNew=genomeLength if (b+plus)>genomeLength else (b + plus)
		else:
			bNew=(b + plus)

		return aNew, bNew


	def find_prom(inFile, outDir):

		with open(inFile, "r") as handle:
			for record in SeqIO.parse(handle, "fasta"):
				# startStop = PromoterAnalysis.split_header(record.id, '_position:', '_length').split('to')

				header = record.id
				sequence = record.seq
				tasID = header.split('_bac:')[0]
				bac = PromoterAnalysis.split_header(header, '_bac:', '_tasys:')
				TAsys = PromoterAnalysis.split_header(header, '_tasys:', '_genome:')
				genomeID = PromoterAnalysis.split_header(header, '_genome:', '_genomeType:')
				genomeType = PromoterAnalysis.split_header(header, '_genomeType:', '_blastType')
				blastType = PromoterAnalysis.split_header(header, '_blastType:', '_position:')
				start, stop = PromoterAnalysis.split_header(header, '_position:', '_genomeLength:').split('to')
				alignLen = abs(int(start) - int(stop))
				genomeLength = (header.split('_genomeLength:')[1])

				# start = startStop[0]
				# stop = startStop[1]

				if stop > start:
					strand = '+'
				else:
					strand = '-'

				fileName = inFile.split('/')[-1].replace('.aln','_extend.csv')
				newRow = [(tasID, bac, genomeType, strand, genomeID, TAsys, 'unknown', blastType, 0, start, stop, alignLen, 0, 0, 0, genomeLength)]
				#AGGAGG or GAGG for SD
				# -35 15-25nt TATA 15-50nt RBS(SD!) 6nt ATG minlenofseq start-stopTA GGAGG
				p1 = 'TTGACA.{15,25}TATAAT.{4,10}GGAGG.{4,10}ATG'
				p2 = 'TTGACA.{15,25}TATAAT.{4,10}GGAGG.{4,10}GTG'
				p3 = 'TATAAT.{4,10}GGAGG.{4,10}ATG'
				p4 = 'TATAAT.{4,10}GGAGG.{4,10}GTG'
				p5 = 'GGAGG.{4,10}ATG'
				p6 = 'GGAGG.{4,10}GTG'
				p7 = '.{40,60}ATG'
				p8 = '.{40,60}GTG'

				p9 = 'TTGACA.{15,60}TATAAT.{4,25}ATG'
				p10 = 'TTGACA.{15,60}TATAAT.{4,25}GAG'

				p11 = 'TATAAT.{4,25}ATG'
				p12 = 'TATAAT.{4,25}GAG'

				p13 = 'TTGACA.{34,60}ATG'
				p14 = 'TTGACA.{34,60}GAG'

				p15 = 'TTGACA.{20,40}GGAGG.{4,10}ATG'
				p16 = 'TTGACA.{20,40}GGAGG.{4,10}GTG'


				#### match 4 promoter elements
				match4 = PromoterAnalysis.search_prom(p1, p2, record, newRow, 'prom_complete_header.txt')

				if not match4:

					#### match 3 promoter elements
					match31 = PromoterAnalysis.search_prom(p3, p4, record, newRow, 'prom_tata_sd_start_header.txt')
					match32 = PromoterAnalysis.search_prom(p9, p10, record, newRow, 'prom_35_tata_start_header.txt')
					match33 = PromoterAnalysis.search_prom(p15, p16, record, newRow, 'prom_35_sd_start_header.txt')

					if not any([match31, match32, match33]):

						#### match 2 promoter elements
						match21 = PromoterAnalysis.search_prom(p5, p6, record, newRow, 'prom_sd_start_header.txt')
						match22 = PromoterAnalysis.search_prom(p11, p12, record, newRow, 'prom_tata_start_header.txt')
						match23 = PromoterAnalysis.search_prom(p13, p14, record, newRow, 'prom_35_start_header.txt')

						if not any([match21, match22, match23]):
							##### match 1 promoter element
							match11 = PromoterAnalysis.search_prom(p7, p8, record, newRow, 'prom_start_header.txt')




	def split_header(string, key1, key2):
		value = string.split(key1)[1].split(key2)[0]
		return value

	def search_prom(prom1, prom2, record, newRow, file):
		result1 = []
		result2 = []
		result1 += re.finditer(prom1, str(record.seq)) or ['Error']
		result2 += re.finditer(prom2, str(record.seq)) or ['Error']
		if result1 or result2:
			OutputProm.write_prom(record.id, f'{outDir}{file}')
			# write prom to csv
			OutputCsv.write_csv(f'{outDir}promoter_summary.csv', newRow)
			return True
		return False


class OutputProm():


	def write_prom(header, outFile):

		with open(outFile, 'a') as hF:
			header=header.lstrip('>')
			hF.write(f'{header}\n')

	def write_fasta(inFile, outDir, fasta=False):
		if fasta:
			# get sequence dirct form file - promoter remains
			print('is fasta')

		else:
			if inFile == f'{outDir}promInput/t_nt_prom.fa':
				OutputProm.get_fasta_without_prom(inFile, outDir, 't_nt')
			elif inFile == f'{outDir}promInput/t_pro_prom.fa':
				OutputProm.get_fasta_without_prom(inFile, outDir, 't_pro')
			elif inFile == f'{outDir}promInput/at_nt_prom.fa':
				OutputProm.get_fasta_without_prom(inFile, outDir, 'at_nt')
			else:
				print('ERROR: no valid input-file ending')
				exit()

	def get_fasta_without_prom(inFile, outDir, tas):


		tasDir = os.path.dirname(os.path.dirname(outDir)) + '/'

		combDir=subprocess.check_output(f'ls -d {tasDir}3_combined_analysis/cTAS/*/ | tail -1', shell=True).decode("utf-8").rstrip('\n') #get dir with cTAS with biggest distance (most cTAS)

		file = subprocess.check_output(f'ls {combDir} | grep "^{tas}*"', shell=True).decode("utf-8").rstrip('\n')

		for hF in ['prom_complete_header.txt', 'prom_tata_sd_start_header.txt', 'prom_35_tata_start_header.txt', 'prom_35_sd_start_header.txt', 'prom_sd_start_header.txt', 'prom_tata_start_header.txt', 'prom_35_start_header.txt', 'prom_start_header.txt']:
			if os.path.isfile(f'{outDir}{hF}'):
				outFile = hF.replace('_header.txt', '.fa')
				os.system(f'seqtk subseq {combDir}{file} {outDir}{hF} >> {outDir}{tas}_ib_{outFile}')


class OutputCsv():

	def write_csv(out, data):
		with open(out, 'ab') as newf:
			np.savetxt(newf, data, delimiter="\t", fmt='%s')

class OutputFastaTas():

	def fasta(group, outDir, databasePath, start, stop, oldStart, oldStop):
		tasID = group['tasID']

		bac = group['bac']
		blastType = group['blastType']
		genomeType = group['genome_chrom']
		genomeID = group['genomeID']
		taSys = group['TAsys']
		taType = '_'.join(taSys.split('_')[:2])

		strand = group['strand']
		genomeLength = group['genomeLength']

		dbfile = f'{databasePath}{OutputFasta.decide_file(blastType, genomeType)}'

		OutputFasta.fastacmd(dbfile, genomeID, start, stop, 1, outDir) if strand == '+' else OutputFasta.fastacmd(dbfile, genomeID, stop, start, 2, outDir) if strand == '-' else print('Error: incorrect strand label') & exit()

		newHeader = f'>{tasID}_bac:{bac}_tasys:{taSys}_genome:{genomeID}_genomeType:{genomeType}_blastType:{blastType}_position:{oldStart}to{oldStop}_genomeLength:{genomeLength}\n'
		saveFasta = f'{outDir}{taType}_prom.fa'
		OutputFasta.write_fasta(saveFasta, newHeader, outDir)

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
		#try:
		#	os.system(f'fastacmd -d {db} -s {subjectID} -S {S} -L {start},{stop} -o {save}{processID}_fastacmd_temp.txt')
		#except Exception as e:
		#	raise
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


args = docopt("""prom.

Usage:
    box4_prom.py [-h] [--t INT] -i FILE -o DIR -d FILE

Options:
    -h --help   show this message and exit.
    -o DIR      Dir for output
    -i FILE     File with msa fasta
    -d FILE     genome database
    --t INT     number of reminders to give. [default: 1]

""")


if __name__ == '__main__':
	print(args)

	inFile = args['-i']
	outDir = args['-o']
	bacDatabase = args['-d']
	threads = int(args['--t'])

	start = time.time()

	if os.path.isfile(inFile):
		if not os.path.exists(outDir):
			os.system(f'mkdir -p {outDir}')

		if PromoterAnalysis.is_fasta(inFile): #check if file is a fasta
			PromoterAnalysis.find_prom(inFile, outDir)
			OutputProm.write_fasta(inFile, outDir, fasta=True)
		else:
			if not os.path.exists(f'{outDir}promInput/'):
				os.system(f'mkdir -p {outDir}promInput/')

			PromoterAnalysis.get_prom_seq(inFile, f'{outDir}promInput/', bacDatabase)
			if os.path.isfile(f'{outDir}promInput/t_nt_prom.fa'):
				PromoterAnalysis.find_prom(f'{outDir}promInput/t_nt_prom.fa', outDir)
				OutputProm.write_fasta(f'{outDir}promInput/t_nt_prom.fa', outDir)

			if os.path.isfile(f'{outDir}promInput/t_pro_prom.fa'):
				PromoterAnalysis.find_prom(f'{outDir}promInput/t_pro_prom.fa', outDir)
				OutputProm.write_fasta(f'{outDir}promInput/t_pro_prom.fa', outDir)

			if os.path.isfile(f'{outDir}promInput/at_nt_prom.fa'):
				PromoterAnalysis.find_prom(f'{outDir}promInput/at_nt_prom.fa', outDir)
				OutputProm.write_fasta(f'{outDir}promInput/at_nt_prom.fa', outDir)


		end = time.time()
		summaryFile=f'{outDir}promoter_summary.csv'
		os.system(f'sort -u {summaryFile} > {outDir}summaryFile.temp') # make summary unique
		os.system(f'mv {outDir}summaryFile.temp {summaryFile}')
		os.system(f'find {outDir} -size 0 -print0 |xargs -0 rm -f --') #remove empty files
		print(f'time to : {end-start}')
	else:
		print(f'ERROR: Input file dose not exists: {inFile}')
		exit()
