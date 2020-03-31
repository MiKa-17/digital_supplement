import csv
import time
from docopt import docopt
import os
import re

class InfernalFilter():

	def filter(databasePath, fileArray, taName, outFolder, locArray):

		tasNr = 0

		for index, inputFile in enumerate(fileArray):
			genChrom = locArray[index]

			with open(inputFile) as openFile:

				for line in openFile.readlines():
					if line.startswith('Query'):
						string = line.split()[1]
						TAsys = re.sub(r'_round\w+\_0001', '', string)
						print(TAsys)
						

					# elif '!' in line or '?' in line:
					elif '!' in line:
						lineArray = line.replace('\n', '').split()

						bac = InfernalFilter.replaceMultiple(('_').join(lineArray[12:]), [' ', '#', ',', '-', '/', '.',':'], '_').replace('__','_').strip('\n')
						
						eValue, _, _, genomeID, genomeStart, genomeStop, strand = line.split()[2:9]

						alignLen = abs(int(genomeStart) - int(genomeStop))

						tasNr += 1
						ta = TAsys.split('_')[0]
						tasID = f'{tasNr}i{ta}'
						### [tasID, bac, genChrom, strand, genomeID, TAsys, TAheader, blastType, eValue, genomeStart, genomeStop, alignLen, startTA, stopTA, lengthTA, genomeLength]
						newLine = [tasID, bac, genChrom, strand, genomeID, TAsys, 'not defined', 'infernal', eValue, genomeStart, genomeStop, alignLen, 0, 0, 0, 0]
						OutputCSV.write_table(f'{outFolder}infernal_summary.csv', newLine)

						### fastacmd
						OutputFasta.fasta([tasID, bac, 'infernal', genChrom, genomeID, TAsys, strand], outFolder, databasePath, genomeStart, genomeStop, taName)

					elif 'Hit alignments:' in line:
						break

					else:
						continue
	
	def replaceMultiple(mainString, toBeReplaces, newString):
		for elem in toBeReplaces :
			mainString = mainString.replace(elem, newString) if elem in mainString else mainString

		return mainString

class OutputCSV:

	def write_table(out, row):
		with open(out, mode='a') as tableFile:
			tableWriter = csv.writer(tableFile, delimiter='\t')
			tableWriter.writerow(row)


class OutputFasta():

	def fasta(group, outFolder, databasePath, start, stop, taName):
		tasID, bac, blastType, genomeType, genomeID, taType, strand = group

		dbfile = f'{databasePath}{OutputFasta.decide_file(blastType, genomeType)}'

		OutputFasta.fastacmd(dbfile, genomeID, start, stop, 1, outFolder) if strand == '+' else OutputFasta.fastacmd(dbfile, genomeID, stop, start, 2, outFolder) if strand == '-' else print('Error: incorrect strand label')

		# newHeader = f'>0_{bac}_{taType}_in_{genomeType}_position_{start}to{stop}\n'
		newHeader = f'>{tasID}_bac:{bac}_tasys:{taType}_genome:{genomeID}_genomeType:{genomeType}_blastType:{blastType}_position:{start}to{stop}_genomeLength:0\n'
		saveFasta = f'{outFolder}/{taType}_infernal.fa'
		OutputFasta.write_fasta(saveFasta, newHeader, outFolder)

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
		os.system(f'fastacmd -d {db} -s {subjectID} -S {S} -L {start},{stop} -o {save}{os.getpid()}_fastacmd_temp.txt')

	def write_fasta(saveResults, newHeader, save):
		with open(f'{save}{os.getpid()}_fastacmd_temp.txt') as cmdTemp:
			with open(saveResults, 'a') as saveFasta:
				[(saveFasta.write(newHeader) if line.startswith('>') else saveFasta.write(line)) for index, line in enumerate(cmdTemp.readlines())]

		os.remove(f'{save}{os.getpid()}_fastacmd_temp.txt')

args = docopt("""infernal.

Usage:
    box2_infernal_filter.py [-h] [--t INT] -l STR -i STR -d DIR -o DIR -n STR

Options:
    -h --help   show this message and exit.
    -o DIR      Folder for output
    -i STR      Infiles - multiple files possible (separate with comma)
    -d DIR      database
    -n STR      name of TA System
    -l STR    location of the input files (genome or chromosome) - multiple files possible (separate with comma)
    --t INT     number of reminders to give. [default: 1]

""")

if __name__ == '__main__':
	print(args)


	fInput = args['-i']
	fOutput = args['-o']
	dInput = args['-d']
	taName = args['-n']
	loc = args['-l']
	threads = int(args['--t'])

	if not os.path.exists(fOutput):
		os.system(f'mkdir -p {fOutput}')

	fileArray = fInput.split(',')
	locArray = loc.split(',')

	InfernalFilter.filter(dInput, fileArray, taName, fOutput, locArray)