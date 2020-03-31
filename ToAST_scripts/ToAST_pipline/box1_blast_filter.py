import time
from docopt import docopt
import numpy as np
# from multiprocessing import Pool
# from functools import partial
# import math
import os
import csv

import pandas as pd

#Ordnerstruk:
#projects/master/ta_system/TASYSTEM/blast/blast_filter

def filter_eValue(databasePath, InputDir, eValue, pident, writeSeq, taName):

	tasNr=0

	for file in os.listdir(InputDir):
		file = os.path.join(InputDir, file)
		if os.stat(file).st_size != 0:
			with open(file) as openFile:
				for line in openFile.readlines():
					lineArray = line.split('\t')
					currentEValue = lineArray[10]
					currentPident = lineArray[2]


					if eValue >= float(currentEValue) and float(currentPident) >= pident:
						### Information from the file name
						fileArray = file.split('/')[-1].split('_')
						blastType = fileArray[-1].split('.')[0]
						genomeType = fileArray[-2]
						taType = '_'.join(fileArray[:2]).split('.')[0] + '_' + taName

						### Information from the blast search
						taID = lineArray[0]
						subjectID = lineArray[1].split('|')[1]
						alignmentLength = lineArray[3]
						startTA, stopTA = int(lineArray[6]) , int(lineArray[7])
						startGenome, stopGenome = int(lineArray[8]), int(lineArray[9])
						lengthTA = int(lineArray[12])
						lengthSubject = int(lineArray[13])
						gaps = int(lineArray[14]) # remove line berak
						bac = replaceMultiple(('_').join(lineArray[15:-1]), [' ', '#', ',', '-', '/', '.',':'], '_').replace('__','_').strip('\n')

						tasNr += 1
						ta = taType.split('_')[0]
						tasID = f'{tasNr}b{ta}'

						dbfile = f'{databasePath}{decide_file(blastType, genomeType)}'
						
						e = '{0:.7f}'.format(eValue)
						saveFasta = f'{save}{taType}_{blastType}_{e}'

						if not os.path.exists(save):
							os.system(f'mkdir -p {save}')

						if startGenome < stopGenome:
							if summaryBool:
								a = [tasID, bac, genomeType, '+', subjectID, taType, taID, blastType, str(currentEValue), str(startGenome), str(stopGenome), str(alignmentLength), str(startTA), str(stopTA), str(lengthTA), str(lengthSubject)]
								write_table(a, f'blast_results.csv')
							if writeSeq:
								fastacmd(dbfile, subjectID, startGenome, stopGenome, 1)

						elif startGenome > stopGenome:
							if summaryBool:
								a = [tasID, bac, genomeType, '-', subjectID, taType, taID, blastType, str(currentEValue), str(startGenome), str(stopGenome), str(alignmentLength), str(startTA), str(stopTA), str(lengthTA), str(lengthSubject)]
								write_table(a, f'blast_results.csv')
							if writeSeq:
								fastacmd(dbfile, subjectID, stopGenome, startGenome, 2)

						if writeSeq:
							newHeader = f'>{tasID}_bac:{bac}_tasys:{taType}_genome:{subjectID}_genomeType:{genomeType}_blastType:{blastType}_position:{startGenome}to{stopGenome}_genomeLength:{lengthSubject}\n'
							write_fasta(f'{saveFasta}.fa', newHeader)

def replaceMultiple(mainString, toBeReplaces, newString):
	for elem in toBeReplaces :
		mainString = mainString.replace(elem, newString) if elem in mainString else mainString

	return mainString
							
def write_table(row, fileName):
	# if type(row) == list:
	# 	print(type(row))

	with open(save + fileName, mode='a') as tableFile:

		# np.savetxt(tableFile, [row], delimiter="\t", fmt='%s %s %s %s %s %s %s %s %s %i %i')
		tableWriter = csv.writer(tableFile, delimiter='\t')
		tableWriter.writerow(row)

def write_fasta(saveResults, newHeader):
	processID = os.getpid()
	with open(f'{save}{processID}_fastacmd_temp.txt') as cmdTemp:
		with open(saveResults, 'a') as saveFasta:
			for index, line in enumerate(cmdTemp.readlines()):
				if line.startswith('>'):
					saveFasta.write(newHeader)
				else:
					saveFasta.write(line)

	os.remove(f'{save}{processID}_fastacmd_temp.txt')

def decide_file(blastType, genomeType):
	if blastType == "blastp":
		if genomeType == "chromosome":
			return('/proteome/full_proteome_bacteria_chr_un.faa')
		elif genomeType == "genome":
			return('/proteome/full_proteome_bacteria_un.faa')
	elif blastType == "blastn" or blastType == "tblastn":
		if genomeType == "chromosome":
			return('/chromosome/full_chromosome_bacteria.fna')
		elif genomeType == "genome":
			return('/genome/full_genome_bacteria.fna')

def fastacmd(db, subjectID, start, stop, S):
	processID = os.getpid()
	os.system(f'fastacmd -d {db} -s {subjectID} -S {S} -L {start},{stop} -o {save}{processID}_fastacmd_temp.txt')

def chunk(l, n):
	for i in range(0, len(l), n):
		yield l[i:i + n]

args = docopt("""blast_filter.
Usage:
    box_1_blast_filter.py [-h] [--evalue INT] [--pident INT ] [--summary] [--seq] [--t INT] -d DIR -i DIR -o DIR -n STR

Options:
    -h --help        show this message and exit.
    -o DIR           Dir for output
    -i DIR           Dir with blast results
    -d DIR           database
    -n STR           Name of the TA-System
    --evalue INT     evalue to filter [default: 0.0000001]
    --pident INT     percentage of identical matches [default: 00.00]
    --seq            save fasta with sequences
    --summary        write summarry
    --t INT          number of reminders to give. [default: 1]

""")


if __name__ == '__main__':

	print(args)
	
	blastDir = args['-i']
	save = args['-o']
	databasePath = args['-d']
	taName = args['-n']
	eValue = float(args['--evalue'])
	pident = float(args['--pident'])

	summaryBool = args['--summary']
	seqBool = args['--seq']
	threads = args['--t']


	#print("start with:")
	#print(blastDir)

	start = time.time()

	filter_eValue(databasePath=databasePath, InputDir=blastDir, eValue=eValue, pident=pident, writeSeq=seqBool, taName=taName)

	# chunksList = list(chunk(os.listdir(blastDir), math.ceil(len(os.listdir(blastDir))/int(threads))))

	# pool = Pool(processes = int(threads))
	# partialBF = partial(filter_eValue, InputDir=blastDir, eValue=eValue, writeSeq=seqBool)
	# pool.map(partialBF, chunksList)

	end = time.time()
	print(f'time to filter')
	print(end-start)
