#!/usr/bin/env python3
import sys
import csv
from itertools import groupby
import time
from docopt import docopt
import os
import pandas as pd
import numpy as np
import itertools
from collections import defaultdict
import re
import multiprocessing as mp
from joblib import Parallel, delayed
from functools import partial
# import math


class CombineSummarys():

	# blast_filter => round => summary
	def merge_summarys(folder):

		current = folder.split('/')[-2]
		filterPath = folder.split('round')[0]
		if current == 'round1':
			return f'{folder}blast_results.csv'
		else:
			prev = 'round' + str(int(current.lstrip('round')) - 1)
			if os.path.isfile(f'{filterPath}{prev}/full_summary_merged.csv'):
				os.system(f'cat {filterPath}{prev}/full_summary_merged.csv {folder}/blast_results.csv >> {folder}full_summary.csv')
			else:
				os.system(f'for round in {filterPath}round*/; do for summary in "$round"blast_results.csv; do echo "$summary"; cat "$summary" >> {folder}full_summary.csv; done; done')

			return f'{folder}full_summary.csv'

class MergeSummary():

	def filter_summary(data, outFolder, databasePath, threads, taName):
		#####group summary
		data = np.sort(np.unique(data, axis=0), order=['bac'])
		a = np.unique(data['bac'], return_counts=True)[1]
		groupData = np.split(data, np.cumsum(a))
		# global groups
		# groups=[]


		for a in groupData:
			a = np.sort(a, order=['genome_chrom'])
			part = np.split(a, np.cumsum(np.unique(a['genome_chrom'], return_counts=True)[1]))
			for b in part:
				b = np.sort(b, order=['strand'])
				p = np.split(b, np.cumsum(np.unique(b['strand'], return_counts=True)[1]))
				for c in p:
					c = np.sort(c, order=['genomeID'])
					p = np.split(c, np.cumsum(np.unique(c['genomeID'], return_counts=True)[1]))
					for d in p:
						d = np.sort(d, order=['TAsys'])
						p =  np.split(d, np.cumsum(np.unique(d['TAsys'], return_counts=True)[1]))
						for f in p:
							f = np.sort(f, order=['blastType'])
							p =  np.split(f, np.cumsum(np.unique(f['blastType'], return_counts=True)[1]))
							if p[0].size != 0:
								# groups.append(p[:-1][0])

								MergeSummary.test_overlapp(p[:-1][0], dist=0.0, itera=None, outFolder=outFolder, databasePath=databasePath, taName=taName)
		# self.data = None
		# retLst = ray.get()
		# partialFunction = partial(MergeSummary.test_overlapp, dist=0.0, itera=None, outFolder=outFolder, databasePath=databasePath, taName=taName)
		# retLst = Parallel(n_jobs=threads)(delayed(partialFunction)([i, group]) for i, group in enumerate(groups))

	def test_overlapp(group, dist, itera, outFolder, databasePath, taName):

		# group = group[1]
		# index = group[0]

		# # np.delete(groups, index)
		# groups.remove(group)
		# # print(len(groups))
		# # exit()


		if itera is None:
			itera = 0

		if len(group) > 1: #group has more than one mamber
			strand = str(group['strand'][0])
			### sort group to start/stop according to strand
			group = np.sort(group, order=['genomeStart', 'genomeStop']) if strand == '+' else np.sort(group, order=['genomeStop', 'genomeStart']) if strand == '-' else print('Error: incorrect strand label')  & exit()
			for i in range(len(group)-1): #for each line in group
				dtype = [('tasID', 'U100'), ('bac', 'U100'), ('genome_chrom', 'U100') , ('strand', 'U100'), ('genomeID', 'U100'), ('TAsys', 'U100'), ('TAheader', 'U100'), ('blastType', 'U100'), ('eValue', 'U100'), ('genomeStart', int), ('genomeStop', int), ('alignLen', int), ('startTA', int), ('stopTA', int), ('lengthTA', int), ('genomeLength', int)]

				twoRows = np.asarray([group[i], group[i+1]], dtype=dtype) #take the two first lines of the group

				start1 = group['genomeStart'][i]
				stop1 = group['genomeStop'][i]
				start2 = group['genomeStart'][i+1]
				stop2 = group['genomeStop'][i+1]

				if strand == '-':
					pos1, pos2, pos3, pos4 = stop2, start2, stop1, start1 #swap positions according to the strand
					newStop, newStart, length = MergeSummary.min_max_len(twoRows, 'genomeStop', 'genomeStart') # claculate the coordinates of the minimal and the maximal positions as newStart and newStop of the hypothetical overlap/match
				elif strand == '+':
					pos1, pos2, pos3, pos4 = start1, stop1, start2, stop2 #swap positions according to the strand
					newStart, newStop, length = MergeSummary.min_max_len(twoRows, 'genomeStart', 'genomeStop') # claculate the coordinates of the minimal and the maximal positions as newStart and newStop of the hypothetical overlap/match
				else:
					print('Error: incorrect strand label')
					exit()

				if MergeSummary.overlaps(pos1, pos2, pos3, pos4, dist): #check overplap for the start/stop of the two rows
					newRow = MergeSummary.new_row(twoRows, i, newStart, newStop, length, dtype) #generate new row for overlap

					newRow = np.array(newRow, dtype=group.dtype)

					###if there are more than two rows in the group, combine newRow with the rest of the group
					newGroup = np.append(newRow, group[i+2:], axis=0) if (len(group)-i)>2 else newRow
					itera += 1
					return MergeSummary.test_overlapp(newGroup, 0.0, itera, outFolder, databasePath, taName)

				else: #if there is no overlap
					if i+1 == len(group)-1: # last two rows of the group without overlap - append both rows
						OutputCsv.write_csv(f'{outFolder}full_summary_merged.csv', twoRows)
						for j in range(2):
							OutputFasta.fasta(twoRows, outFolder, databasePath, twoRows['genomeStart'][j], twoRows['genomeStop'][j], taName)

					else: #if there are more than 2 rows in group, and the first two rows have no overlap
						prevRow = np.asarray([group[i]], dtype=dtype) # only write first row to file
						OutputCsv.write_csv(f'{outFolder}full_summary_merged.csv', prevRow)
						OutputFasta.fasta(prevRow, outFolder, databasePath, int(prevRow['genomeStart']), int(prevRow['genomeStop']), taName)

		else: #if there is only one mamber in group, write member to files
			OutputCsv.write_csv(f'{outFolder}full_summary_merged.csv', group)
			OutputFasta.fasta(group, outFolder, databasePath, group['genomeStart'][0], group['genomeStop'][0], taName)



	def new_row(data, index, start, stop, length, dtype):
		return[(data['tasID'][0], data['bac'][0], data['genome_chrom'][0], data['strand'][0], data['genomeID'][0], data['TAsys'][0], 'NaN - rows are merged', data['blastType'][0], 'NaN - rows are merged', start, stop, length, data['startTA'][0], data['stopTA'][0], data['lengthTA'][0], data['genomeLength'][0])]

	def min_max_len(group, columnA, columnB):
		minPos = np.amin(group[columnA])
		maxPos = np.amax(group[columnB])
		length = str((maxPos)-(minPos)+1)

		return [minPos, maxPos, length]

	def overlaps(a, b, c, d, dist):

		overlap = True
		### stopOne<startTwo or stopTwo<startOne for + strand
		### startTwo<stopOne or startOne<stopTwo for - strand
		if (b+dist)<c or (d+dist)<a:
			overlap = False
		return overlap

class OutputCsv():

	def write_csv(out, data):
		with open(out, 'ab') as f:
			np.savetxt(f, data, delimiter="\t", fmt='%s')

class OutputFasta():

	def fasta(group, outFolder, databasePath, start, stop, taName):
		tasID = group['tasID'][0]
		bac = group['bac'][0]
		blastType = group['blastType'][0]
		genomeType = group['genome_chrom'][0]
		genomeID = group['genomeID'][0]
		taType = group['TAsys'][0]
		strand = group['strand'][0]
		genomeLength = group['genomeLength'][0]

		dbfile = f'{databasePath}{OutputFasta.decide_file(blastType, genomeType)}'

		# OutputFasta.fastacmd(dbfile, genomeID, start, stop, 1, outFolder) if strand == '+' else OutputFasta.fastacmd(dbfile, genomeID, stop, start, 2, outFolder) if strand == '-' else print('Error: incorrect strand label') & exit()

		OutputFasta.fastacmd(dbfile, genomeID, start, stop, 1, outFolder) if (stop > start) else OutputFasta.fastacmd(dbfile, genomeID, stop, start, 2, outFolder) if(stop < start) else print('Error: incorrect strand label') & exit()

		newHeader = f'>{tasID}_bac:{bac}_tasys:{taType}_genome:{genomeID}_genomeType:{genomeType}_blastType:{blastType}_position:{start}to{stop}_genomeLength:{genomeLength}\n'
		saveFasta = f'{outFolder}/{taType}_{blastType}_merged.fa'
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



#################### FIND complete TAS

class StrToBool():

	def str_to_bool(s):
		if s == 'True':
			return True
		elif s == 'False':
			return False
		else:
			raise ValueError # evil ValueError that doesn't tell you what the wrong value was


class FindCompleteTAS():

	def raise_distance(data, outFolder, databasePath, comb, taName, numbercTAS, maxNumbercTAS, distTAS, startDist, maxDist, firstOc, saveFasta, ncTAS=0, itera=0):


		if (numbercTAS >= ncTAS) and (numbercTAS < maxNumbercTAS) and (distTAS <= maxDist) and (itera < 550):
			ncTAS = numbercTAS
			numbercTAS = FindCompleteTAS.find_tas(data, f'{outFolder}cTASv{distTAS}/', distTAS, databasePath, comb, taName, iteraV=True, saveCSV=False, saveFasta=False)
			if (numbercTAS > ncTAS):
				firstOc = distTAS
				itera = 0
				print(f'Number_of_cTAS:{numbercTAS}\nwith_maximal_distance_of:{distTAS}\n')

			else:
				itera += 1
			distTAS += 3
			FindCompleteTAS.raise_distance(data, outFolder, databasePath, comb, taName, numbercTAS, maxNumbercTAS, distTAS, startDist, maxDist, firstOc, saveFasta, ncTAS=ncTAS, itera=itera)
		else:
			# save cTAS to file
			numbercTASstart = FindCompleteTAS.find_tas(data, f'{fOutput}cTASv{startDist}/', startDist, dInput, comb, taName, saveFasta=saveFasta, saveCSV=True)
			if firstOc > startDist:
				if not os.path.exists(f'{fOutput}cTASv{firstOc}/'):
					os.system(f'mkdir -p {fOutput}cTASv{firstOc}/')

				numbercTAS = FindCompleteTAS.find_tas(data, f'{outFolder}cTASv{firstOc}/', firstOc, databasePath, comb, taName, saveCSV=True, saveFasta=saveFasta)
				print(f'max_Number_of_cTAS:{numbercTAS}\nmax_Number_of_cTAS_was_reached_with_maximal_distance_of:{firstOc}\n')
				# write cTAS to statistic
			else:
				print(f'max_Number_of_cTAS:{numbercTASstart}\nmax_Number_of_cTAS_was_reached_with_maximal_distance_of:{startDist}\n')



	def find_tas(data, outFolder, distanceTAS, databasePath, combBI, taName, iteraV=False, saveCSV=True, saveFasta=True):
		# data['TAsys'] = data['TAsys'].map(lambda x: x.rstrip('_round5_0001'))

		groupData = data.groupby(['bac', 'genome_chrom'])

		n = 0
		for key, group in groupData:
			# 	if len(groupData) > 1:
			if len(group) > 1:
			# for key, group in groupData:
				group = group.sort_values(by=['genomeStart'])
				subGroup = group[group.blastType != 'blastp'] #remove blastp hits, because it's found in proteome, not genome
				if len(subGroup) > 1:
					toCompare = subGroup.index.values.tolist()
					#####wenn man dist erhoeht, 'verschmelzen' die tas, die auf einem Genom gefunden werden
					overlappingRows = FindCompleteTAS.test_overlapp_tas(toCompare, subGroup, distanceTAS)
					nRows = FindCompleteTAS.merge(overlappingRows)

					a = 0
					genome = subGroup['bac'].iloc[0]
					if len(nRows) > 0:
						for combination in nRows:
							sGroup = group.loc[combination]

							################################# ich entferne blastp in zeile 203 und suche trotzdem danach ############################# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ###############

							checkSubTAS = [False, False, False, False]# at_nt, t_nt, t_pro_blastp, t_pro_tblastn/infernal

							tasDict = defaultdict(list) # Use the sequence as the key and then have a list of id's as the value

							for index, row in sGroup.iterrows():
								checkSubTAS[FindCompleteTAS.check_tas(row['TAsys'], row['blastType'])] = True
								tasDict[row['TAsys']+'_'+row['blastType']].append([row['strand'], row['genomeStart'], row['genomeStop']])

							if checkSubTAS[0] and (checkSubTAS[1] or checkSubTAS[3]): # if there is an AT and T (blastn or tblastn) -> cTAS is found

								for key in tasDict.keys(): ## are there multiple sequences for one Type (t / at) in one complete TAS ?
									if len(tasDict[key]) > 1:
										if iteraV:
											continue
										else:
											if saveCSV:
												print(f'{genome}: {tasDict}')


								#### check if cTAS is found with blast and infernal
								"""
								check if cTAS is found with both: blast and infernal

								group dataframe accoring to bac, genome_chrom, blastType???

								make a set each for blast and infernal  (t_nt, at_nt) == (t_nt, at_nt)
								compare both sets
								=> VERY STRICT

								get common elements of set -> 2 common both are found with both
								-> 1 only one is found with both
								-> 0 no element is found with both => do not write it to file

								"""

								if combBI:
									combInf = set()
									combBlast = set()
									for index, row in sGroup.iterrows():
										if row['blastType'] == 'infernal':
											combInf.add(row['TAsys'].split('_')[0])
										else:
											combBlast.add(row['TAsys'].split('_')[0])


									if len(combInf&combBlast) >= 1: #if len*combInf&combBlast) ==1 -> AT or T found with both; if ==2 -> AT and T are found with both

										n += 1
										a += 1
										if saveCSV:
											if not os.path.isfile(f'{outFolder}summary_cTAS.csv'):
												sGroup.to_csv(f'{outFolder}summary_cTAS.csv', sep='\t', header=False, index=False)
											else:
												sGroup.to_csv(f'{outFolder}summary_cTAS.csv', sep='\t', header=False, index=False, mode='a')
										if saveFasta:
											for index, row in sGroup.iterrows():
												OutputFastaTas.fasta(row, outFolder, databasePath, int(row['genomeStart']), int(row['genomeStop']), comb=True)

								else:
									n +=1 #count complete tas in summary
									a += 1 #count complete tas on the same genome
									if saveCSV:
										if not os.path.isfile(f'{outFolder}summary_cTAS.csv'):
											sGroup.to_csv(f'{outFolder}summary_cTAS.csv', sep='\t', header=False, index=False)
										else:
											sGroup.to_csv(f'{outFolder}summary_cTAS.csv', sep='\t', header=False, index=False, mode='a')
									if saveFasta:
										for index, row in sGroup.iterrows():
											OutputFastaTas.fasta(row, outFolder, databasePath, int(row['genomeStart']), int(row['genomeStop']))

					if a > 1:
						if saveCSV:
							print(f'{a} complete TAS are found on one genome: {genome}')

		return n

	def test_overlapp_tas(toCompare, group, dist):
		overlappingRows = []

		for pair in itertools.combinations(toCompare, 2):
			start1, stop1, start2, stop2 = FindCompleteTAS.start_stop(group, pair)
			if (start1 > stop1) and (start2 < stop2): # seq1 on minusStrand, seq2 on plusStrand
				if MergeSummary.overlaps(stop1, start1, start2, stop2, dist):
					overlappingRows.append([pair[0], pair[1]])
			elif (start2 > stop2) and (start1 < stop1): #seq2 on minusStrand, seq1 on plusStrand
				if MergeSummary.overlaps(start1, stop1, stop2, start2, dist):
					overlappingRows.append([pair[0], pair[1]])
			elif (start1 > stop1) and (start2 > stop2): #seq1 and seq2 on minusStrand
				if MergeSummary.overlaps(stop1, start1, stop2, start2, dist):
					overlappingRows.append([pair[0], pair[1]])
			else: #seq1 and seq2 on plusStrand
				if MergeSummary.overlaps(start1, stop1, start2, stop2, dist):
					overlappingRows.append([pair[0], pair[1]])

		return overlappingRows

	def start_stop(group, pair):
		start1 = group.loc[pair[0]]['genomeStart']
		stop1 = group.loc[pair[0]]['genomeStop']
		start2 = group.loc[pair[1]]['genomeStart']
		stop2 = group.loc[pair[1]]['genomeStop']

		return start1, stop1, start2, stop2

	def merge(overlappingRows):

		nRows = overlappingRows
		n = 0
		while n < len(nRows)-1:

			for pair in itertools.combinations(nRows, 2):
				oIndex = np.intersect1d(pair[0], pair[1]) #findet duplikate zwischen zwei listen
				if len(oIndex) > 0:  #es existieren duplikate
					n=0
					cRows = pair[0] + list(set(pair[1])-set(pair[0])) #entfernen der duplikate
					nRows[nRows.index(pair[0])] = cRows
					nRows.pop(nRows.index(pair[1]))
					break
				else:
					n += 1

		return nRows

	def check_tas(hit, blastType):

		if hit.split('_')[0] == 'at':
			return 0
		elif hit.split('_')[0] == 't':
			if hit.split('_')[1] == 'nt':
				return 1
			elif hit.split('_')[1] == 'pro':
				if blastType == 'blastp':
					return 2
				elif blastType == 'tblastn':
					return 3
				elif blastType == 'infernal':
					return 3
				else:
					print('ERROR: incorrect blastType')
					exit()
			else:
				print('ERROR: incorrect nt/pro type')
				exit()
		else:
			print('ERROR: incorrect at/t type')
			exit()



class OutputFastaTas():

	def fasta(group, outFolder, databasePath, start, stop, comb=False):
		tasID = group['tasID']
		bac = group['bac']
		blastType = group['blastType']
		genomeType = group['genome_chrom']
		genomeID = group['genomeID']
		taType = group['TAsys']
		strand = group['strand']
		genomeLength = group['genomeLength']

		dbfile = f'{databasePath}{OutputFasta.decide_file(blastType, genomeType)}'

		# OutputFasta.fastacmd(dbfile, genomeID, start, stop, 1, outFolder) if strand == '+' else OutputFasta.fastacmd(dbfile, genomeID, stop, start, 2, outFolder) if strand == '-' else print('Error: incorrect strand label') & exit()


		OutputFasta.fastacmd(dbfile, genomeID, start, stop, 1, outFolder) if (stop > start) else OutputFasta.fastacmd(dbfile, genomeID, stop, start, 2, outFolder) if (start > stop) else print('Error: incorrect strand label') & exit()

		newHeader = f'>{tasID}_bac:{bac}_tasys:{taType}_genome:{genomeID}_genomeType:{genomeType}_blastType:{blastType}_position:{start}to{stop}_genomeLength:{genomeLength}\n'
		if comb:
			# if (blastType == 'infernal') or (blastType == 'blastn'):
			saveFasta = f'{outFolder}/{taType}_ib_cTAS.fa'
			# else:

			# 	saveFasta = f'{outFolder}/{taType}_{blastType}_cTAS.fa'
		else:
			saveFasta = f'{outFolder}/{taType}_{blastType}_cTAS.fa'
		OutputFasta.write_fasta(saveFasta, newHeader, outFolder)


args = docopt("""merge_findctas.

Usage:
    merge_findctas.py [-h] [--s INT] [--sMax INT] [--rS] [--tas] [--m] [--t INT] [--comb] [--saveFasta BOOL] -d DIR -i FILE -o DIR -n STR

Options:
    -h --help           show this message and exit.
    -o DIR              Folder for output
    -i FILE             File with summary oder folder with multiple summarys
    -d DIR              database
    -n STR              name of TA System
    --s INT             (start) space between two TA-sequences [default: 0]
    --sMax INT          maximal space between two TA-seqences, only usable with parameter --rS [default: 2000]
    --rS                raise v automatic
    --tas               find complete Toxin-Antitoxin-Systems
    --m                 merge summary
    --comb              is the input a combined file between blast and infernal - output only cTAS that are found with both methods
    --saveFasta BOOL    True/False save found cTAS sequences as FastaFormat. [default: True]
    --t INT             number of reminders to give. [default: 1]

""")

if __name__ == '__main__':
	print(args)

	fInput = args['-i']
	fOutput = args['-o']
	dInput = args['-d']
	taName = args['-n']
	findTAS = args['--tas']
	mergeSum = args['--m']
	distTAS = int(args['--s'])
	maxDist = int(args['--sMax'])
	threads = int(args['--t'])
	comb = args['--comb']
	rV = args['--rS']
	saveFasta = StrToBool.str_to_bool(args['--saveFasta'])

	start = time.time()

	if not mergeSum and not findTAS:
		print('Error: please give at least one method (--m or --tas)')
		exit()

	if not os.path.exists(fOutput):
		os.system(f'mkdir -p {fOutput}')

	if os.path.isdir(fInput):
		mSummary = CombineSummarys.merge_summarys(fInput)
	elif os.path.isfile(fInput):
		mSummary = fInput
	else:
		print('Error: please insert Input-Folder or Input-File for merging - path not found')
		exit()


	sys.setrecursionlimit(10000)
	if mergeSum:
		start = time.time()
		if os.path.isfile(f'{fOutput}full_summary_merged.csv'):
			print(f'ERROR: merged output file already exists, {fOutput}full_summary_merged.csv')
			exit()
		else:
			dtype = [('tasID', 'U100'), ('bac', 'U100'), ('genome_chrom', 'U100') , ('strand', 'U100'), ('genomeID', 'U100'), ('TAsys', 'U100'), ('TAheader', 'U100'), ('blastType', 'U100'), ('eValue', 'U100'), ('genomeStart', int), ('genomeStop', int), ('alignLen', int), ('startTA', int), ('stopTA', int), ('lengthTA', int), ('genomeLength', int)]
			data = np.loadtxt(mSummary, delimiter='\t', dtype=dtype)
			MergeSummary.filter_summary(data, fOutput, dInput, threads, taName)
			end = time.time()
			print(f'time to merge: {end-start}')

	if findTAS:
		if os.path.isfile(f'{fOutput}summary_cTAS.csv'):
			print(f'ERROR: find-cTAS-output file already exists:{fOutput}summary_cTAS.csv')
			exit()

		startTAS = time.time()
		if mergeSum:
			data = pd.read_csv(f'{fOutput}full_summary_merged.csv', header=None, names=['tasID', 'bac', 'genome_chrom', 'strand', 'genomeID', 'TAsys', 'TAheader', 'blastType', 'eValue', 'genomeStart', 'genomeStop', 'alignLen', 'startTA', 'stopTA', 'lengthTA', 'genomeLength'], sep='\t', dtype={'eValue': object})

		else:
			data = pd.read_csv(mSummary, header=None, names=['tasID', 'bac', 'genome_chrom', 'strand', 'genomeID', 'TAsys', 'TAheader', 'blastType', 'eValue', 'genomeStart', 'genomeStop', 'alignLen', 'startTA', 'stopTA', 'lengthTA', 'genomeLength'], sep='\t', dtype={'eValue': object, 'genomeStart': int, 'genomeStop': int, 'startTA': int, 'stopTA': int, 'lengthTA': int, 'genomeLength': int})
		data = data.sort_values(by=['bac'])

		if rV:
			"""
			raise distance between T and AT (v) automatically

			cound number of at_nt and t_nt
			min number is number of max possible cTAS

			"""

			if not os.path.exists(f'{fOutput}cTASv{distTAS}/'):
				os.system(f'mkdir -p {fOutput}cTASv{distTAS}/')

			atNr = data.TAsys.str.count('^at_nt').sum()
			tNr = data.TAsys.str.count('^t_nt').sum()
			tproNr = data.apply(lambda x: True if 't_pro' in x['TAsys'] and x['blastType'] == 'tblastn' else False , axis=1)
			tproNr = len(tproNr[tproNr == True].index)

			maxcTASNr = min(atNr, tNr+tproNr)

			FindCompleteTAS.raise_distance(data=data, outFolder=fOutput, comb=comb, taName=taName, databasePath=dInput, numbercTAS=0, maxNumbercTAS=maxcTASNr, distTAS=distTAS, startDist=distTAS, maxDist=maxDist, firstOc=0, saveFasta=saveFasta)

			endTAS = time.time()
			print(f'time to find complete TAS: {endTAS-startTAS}')
		else:
			if not os.path.exists(f'{fOutput}cTASv{distTAS}/'):
				os.system(f'mkdir -p {fOutput}cTASv{distTAS}/')

			numbercTAS = FindCompleteTAS.find_tas(data, fOutput, distTAS, dInput, comb, taName, saveFasta=saveFasta)
			print(f'max_Number_of_cTAS:{numbercTAS}\nmax_Number_of_cTAS_was_reached_with_maximal_distance_of:{distTAS}\n')

			endTAS = time.time()
			print(f'time to find complete TAS: {endTAS-startTAS}')
