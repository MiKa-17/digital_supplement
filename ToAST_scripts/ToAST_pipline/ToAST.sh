#!/bin/bash

set -o errexit
set -o pipefail

cluster () {
	for file in $1*blastn*.fa; do
		fname=$(basename ${file##*/} .fa)
		cd-hit-est -i "$file" -o $2"$fname"_cdhit1.fa -c 1 >> $2log_"$fname"_cdhit1.txt
	done

	if [ -f $1*blastp*.fa ]; then
		for file in $1*blastp*.fa; do
	 		fname=$(basename ${file##*/} .fa)
	 		cd-hit -i "$file" -o $2"$fname"_cdhit1.fa -c 1 >> $2log_"$fname"_cdhit1.txt
		done
	fi

}

translate (){
	transeq -sequence $1*_tblastn_merged_cdhit1.fa -outseq $1t_pro_"$taName"_tblastn_cdhit1_protein.fa -clean -table 11 -sformat pearson >> $1temp.txt
	rm $1temp.txt
}

msa (){
	for file in $1*.fa; do
		if [ -f "$file" ]; then
			fname=$(basename ${file##*/} .fa)
			if [ $(grep '^>' $file -c) -ge 2 ]; then #more than one sequence in fasta-file
				echo "$fname"
				kalign -i $file -o $2"$fname".aln -f fasta -c tree -tgpe 0 -q #save in fasta-format and sort according to tree
				# muscle -in "$file" -out $2"$fname".aln -diags
			else
				cp $file $2"$fname".aln #just copy file with one sequence
			fi
		fi
	done
}

fold (){

	baseDir="$3"
	msaDir="$1"
	fileEnd="$2"
	newEnd="$4"

	scriptpath=$(pwd)/
	cd $3fold/
	for file in $1; do
		if [ -f "$file" ]; then
			fname=$(basename ${file##*/} $2)
			RNAalifold --color --id-prefix="$fname"_$4 --aln-stk="$fname"_$4 -q "$file"
		fi
	done
	cd $scriptpath


}

blast_pipeline () {

	prev=round"$(($r - 1))"
	round=round$r
	echo ">>>>>>>>>> start sequence analysis ""$round"
	# ta=$(basename $taDir)
	startSpace="$1"

	if [ $r = 1 ]; then
		toxinSeq="$tNtSeq"
		proteinSeq="$tProSeq"
		antitoxinSeq="$atNtSeq"
	else
		seqDir="$filterDir""$prev"/cluster/
		toxinSeq="$seqDir"t_nt*.fa
		antitoxinSeq="$seqDir"at_nt*.fa
		if [ -f "$seqDir"t_pro_"$taName"_blastp_tblastn_comb_cdhit1.fa ]; then
			proteinSeq="$seqDir"t_pro_"$taName"_blastp_tblastn_comb_cdhit1.fa
		elif [ -f "$seqDir"t_pro_"$taName"_blastp_merged_cdhit1.fa ]; then
			proteinSeq="$seqDir"t_pro_"$taName"_blastp_merged_cdhit1.fa
		elif [ -f "$seqDir"t_pro_"$taName"_tblastn_cdhit1_protein.fa ]; then
			proteinSeq="$seqDir"t_pro_"$taName"_tblastn_cdhit1_protein.fa
		else
			proteinSeq="no_file"
		fi
	fi

	echo ">>>>> start blast search"
	python3 box1_blast_search.py --tnt $toxinSeq --tpro "$proteinSeq" --atnt $antitoxinSeq --t "$threads" -o "$blastDir""$round"/ -d $bacdatabase -n "$taName"

	echo ">>>>> start to filter blast results"
	nice python3 box1_blast_filter.py -d $bacdatabase -i "$blastDir""$round"/ -o "$filterDir""$round"/ -n $taName --summary --evalue "$evalue" --pident "$pident"

	echo ">>>>> start to merge blast results at the same position and search for complete tas"
	mkdir -p "$taDir"1_sequence_analysis/cTAS/"$round"/
	#### merge summary
	nice python3 merge_findctas.py -i "$filterDir""$round"/ -o "$filterDir""$round"/ -d $bacdatabase --m --t 1 -n "$taName"
	#### find cTAS
	if [ "$recS" == true ]; then
		nice python3 merge_findctas.py -i "$filterDir""$round"/full_summary_merged.csv -o  "$taDir"1_sequence_analysis/cTAS/"$round"/ -d $bacdatabase --tas --t 1 -n "$taName" --rS --sMax "$maxS" --s "$startSpace" >> "$taDir"1_sequence_analysis/cTAS/"$round"/log_cTAS.txt
	else
		nice python3 merge_findctas.py -i "$filterDir""$round"/full_summary_merged.csv -o  "$taDir"1_sequence_analysis/cTAS/"$round"/ -d $bacdatabase --tas --t 1 -n "$taName" --s "$startSpace" > "$taDir"1_sequence_analysis/cTAS/"$round"/log_cTAS.txt
	fi

	for dir in "$taDir"1_sequence_analysis/cTAS/"$round"/*/; do
		if [ -z "$(ls -A "$dir")" ]; then
			rmdir $dir
		fi
	done

	### create Dir
	mkdir -p "$filterDir""$round"/cluster/ "$taDir"/1_sequence_analysis/fold/ "$filterDir""$round"/extend/cluster/ "$taDir"1_sequence_analysis/msa/"$round"/extend/ "$taDir"1_sequence_analysis/cTAS/"$round"/extend/

	echo ">>>>> start clustering the blast results."
	clusterDir="$filterDir""$round"/cluster/
	cluster $filterDir$round/ $clusterDir

	# combine blastp and tblastn search fasta:
	echo ">>>>> start combine blastp and tblastn search"

	if [ -f "$clusterDir"*_tblastn_merged_cdhit1.fa ]; then
		translate "$clusterDir"
		if [ -f "$clusterDir"*_blastp_merged_cdhit1.fa ]; then
			cat "$clusterDir"t_pro_"$taName"_tblastn_cdhit1_protein.fa "$clusterDir"*_blastp_merged_cdhit1.fa >> "$clusterDir"t_pro_"$taName"_blastp_tblastn_comb.fa
			cd-hit -i "$clusterDir"t_pro_"$taName"_blastp_tblastn_comb.fa -o "$clusterDir"t_pro_"$taName"_blastp_tblastn_comb_cdhit1.fa -c 1 >> "$clusterDir"t_pro_"$taName"_blastp_tblastn_comb_cdhit.txt
			rm "$clusterDir"t_pro_"$taName"_blastp_tblastn_comb.fa
		fi
	fi


	echo ">>>>> start msa with kalign"
	msa "$clusterDir" "$taDir"1_sequence_analysis/msa/"$round"/ #msa IN OUT

	echo ">>>>> start predict sec structure with RNAalifold"
	############# exclude comb results and tblastn_protein
	fold "$taDir"1_sequence_analysis/"msa/"$round"/*_nt_*blastn*.aln" "_blastn_merged_cdhit1.aln" "$taDir"1_sequence_analysis/ "$round"
	fold "$taDir"1_sequence_analysis/"msa/"$round"/t_pro*tblastn_merged_cdhit1.aln" "_tblastn_merged_cdhit1.aln" "$taDir"1_sequence_analysis/ "$round"

	echo ">>>>> start to extend shorter sequences"
	############# exclude comb results and tblastn_protein
	extendDir="$filterDir""$round"/extend/

	for msafile in "$taDir"1_sequence_analysis/msa/"$round"/*_nt_*blastn*.aln; do
		if [ -f "$msafile" ]; then
			python3 box1_extend.py -i $msafile -o "$extendDir" -n "$taName" -d $bacdatabase
		fi
	done
	for msafile in "$taDir"1_sequence_analysis/msa/"$round"/t_pro*tblastn_merged_cdhit1.aln; do
		if [ -f "$msafile" ]; then
			python3 box1_extend.py -i $msafile -o "$extendDir" -n "$taName" -d $bacdatabase
		fi
	done

	cluster "$extendDir" "$extendDir"cluster/

	echo '>>>>> start msa extend with kalign'
	msa "$extendDir"cluster/ "$taDir"1_sequence_analysis/msa/"$round"/extend/
	echo '>>>>> start predict extend sec structure with RNAalifold'

	fold "$taDir"1_sequence_analysis/msa/"$round"/extend/"*nt_*blastn*.aln" _blastn_extend_cdhit1.aln "$taDir"1_sequence_analysis/ "$round"_extend
	fold "$taDir"1_sequence_analysis/msa/"$round"/extend/"*tblastn*.aln" _tblastn_extend_cdhit1.aln "$taDir"1_sequence_analysis/ "$round"_extend

	if [ "$recS" == true ]; then
		nice python3 merge_findctas.py -i "$extendDir"extend_summary.csv -o "$taDir"1_sequence_analysis/cTAS/"$round"/extend/ -d "$bacdatabase" --tas -n "$taName" --rS --sMax "$maxS" --s "$startSpace" >> "$taDir"1_sequence_analysis/cTAS/"$round"/extend/log_cTAS.txt
	else
		nice python3 merge_findctas.py -i "$extendDir"extend_summary.csv -o "$taDir"1_sequence_analysis/cTAS/"$round"/extend/ -d "$bacdatabase" --tas -n "$taName" --s "$startSpace" >> "$taDir"1_sequence_analysis/cTAS/"$round"/extend/log_cTAS.txt
	fi

	for dir in "$taDir"1_sequence_analysis/cTAS/"$round"/extend/*/; do
		if [ -z "$(ls -A "$dir")" ]; then
			rmdir $dir
		fi
	done

}


process_results () {

	inDir="$1"
	clusterDir="$2"cluster/
	outDir="$3"
	foldPrefix="$4"


	for fasta in "$inDir"*.fa; do
		echo "$fasta"

		fname=$(basename ${fasta##*/} .fa)
		echo '>>> start clustering'
		cd-hit-est -i "$fasta" -o "$clusterDir""$fname"_cdhit1.fa -c 1 >> "$clusterDir"log_"$fname"_cdhit1.txt

		echo '>>> start kalign'
		if [ $(grep '^>' "$clusterDir""$fname"_cdhit1.fa -c) -ge 2 ]; then #more than one sequence in fasta-file
			kalign -i "$clusterDir""$fname"_cdhit1.fa -o "$outDir"msa/"$fname".aln -f fasta -c tree -tgpe 0 -quiet #save in fasta-format and sort according to tree
		else
			cp "$outDir"cluster/"$fname"_cdhit1.fa "$outDir"msa/"$fname".aln #just copy file with one sequence
		fi
		echo '>>> start kalign full'

		if [ $(grep '^>' "$fasta" -c) -ge 2 ]; then #more than one sequence in fasta-file
			kalign -i "$fasta" -o "$outDir"msa/"$fname"_full.aln -f fasta -c tree -tgpe 0 -quiet #save in fasta-format and sort according to tree
		else
			cp $fasta "$outDir"msa/"$fname".aln #just copy file with one sequence
		fi

	done

	echo '>>> start folding'
	fold "$outDir""msa/*.aln" .aln "$outDir" "$foldPrefix"

}

infernal_pipeleine () {

	echo '>>>>> start infernal'

	takeExtend="$1"
	startSpace="$2"
	infRound="$3"

	if [ "$takeExtend" == false ]; then
		strucDir="$taDir"2_structure_analysis/
		strucFiles="$taDir"1_sequence_analysis/fold/*round"$infRound".stk
		fileEnd=_round"$infRound".stk
	else
		strucDir="$taDir"extend_analysis_2_3_4/2_structure_analysis/
		strucFiles="$taDir"1_sequence_analysis/fold/*round"$infRound"_extend.stk
		fileEnd=_round"$infRound"_extend.stk
	fi

	searchDir="$strucDir"infernal_search/
	filterDir="$strucDir"infernal_filter/
	mkdir -p "$filterDir"cluster/ "$strucDir"msa/ "$strucDir"fold/ "$strucDir"cTAS/ "$searchDir"

	for file in $strucFiles; do
		f=$(basename $file)
		# ff="${f%.*}"
		ff=$(basename ${file##*/} "$fileEnd")
		infernal "$searchDir""$ff" "$file" "$bacdatabase"
		nice python3 box2_infernal_filter.py -i "$searchDir""$ff"_cmsearch_chrom.txt,"$searchDir""$ff"_cmsearch_genome.txt -d "$bacdatabase" -o "$filterDir" -n "$taName" -l chromosome,genome
		echo 'infernal done with '"$file"
	done

	echo '>>>>> start cluster infernal results with cdhit and get msa with kalign and predict secundary structure'
	process_results "$filterDir" "$filterDir" "$strucDir" round"$infRound" # IN OUT SearchFor roundNr
	echo '>>>>> start find cTAS of Infernal'

	if [ "$recS" == true ]; then
		nice python3 merge_findctas.py -i "$filterDir"infernal_summary.csv -o "$strucDir"cTAS/ -d $bacdatabase --tas -n $taName --rS --sMax "$maxS" --s "$startSpace" >> "$strucDir"cTAS/log_cTAS.txt
	else
		nice python3 merge_findctas.py -i "$filterDir"infernal_summary.csv -o "$strucDir"cTAS/ -d $bacdatabase --tas -n $taName --s "$startSpace"  >> "$strucDir"cTAS/log_cTAS.txt
	fi

	for dir in "$strucDir"cTAS/*/; do
		if [ -z "$(ls -A "$dir")" ]; then
			rmdir $dir
		fi
	done
}

infernal () {
	nice cmbuild $1_cmbuild.txt $2

	nice cmcalibrate --cpu "$threads" $1_cmbuild.txt
	# notextw -> genom namen bleiben erhalten und werden nicht abgeschnitten
	nice cmsearch -o $1_cmsearch_chrom.txt --cpu "$threads" --notextw $1_cmbuild.txt $3chromosome/full_chromosome_bacteria.fna

	nice cmsearch -o $1_cmsearch_genome.txt --cpu "$threads" --notextw $1_cmbuild.txt $3genome/full_genome_bacteria.fna
}


combine_blast_infernal () {

	takeExtend="$1"
	startSpace="$2"
	blastRound="$3"

	if [ "$takeExtend" == false ]; then
		combDir="$taDir"3_combined_analysis/
		seqFile="$taDir"1_sequence_analysis/blast_filter/round"$blastRound"/full_summary_merged.csv
		strucDir="$taDir"2_structure_analysis/infernal_filter/
	else
		combDir="$taDir"extend_analysis_2_3_4/3_combined_analysis/
		strucDir="$taDir"extend_analysis_2_3_4/2_structure_analysis/infernal_filter/
		seqFile="$taDir"1_sequence_analysis/blast_filter/round"$blastRound"/extend/extend_summary.csv

	fi

	echo '>>>>> start combine results of sequence and structure search'
	mkdir -p "$combDir"cTAS/

	## combine summarys of blast and infernal
	if ! [ -f "$combDir"blast_infernal_summary.csv ]; then
		cat "$seqFile" "$strucDir"infernal_summary.csv >> "$combDir"blast_infernal_summary.csv
	fi

	### find cTAS in combined data
	echo '>>>>> start find combined cTAS of Infernal and Blast'


	if [ "$recS" == true ]; then
		python3 merge_findctas.py -i "$combDir"blast_infernal_summary.csv -o "$combDir"cTAS/ -n $taName --tas --comb --s "$startSpace" -d $bacdatabase --rS --sMax "$maxS" >> "$combDir"cTAS/log_cTAS.txt
	else
		python3 merge_findctas.py -i "$combDir"blast_infernal_summary.csv -o "$combDir"cTAS/ -n $taName --tas --comb --s "$startSpace" -d $bacdatabase >> "$combDir"cTAS/log_cTAS.txt
	fi

	### cluster and msa

	for dir in "$combDir"cTAS/*/; do
		if [ -z "$(ls -A "$dir")" ]; then
			rmdir $dir
		fi
	done

	subdircount=$(find "$combDir"cTAS/ -maxdepth 1 -type d | wc -l)

	if [ $subdircount -eq 1 ]; then
		echo 'no cTAS in combined analysis'
		echo -e 'no cTAS in combined analysis\nmax_Number_of_cTAS:0\nmax_Number_of_cTAS_was_reached_with_maximal_distance_of:'"$maxS" >> "$combDir"cTAS/log_cTAS.txt

	else
		maxDist=$(ls -d "$combDir"cTAS/*/ | tail -1)

		mkdir -p "$combDir"cluster/ "$combDir"msa/ "$combDir"fold/
		echo '>>>>> start only take sequences with most cTAS for cluster, msa and predict secundary structure'
		##### only take sequences with most cTAS for cluster, msa and folding
		process_results "$maxDist" "$combDir" "$combDir" "comb"

	fi


}

find_promoter () {

	takeExtend="$1"
	startSpace="$2"

	if [ "$takeExtend" == false ]; then
		combDir="$taDir"3_combined_analysis/
		promDir="$taDir"4_promoter_analysis/

	else
		combDir="$taDir"extend_analysis_2_3_4/3_combined_analysis/
		promDir="$taDir"extend_analysis_2_3_4/4_promoter_analysis/

	fi

	echo '>>>>> start find promoter motives in combined cTAS'

	dir=$(ls -d "$combDir"cTAS/*/ | tail -1)

	mkdir -p "$promDir" "$promDir"/cTAS/
	echo $dir
	echo $promDir
	python3 box4_prom.py -i "$dir"summary_cTAS.csv -o "$promDir" -d $bacdatabase

	##### cTAS (without promoter!!!!)
	echo '>>>>> start find combined cTAS of comb Hits with prom'

	if [ "$recS" == true ]; then
		python3 merge_findctas.py -i "$promDir"promoter_summary.csv -o "$promDir"cTAS/ -n "$taName" -d "$bacdatabase" --tas --comb --s "$startSpace"  --rS --sMax "$maxS" >> "$promDir"/cTAS/log_cTAS.txt
	else
		python3 merge_findctas.py -i "$promDir"promoter_summary.csv -o "$promDir"cTAS/ -n "$taName" -d "$bacdatabase" --tas --comb --s "$startSpace" >> "$promDir"/cTAS/log_cTAS.txt
	fi

	## cluster and msa
	echo '>>>>> start only take sequences with most cTAS for cluster, msa and predict secundary structure'

	maxDist=$(ls -d "$promDir"cTAS/*/ | tail -1)
	echo "$maxDist"

	mkdir -p "$promDir"cluster/ "$promDir"msa/ "$promDir"fold/

	##### only take sequences with most cTAS for cluster, msa and folding
	process_results "$maxDist" "$promDir" "$promDir" "comb"

}

bacteria () {

	takeExtend="$1"
	startSpace="$2"

	if [ "$takeExtend" == false ]; then
		bacDir="$taDir"5_bacteria_with_cTAS/
		combDir="$taDir"3_combined_analysis/
	else
		combDir="$taDir"extend_analysis_2_3_4/3_combined_analysis/
		bacDir="$taDir"extend_analysis_2_3_4/5_bacteria_with_cTAS/

	fi
	mkdir -p "$bacDir"
	python3 box5_bac.py -i "$combDir" -o "$bacDir" -d "$bacdatabase"


	for bac in "${bacDir[@]}"*/; do
		mkdir -p "$bac"/cluster/ "$bac"/msa/ "$bac"fold/ "$bac"/cTAS/
		process_results "$bac" "$bac" "$bac"
		if [ "$recS" == true ]; then
			python3 merge_findctas.py -i "$bac"*_summary.csv -o "$bac"cTAS/ -n "$taName" -d "$bacdatabase" --tas --comb --s "$startSpace"  --rS --sMax "$maxS" >> "$bac"cTAS/log_cTAS.txt
		else
			python3 merge_findctas.py -i "$bac"*_summary.csv -o "$bac"cTAS/ -n "$taName" -d "$bacdatabase" --tas --comb --s "$startSpace" >> "$bac"cTAS/log_cTAS.txt
		fi
	done

}


statistic () {

	takeExtend="$1"
	combcTAS="$2"

	if [ "$takeExtend" == false ]; then
		baseDir="$taDir"
		statName=statistic
 		strucQuery='structure of round'"$roundNr"
 		combQuery='BlastRound'"$roundNr"'&Infernal'

		# seqFile="$taDir"1_sequence_analysis/blast_filter/round"$roundNr"/full_summary_merged.csv

	else
		baseDir="$taDir"extend_analysis_2_3_4/
		statName='extend statistic'
		strucQuery='extend structure of round'"$roundNr"
		combQuery='ExtendBlastRound'"$roundNr"'&Infernal'
		# seqFile="$taDir"1_sequence_analysis/blast_filter/round"$roundNr"/extend/extend_summary.csv
	fi

	strucDir="$baseDir"2_structure_analysis/
	combDir="$baseDir"3_combined_analysis/
	promDir="$baseDir"4_promoter_analysis/
	bacDir="$baseDir"5_bacteria_with_cTAS/

	echo '>>>>> start writing statistic'
	##### TA Name and header
	echo -e '>>>>> '"$statName"' <<<<<\tNameOfTAS:\t'"$taName"'\t\nmethod\tRound/Query\tNrOfCompleteTAS\tMaxDistanceForCTAS\tAtSeq\tAtSeqUnique\tTntSeq\tTntSeqUnique\tTpro(blastp)\tTpro(blastp)Unique\tTpro(tblastn)\tTpro(tblastn)Unique\tTpro(blastp&tblastn)Unique'>> "$taDir"statistic.csv

	##### 1. sequence analyis

 	echo '>>>>> start writing seq statistic'
	for ((r=1; r<=$roundNr; r++)); do
		echo -e '1. sequence analysis\tround'"$r"'\t' | tr -d '\n' >> "$taDir"statistic.csv
		# echo -e "$1\t" | tr -d '\n' >> "$taDir"statistic.csv
	 	if [ "$takeExtend" == false ]; then
	 		tasDir="$taDir"1_sequence_analysis/cTAS/round"$r"/
	 		secDir="$taDir"1_sequence_analysis/blast_filter/round"$r"/
			# seqFile="$taDir"1_sequence_analysis/blast_filter/round"$roundNr"/full_summary_merged.csv

		else
			tasDir="$taDir"1_sequence_analysis/cTAS/round"$r"/extend/
			secDir="$taDir"1_sequence_analysis/blast_filter/round"$r"/extend/
			# seqFile="$taDir"1_sequence_analysis/blast_filter/round"$roundNr"/extend/extend_summary.csv
		fi

		grep_ctas_number "$tasDir"

	 	fileNames=("at*.fa" "t_nt*.fa" "t_pro*blastp_merged*.fa" "t_pro*tblastn_merged*.fa" "t_pro*blastp_tblastn_comb*.fa")
	 	for file in "${fileNames[@]}"; do
	 		grep_seq_number "$secDir" $file "$secDir"
	 	done
	 	echo -e '\t' >> "$taDir"statistic.csv
	done


 	### 2. structure analysis
 	echo '>>>>> start writing stuct statistic'
	echo -e '2. structure analysis\t'"$strucQuery"'\t' | tr -d '\n' >> "$taDir"statistic.csv
 	grep_ctas_number "$strucDir"cTAS/
 	echo -e '---\t---\t---\t---\t---\t---\t---\t---\t---\t' >> "$taDir"statistic.csv
	structure_statistic "at*cmsearch_genome.txt" "at*infernal*.fa" 2 7  "$strucDir"
	structure_statistic "t_nt*cmsearch_genome.txt" "t_nt*infernal*.fa" 4 5  "$strucDir"
	structure_statistic "t_pro*cmsearch_genome.txt" "t_pro*infernal*.fa" 8 1  "$strucDir"

 	##### 3. combined analysis
 	echo '>>>>> start writing comb statistic'
	echo -e '3. combined analysis\t'"$combQuery"'\t' | tr -d '\n' >> "$taDir"statistic.csv
	grep_ctas_number "$combDir"cTAS/

	if [ "$combcTAS" == true ]; then
	 	dir=$(ls -d "$combDir"cTAS/*/ | tail -1)
	 	combFiles=("at*ib*.fa" "t_nt*ib*.fa" "NoFile" "t_pro*ib*.fa" "t_pro*blastp_tblastn_comb*.fa")
	 	for combFile in "${combFiles[@]}"; do
			grep_seq_number "$dir" $combFile "$combDir"
	 	done
	 	echo -e '\t' >> "$taDir"statistic.csv

		##### 4. promoter analysis

		echo '>>>>> start writing prom statistic'
		echo -e '4. promoter analysis (exact match)\tcTAS seq of combined analysis\t' | tr -d '\n' >> "$taDir"statistic.csv

		grep_ctas_number "$promDir"cTAS/

		promKind=('prom_complete' 'prom_tata_sd_start' 'prom_35_tata_start' 'prom_35_sd_start' 'prom_sd_start' 'prom_tata_start' 'prom_35_tata' 'prom_start')

		for prom in "${promKind[@]}"; do

			promFiles=("at_nt_ib_""$prom"".fa" "t_nt_ib_""$prom"".fa" "NoFile" "t_pro_ib_""$prom"".fa" "LastFile")

			echo -e '\t' >> "$taDir"statistic.csv
			echo -e '4. promoter analysis (exact match)\t'"$prom"'\t' | tr -d '\n' >> "$taDir"statistic.csv
			echo -e "--\t--\t" | tr -d '\n' >> "$taDir"statistic.csv

			for promFile in "${promFiles[@]}"; do
				grep_seq_number_prom "$promDir" "$promFile"
			done

		done


		echo -e '\t\n' >> "$taDir"statistic.csv

		##### 5.
		## in which bacteria are cTAS found
		echo '>>>>> start writing bac statistic'
		echo -e '5. cTAS are found in: (query: cTAS seq of combined analysis)' | tr -d '\n' >> "$taDir"statistic.csv

		bacFiles=("at_nt*.fa" "t_nt*.fa" "NoFile" "t_pro*.fa" "t_pro*blastp_tblastn_comb*.fa")
		for bac in "${bacDir[@]}"*/; do
			bacName=$(basename $bac)
			echo -e '\t'"$bacName"'\t'| tr -d '\n' >> "$taDir"statistic.csv
			grep_ctas_number "$bac"cTAS/
			for bacFile in "${bacFiles[@]}"; do
				grep_seq_number "$bac" $bacFile "$bac"
			done
			echo -e '' >> "$taDir"statistic.csv
		done

		echo -e 'multiple cTAS are found in: (query: cTAS seq of combined analysis)' | tr -d '\n' >> "$taDir"statistic.csv


		if [ $(grep 'found on one' "$bac"cTAS/log_cTAS.txt | wc -l) -eq 0 ]; then
			echo -e '\tNo multiple cTAS on one Genome\t' >> "$taDir"statistic.csv
		else
			grep 'found on one' "$bac"cTAS/log_cTAS.txt | sort | while read -r line; do
				IFS=' '
				read -ra ADDR <<< $line
				numberTAS=${ADDR[0]}
				nameTAS=${ADDR[-1]}
				echo -e '\t'"$nameTAS"'\t'"$numberTAS"'\t' >> "$taDir"statistic.csv
			done
		fi
	else
		echo -e '\nno cTAS in combined analysis leads to no statistic for promoter and bacteria analysis' >> "$taDir"statistic.csv
	fi


}

grep_ctas_number () {
	grep 'max_Number_of_cTAS:' "$1"log_cTAS.txt | awk -F ":" '{print $2}' | tr -d '\n' >> "$taDir"statistic.csv && echo -e "\t" | tr -d '\n' >> "$taDir"statistic.csv
	grep 'max_Number_of_cTAS_was_reached_with' "$1"log_cTAS.txt | awk -F ":" '{print $2}' | tr -d '\n' >> "$taDir"statistic.csv && echo -e "\t" | tr -d '\n' >> "$taDir"statistic.csv
}

grep_seq_number () {

	tasDir="$1"
	inFile=$2
	baseDir="$3"

	# __base=$(echo $inFile | awk -F "*" '{print $inFile}')

	if [[ "$inFile" == "t_pro*blastp_tblastn_comb*.fa" ]]; then
		if [ -f "$baseDir"cluster/$inFile ]; then
			grep '^>' "$baseDir"cluster/$inFile -c | tr -d '\n' >> "$taDir"statistic.csv && echo -e "\t" | tr -d '\n' >> "$taDir"statistic.csv
		else
			echo -e "0\t" | tr -d '\n' >> "$taDir"statistic.csv
		fi
	else
		if [ -f "$tasDir"$inFile ]; then
			grep '^>' "$tasDir"$inFile -c | tr -d '\n' >> "$taDir"statistic.csv && echo -e "\t" | tr -d '\n' >> "$taDir"statistic.csv
			grep '^>' "$baseDir"cluster/$inFile -c | tr -d '\n' >> "$taDir"statistic.csv && echo -e "\t" | tr -d '\n' >> "$taDir"statistic.csv
		else
			echo -e "0\t0\t" | tr -d '\n' >> "$taDir"statistic.csv
		fi
	fi
}

grep_seq_number_prom () {

	promDir="$1"
	inFile=$2

	if [[ "$inFile" == "LastFile" ]]; then
		echo -e "--\t" | tr -d '\n' >> "$taDir"statistic.csv
	elif [[ "$inFile" == "NoFile" ]]; then
		echo -e "--\t--\t" | tr -d '\n' >> "$taDir"statistic.csv
	else
		if [ -f "$promDir"$inFile ]; then
			grep '^>' "$promDir"$inFile -c | tr -d '\n' >> "$taDir"statistic.csv
			echo -e "\t--\t" | tr -d '\n' >> "$taDir"statistic.csv
		else
			echo -e "0\t--\t" | tr -d '\n' >> "$taDir"statistic.csv
		fi
	fi
}

structure_statistic () {

	if [ -f "$5"infernal_search/$1 ]; then
		echo -e '2. structure analysis\t' | tr -d '\n' >> "$taDir"statistic.csv
		grep '^Query:' "$5"infernal_search/$1 | awk -F " " '{print $2}' | tr -d '\n' >> "$taDir"statistic.csv
		echo -e '\t' | tr -d '\n' >> "$taDir"statistic.csv
		mySting=$(printf "%$3s"); echo -e ${mySting// /'---\t'} | tr -d '\n' >> "$taDir"statistic.csv
		grep_seq_number "$5"infernal_filter/ $2 "$5"infernal_filter/
		mySting=$(printf "%$4s"); echo -e ${mySting// /'---\t'} | tr -d '\n' >> "$taDir"statistic.csv
		echo -e '\t' >> "$taDir"statistic.csv
	fi

}


function usage() {
	echo ""
	echo "Pipeline to detect Toxin-Antitoxin-Systems. Using sequence and structure analysis and promoter detection."
	echo ""
	echo "Usage:"
	echo "	./pipeline_tas.sh [optional_arguments] -n=STR -t=FILE -p=FILE -a=FILE "
	echo ""
	echo "-h | --help			show usage and exit."
	echo ""
	echo ">>>>> Required arguments:"
	echo "-n=STR | --name=STR 		Name of the TA-System (toxin_antitoxin)"
	echo "-t=FILE | --tNt=FILE		Fasta-File with Toxin nucleotide sequence"
	echo "-p=FILE | --tPr=FILE 		Fasta-File with Toxin aminoacid sequence"
	echo "-a=FILE | --aNt=FILE		Fasta-File with Antitoxin nucleotied sequence"
	echo "-d=STR | --db=STR 		Directory-Path to the bacteria-database (blast-db) "
	echo ""
	echo ">>>>> Optional arguments: "
	echo "-o=STR | --out=STR 		Directory-Path to save output [default: current dir]"
	echo "-s=INT | --subp=INT 		Number of subProgams to be execut [default: 12345]"
	echo "					1 = Seqeuence analysis (Blast)"
	echo "					2 = Structrue analysis (Infernal)"
	echo "					3 = Combined analysis"
	echo "					4 = Promoter analysis"
	echo "					5 = Summary writing"
	echo "-e=INT | --evalue=INT		Evalue-threshold to filter the results of blast [default: 0.0000001]"
	echo "-i=INT | --pident=INT		Evalue-threshold to filter the results of blast [default: 0.0000001]"
	echo "-c=INT | --cores=INT		number of used threads for blast [default: 1]"
	echo "-r=INT | --rounds=INT		number of blast iterations [default: 5]"
	echo "-b=INT | --bround=INT		take sec. structrue for infernal of this round of blast [default: highest number of blast iterations]"
	echo "--extend			use extend structure for infernal [default: false]"
	echo "--recS				use rising Space between T and AT to predict cT [default: false]"
	echo "--startS=INT			Space between T and AT to predict cTAS"
	echo "--maxS=INT			maximal permitted Space between T and AT to predict cTAS"
	echo ""
}


while [ "$1" != "" ]; do
	PARAM=`echo $1 | awk -F= '{print $1}'`
	VALUE=`echo $1 | awk -F= '{print $2}'`
	case $PARAM in
		-h | --help) usage && exit;;
		-n | --name) tasName=$VALUE;;
		-t | --tNt) Tnt_SEQ=$VALUE;;
		-p | --tPr) Tpro_SEQ=$VALUE;;
		-a | --aNt) Ant_SEQ=$VALUE;;
		-d | --db) bacdatabase=$VALUE;;
		-o | --out) outDir=$VALUE;;
		-s | --subp) subp=$VALUE;;
		-e | --evalue) evalue=$VALUE;;
		-i | --pident) pident=$VALUE;;
		-r | --rounds) rounds=$VALUE;;
		-b | --bround) blastRound=$VALUE;;
		--startS) startSpace=$VALUE;;
		--maxS) maxSpace=$VALUE;;
		--recS) recS=true;;
		--extend) extend=true;;
		-c | --cores) threads=$VALUE;;
		-* | --*= | --*=*--*) # unsupported flags
			echo "Error: Unsupported flag \"$PARAM\""
			usage
			exit 1
			;;
		*)
			echo "ERROR: unknown parameter \"$PARAM\""
			usage
			exit 1
			;;
	esac
	shift
done

if [ ! -f "$Tnt_SEQ" ] && [ ! -f "$Tpro_SEQ" ] && [ ! -f "$Ant_SEQ" ]; then
	echo ""
	echo "ERROR: Could not find sequence files; Please, give at least one Toxin-Antitoxin-Sequence"
	usage
	exit
fi
if [ ! -d "$bacdatabase" ]; then
	echo ""
	echo "ERROR: Could not find diretory for bacteria database"
	usage
	exit
fi
# bacdatabase=${dbPath:-"/mnt/fass1/genomes/new_bacteria/bacteria_blast_db/"}

#### setting default parameter ####

tNtSeq=${Tnt_SEQ:-"CAUTION:NoToxinNtSequence"}
tProSeq=${Tpro_SEQ:-"CAUTION:NoToxinProteinSequence"}
atNtSeq=${Ant_SEQ:-"CAUTION:NoAntitoxinNtSequence"}
taName=${tasName:-"toast_output"} ################ geht das mit den INPUT NAMEN ????????????????????????????????????/
outDir=${outDir:-"$(pwd)"}
subProg=${subp:-123456}
evalue=${evalue:-0.0000001}
pident=${pident:-0.0}
threads=${threads:-1}
roundNr=${rounds:-5}
blastRound=${blastRound:-$roundNr}
extend=${extend:-false}
recS=${recS:-false}
startS=${startSpace:-0}
maxS=${maxSpace:-300}

#### write the choosen parameters ####

echo ">>>>> start to search for toxin antitoxin systems with these arguments:"
echo "Fasta-File with Toxin nucleotide sequence-------------------"$tNtSeq
echo "Fasta-File with Toxin aminoacid sequence--------------------"$tProSeq
echo "Fasta-File with Antitoxin nucleotied sequence---------------"$atNtSeq
echo "Name of the TA-System---------------------------------------"$taName
echo "Directory-Path to save output-------------------------------"$outDir
echo "Directory-Path to the bacteria-database---------------------"$bacdatabase
echo "Number of subProgams to be executed-------------------------"$subProg
echo "Evalue-threshold to filter the results of blast-------------"$evalue
echo "percent identity threshold to filter the results of blast---"$pident
echo "number of used threads for blast----------------------------"$threads
echo "number of blast iterations----------------------------------"$roundNr
echo "take structure of round x for infernal----------------------"$blastRound
echo "use extend structure for infernal---------------------------"$extend
echo "use rising Space between T and AT to predict cTAS-----------"$recS
echo "start Space between T and AT to predict cTAS----------------"$startS
echo "maximal permitted Space between T and AT to predict cTAS:---"$maxS
echo -e "\n"


taDir=$(readlink -f "$outDir")/"$tasName"/
mkdir -p "$taDir"

totalStart=$SECONDS

###################### Sequence analysis ######################
if [[ " ${subProg[*]} " == *"1"* ]]; then
	START=$SECONDS
	echo "############################################ start to analyse according to sequence"
	blastDir="$taDir"1_sequence_analysis/blast_search/
	filterDir="$taDir"1_sequence_analysis/blast_filter/
	mkdir -p "$filterDir"
	for ((r=1; r<=$roundNr; r++)); do
		START_TIME=$SECONDS
		blast_pipeline "$startS"
		ELAPSED_TIME=$(($SECONDS - $START_TIME))
		echo "###################### finish $round runtime for sequence analysis: $ELAPSED_TIME "
	done
	END=$(($SECONDS - $START))
	echo "############################################ finish sequence analysis runtime: ""$END"
	echo "#"
else
	echo "############################################ CAUTION! Sequence-search and all processing steps are skiped"
fi
echo '#'

###################### Structure analysis ######################
if [[ " ${subProg[*]} " == *"2"* ]]; then
	startInfernal=$SECONDS
	echo "############################################ start to analyse according to structure "
	infernal_pipeleine "$extend" "$startS" "$blastRound"
	endInfernal=$(($SECONDS - $startInfernal))
	echo "############################################ finish structrue analysis runtime: ""$endInfernal"
else
	echo "############################################ CAUTION! Structure search and all processing steps are skiped"
fi


# if [[ " ${subProg[*]} " == *"2"* ]]; then
# 	echo "############################################ start to analyse according to structure "
# 	startInfernal=$SECONDS
# 	infernal_pipeleine "$extend" "$startS" "$blastRound"
# 	endInfernal=$(($SECONDS - $startInfernal))
# 	echo "############################################ finish structrue analysis runtime: ""$endInfernal"
# 	echo "#"
# else
# 	echo "############################################ CAUTION! Structure-search and all processing steps are skiped"
# fi
# echo "#"

###################### Combined analysis ######################
if [[ " ${subProg[*]} " == *"3"* ]]; then
	echo "############################################ start combine the results of structure and sequence results "
	combine_blast_infernal "$extend" "$startS" "$blastRound"
else
	echo "############################################ CAUTION! Combined analysis and all processing steps are skiped"
fi

if [ "$extend" == false ]; then
	combDir="$taDir"3_combined_analysis/
else
	combDir="$taDir"extend_analysis_2_3_4/3_combined_analysis/
fi

subdircount=$(find "$combDir"cTAS/ -maxdepth 1 -type d | wc -l)

if [ $subdircount -eq 1 ]; then
	echo -e "no cTAS in combined analysis\npromoter analysis (4) and bacteria analysis (5) are skiped\ntry to change parameter for better results\ne.g e-value for blast filter or maximal distance between toxin/antitoxin sequences"
	combcTAS=false
else
	###################### Promoter analysis ######################
	if [[ " ${subProg[*]} " == *"4"* ]]; then
		echo "############################################ start to analyse promoter in combineds complete TAS "
		find_promoter "$extend" "$startS"
	else
		echo "############################################ CAUTION! promoter analysis and all processing steps are skiped"

	fi

	###################### Bacteria ######################
	if [[ " ${subProg[*]} " == *"5"* ]]; then
		echo "############################################ start to analyse in which bacteria cTAS are found "
		bacteria "$extend" "$startS"
	else
		echo "############################################ CAUTION! Show in witch bacteria cTAS are found is skiped"

	fi
	combcTAS=true
fi

###################### Statistic ######################
if [[ " ${subProg[*]} " == *"6"* ]]; then
	echo "############################################ start to write statistic"
	statistic "$extend" "$combcTAS"
else
	echo "############################################ CAUTION! writing the statistic is skiped"

fi

totalEnd=$(($SECONDS - $totalStart))

echo "####################################################### finish complete TAS analysis runtime: ""$totalEnd"
