#!/bin/bash

set -o errexit
set -o pipefail



bacteria () {

	startSpace="$2"
	recS="$1"
	bacDir="$3"
	maxS="$4"
	bacdatabase="$5"
	recS="$6"

	for bac in "${bacDir[@]}"*/; do
		echo $bac
		echo $bacDir

		mkdir -p "$bac"/cluster/ "$bac"/msa/ "$bac"fold/ "$bac"/cTAS/
		process_results "$bac" "$bac" "$bac"
		if [ "$recS" == true ]; then
			python3 merge_findctas.py -i "$bac"*_summary.csv -o "$bac"cTAS/ -n bacteria -d "$bacdatabase" --tas --comb --s "$startSpace" --rS --sMax "$maxS" >> "$bac"cTAS/log_cTAS.txt
		else
			python3 merge_findctas.py -i "$bac"*_summary.csv -o "$bac"cTAS/ -n bacteria -d "$bacdatabase" --tas --comb --s "$startSpace" >> "$bac"cTAS/log_cTAS.txt
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

tasAll=(aapa_isoa/aapa_isoa_sr lpt_rna2 zoro_orzo/zoro_orzo_sr cjpt_cjra yont_sr6 cds_rcd dinq_agrb/dinq_agrb_sr fst_rna1 hok_sok_ena hok_sok_tadb/hok_sok_tadb_sr ibs_sib ldr_rdl shob_ohsc/shob_ohsc_sr tisb_istr/tisb_istr_sr txpa_rata_bacillus/txpa_rata_bacillus_sr txpa_rata_enterococcus ralr_rala srnb_srnc syme_symr sprg_sprf/sprg_sprf_sr75 pnda_pndb)


for tas in "${tasAll[@]}"; do
 	echo "start with filter bacteria in ""$tas"
  	baseDir=$(readlink -f "../results/")
  	taDir="$baseDir"/ta_systems/"$tas"/
  	bacDir="$baseDir"/bacteria/
  	mkdir -p "$bacDir"
  	combDir="$taDir""3_combined_analysis/"
  	python3 box5_bac.py -i "$combDir" -o "$bacDir" -d "/mnt/fass1/genomes/new_bacteria/bacteria_blast_db/"
done

echo ">>>>>>>>start_processing"

bacDir=$(readlink -f "../results/bacteria/")


bacteria false 0 "$bacDir"/ 350 /mnt/fass1/genomes/new_bacteria/bacteria_blast_db/ true
