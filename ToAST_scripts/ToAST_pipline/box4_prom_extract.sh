

detect_prom () {
	taDir=/home/co78wid/projects/master/ta_systems/$1/
	promDir="$taDir"4_promoter_analysis/

	mkdir -p "$promDir"

	python3 box4_prom.py -i "$taDir"combine_bi/cTASv$2/full_summary_merged_complete_tas.csv -o "$promDir" -d $bacdatabase >> "$promDir"log_prom.txt



}

### all tas with cTAS with v=0
# tas=(cds2517_rcd8 ralr_rala syme_symr)
# # yont_sr6 cds2517_rcd8 ralr_rala syme_symr aapa_isoa)

bacdatabase=/mnt/fass1/genomes/new_bacteria/bacteria_blast_db/

# for taName in "${tas[@]}"; do
# 	echo $taName
# 	detect_prom $taName 0
# done

# detect_prom dinq_agrb 241

# detect_prom pnda_pndb 49

# detect_prom ldr_rdl 57

# detect_prom shob_ohsc 2900

# detect_prom tisb_istr 1650

# detect_prom zoro_orzo 280

# detect_prom srnb_srnc 1400

detect_prom cjpt_cjra 15

detect_prom fst_rna1 41