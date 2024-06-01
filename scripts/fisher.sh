#!/usr/bin/env bash

set -e
set -o pipefail


SRC_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
LIBCOMMON="${SRC_DIR%/*}/lib/libcommon.sh"
source "$LIBCOMMON"

function polysolver_faster () {
  local bam="$1"
  local tag_file="$2"
	local hla_bed="$3"
	local out_fished_rids="$4"

	local outdir="${out_fished_rids%/*}"
  info "$0" ${LINENO} "Fish TAG reads from alignments using faster mode"

	info "$0" ${LINENO} "Fish TAG reads from alignments on chromosome 6"
  local fished_chr6_rids="${outdir}/chr6.tag.ids"
	local chrom6=""
	chrom6=$( samtools view -H "$bam" \
		| grep "SN:6\|SN:chr6\|SN:NC000006\|SN:CM000668" \
		| cut -f2 \
		| sed 's/SN://g'
	)
	if [ -z "$chrom6" ]; then
    die "$0" "$LINENO" "Failed to extract ref name for chromosome 6 from header"
	fi
	samtools view -@"$nproc" "$bam" "$chrom6" \
		| cut -f1,10 \
		| grep -f "$tag_file" \
		| cut -f1 > "$fished_chr6_rids" \
    || die "$0" "$LINENO" "Failed to fish reads from alignments on chromosome 6"

	info "$0" ${LINENO} "Fish TAG reads from unplaced alignments"
  local fished_unmapped_rids="${outdir}/unmapped.tag.ids"
	#awk 'BEGIN{print "4\n8\n"}' | \
	#	xargs -I{} bash -c "samtools view -@$nproc -f{} -F2 $bam | cut -f1,10 | grep -f $tag_file | cut -f1 >> $fished_unmapped_rids"
	samtools view -@"$nproc" -f4 "$bam" \
		| cut -f1,10 \
		| grep -f "$tag_file" \
		| cut -f1 >> "$fished_unmapped_rids"

  local fished_hla_rids="${outdir}/fished.hla.aln.ids"
	fish_from_hla_aln "$bam" "$hla_bed" "$fished_hla_rids"
  info "$0" ${LINENO} "Get reads mapped to HLA regions in BAM [DONE]" 

	cat "$fished_hla_rids" "$fished_unmapped_rids" "$fished_chr6_rids" \
		| sort \
		| uniq > "$out_fished_rids" \
		|| die "$0" "$LINENO" "Failed to merge TAG and HLA read IDs"

	rm -f "$fished_chr6_rids" "$fished_unmapped_rids" "$fished_hla_rids"
}


function fish_from_hla_aln () {
  local bam="$1"
  local hla_bed="$2"
	local out_fished_hla_rids="$3"

  #getting chr6 region
  info "$0" ${LINENO} "Get reads mapped to HLA regions in BAM" 
  #fished_hla_rids="${outdir}/${sample}.fished.hla.aln.ids"
  samtools view -ML "$hla_bed" "$bam" \
    | cut -f1 \
    | sort \
    | uniq > "$out_fished_hla_rids" \
    || die "$0" "$LINENO" "Failed to get reads mapped to HLA regions in BAM"
}

function fisher_usage () {
	local program
  program="${0##*/}"
	cat << EO
fisher_usage: $program [options]
Options:
EO
	cat << EO | column -s\& -t
	-b or --bam    & Specify the path to the BAM file [Required]	
	-t or --tag    & Specify the TAG file, e.g. abc_v14.uniq [Required]
	--bed    & Specify the path to the HLA region defined in BED [Required]
	--sample    & Specify the sample name [Required]
	--out    & Specify the path to file hosting fished fq file list [Required]
	--nproc    & Specify the number of CPUs used [8]
EO
}

function fisherman () {
	local bam=
	local sample=
	local tag_file=
	local hla_bed=
	local out=
	local mode="polysolver_faster"
	local nproc=8

	if [ "$#" -le 1 ];
	then
		fisher_usage
		exit 1
	fi

	while [ $# -gt 0 ]; do
		case $1 in
			-h|--help)
				fisher_usage
				exit 0;;
			--sample)
				shift; sample="$1";;
			-b|--bam)
				shift; bam=$(parse_path "$1");;
			-t|--tag)
				shift; tag_file=$(parse_path "$1");;
			--bed)
				shift; hla_bed=$(parse_path "$1");;
			--out)
				shift; out="$1";;
			--mode)
				shift; mode="$1";;
			-p|--nproc)
				shift; nproc="$1";;
			--)
				shift; break;;
			*)
				echo "Invalid option: $1" 1>&2
				fisher_usage; exit 0;
				;;
		esac
		shift
	done

	if [ -z "$bam" ]; then
		error "$0" ${LINENO} "polysolver requires a BAM file to work with" 
		fisher_usage
		exit 1
	fi

	if [ -z "$sample" ]; then
		error "$0" ${LINENO} "polysolver needs a sample name" 
		fisher_usage
		exit 1
	fi

	if [ -z "$tag_file" ]; then
		error "$0" ${LINENO} "polysolver requires a TAG sequence file, such as abc_v14.uniq" 
		fisher_usage
		exit 1
	fi

	if [ -z "$hla_bed" ]; then
		error "$0" ${LINENO} "polysolver requires the HLA region file in BED" 
		fisher_usage
		exit 1
	fi

	info "$0" ${LINENO} "Run fisher for HLA-related read candidates" 
	local start_time
	start_time=$(date +%s)

	local outdir="${out%/*}"
	outdir=$( make_dir "$outdir" )
	local logdir=""
	logdir=$( make_dir "$outdir/log")

	# getting matching tag sequences
	local fish_done="${logdir}/${sample}.fish.done"
	if [ -f "$fish_done" ]; then
		info "$0" ${LINENO} "Fishing has previously done."
		exit 0
	fi

	info "$0" ${LINENO} "Start fishing" 
	local out_fished_rids="${outdir}/${sample}.fished.ids"
	if [ "$mode" = "polysolver_faster" ]; then
		polysolver_faster "$bam" "$tag_file" "$hla_bed" "$out_fished_rids"
	fi
	if [ -z "$out_fished_rids" ] || [ ! -f "$out_fished_rids" ]; then
		die "$0" "$LINENO" "Failed to find fished read id file $out_fished_rids"
	fi
	local n_fished_reads=0
	n_fished_reads=$( wc -l < "$out_fished_rids" )
	if(( "$n_fished_reads" == 0 )); then
		die "$0" "$LINENO" "Fisher failed to find any reads. Exit"
	fi

	# cash in actual reads given these ids
	info "$0" ${LINENO} "Cash in actual reads given fished IDs" 
	local fished_R1="${outdir}/${sample}.fished.R1.fastq"
	local fished_R2="${outdir}/${sample}.fished.R2.fastq"
	samtools view -@"$nproc" -bh -N "$out_fished_rids" "$bam" \
		| samtools sort -n \
		| samtools fastq -n -1 "$fished_R1" -2 "$fished_R2" -0 /dev/null -s /dev/null \
		|| die "$0" "$LINENO" "Failed to cash in actual reads given ids"

	if [ ! -f "$fished_R1" ] || [ ! -f "$fished_R2" ]; then
		die "$0" "$LINENO" "Failed to find either $fished_R1 or $fished_R2 fastq file"
	fi
	echo -e "$fished_R1\n$fished_R2" > "$out"

	local end_time
	local runtime
	end_time=$(date +%s)
	runtime=$(( end_time - start_time ))
	runtime=$( date -u -d @"$runtime" +'%M.%S')
	# collect the number of fished reads
	fisher_stat_file="${outdir}/${sample}.fisher.stat.tsv"
	printf "%s\t%s\t%s\n" "SampleID" "NumFished" "FishTime" > "$fisher_stat_file"
	printf "%s\t%s\t%s\n" \
		"$sample" \
		"$n_fished_reads" \
		"${runtime}m" >> "$fisher_stat_file"

	touch "$fish_done"
	info "$0" ${LINENO} "Fish done, time to go home" 
}

fisherman "$@"