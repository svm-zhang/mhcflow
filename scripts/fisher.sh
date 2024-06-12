#!/usr/bin/env bash

set -e
set -o pipefail


SRC_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
LIBCOMMON="${SRC_DIR%/*}/lib/libcommon.sh"
source "$LIBCOMMON"

function fish () {
  local bam="$1"
  local tag_file="$2"
	local hla_bed="$3"
	local out_fished_rids="$4"
	local mode="$5"

	local outdir="${out_fished_rids%/*}"
  local fished_tag_rids="${outdir}/fished.tag.ids"
	case "$mode" in
		faster) fish_tag_faster "$bam" "$tag_file" "$fished_tag_rids";;
		fast) fish_tag_fast "$bam" "$tag_file" "$fished_tag_rids";;
		*) die "${FUNCNAME[0]}" "$LINENO" "Unrecognized fishing mode $mode"
	esac

  local fished_hla_rids="${outdir}/fished.hla.aln.ids"
	fish_from_hla_aln "$bam" "$hla_bed" "$fished_hla_rids"

	cat "$fished_hla_rids" "$fished_tag_rids" \
		| sort \
		| uniq > "$out_fished_rids" \
		|| die "${FUNCNAME[0]}" "$LINENO" "Failed to merge TAG and HLA read IDs"

	rm -f "$fished_hla_rids" "$fished_tag_rids"
}

function fish_tag_faster () {
  local bam="$1"
  local tag_file="$2"
	local out="$3"

  info "${FUNCNAME[0]}" ${LINENO} \
		"Fish TAG reads from alignments using faster mode"
	info "${FUNCNAME[0]}" ${LINENO} \
		"Fish TAG reads from alignments on chromosome 6"
	local chrom6=""
	chrom6=$( samtools view -H "$bam" \
		| { grep "SN:6\|SN:chr6\|SN:NC000006\|SN:CM000668" || test "$?" = 1; } \
		| cut -f2 \
		| sed 's/SN://g' \
    || die "${FUNCNAME[0]}" "$LINENO" "Failed to grep chrom6 ref in BAM header"
	)
	if [ -z "$chrom6" ]; then
    die "${FUNCNAME[0]}" "$LINENO" \
			"Failed to extract ref name for chromosome 6 from header"
	fi
	samtools view -@"$nproc" "$bam" "$chrom6" \
		| cut -f1,10 \
		| { grep -F -f "$tag_file" || true; } \
		| cut -f1 > "$out" \
    || die "${FUNCNAME[0]}" "$LINENO" \
			"Failed to fish reads from alignments on chromosome 6"

	info "${FUNCNAME[0]}" "$LINENO" "Fish TAG reads from unplaced alignments"
	samtools view -@"$nproc" -f4 "$bam" \
		| cut -f1,10 \
		| { grep -F -f "$tag_file" || true; } \
		| cut -f1 >> "$out"
  info "${FUNCNAME[0]}" "$LINENO" \
		"Fish TAG reads from alignments using faster mode [DONE]"
}

function fish_tag_fast () {
	local bam="$1"
  local tag_file="$2"
	local out="$3"

  info "${FUNCNAME[0]}" "$LINENO" \
		"Fish TAG reads from alignments using fast mode"
	local bam_size=0
	bam_size=$( du -bs "$bam" | awk '{print $1/2^30}' )
	if [ -z "$bam_size" ]; then
		die "${FUNCNAME[0]}" "$LINENO" "Failed to obtaint input BAM file size"
	fi
	info "${FUNCNAME[0]}" "$LINENO" \
		"Obtain input BAM file size: ${bam_size}G" 
	info "${FUNCNAME[0]}" "$LINENO" \
		"Split BAM into parts for faster fishing reads"
	local outdir="${out%/*}"
	local split_bam_dir="${outdir}/bams"
	split_bam_dir=$( make_dir "$split_bam_dir" )
	bam_size_int=$( printf "%.0f" "$bam_size" )
	split_size=$(( (bam_size_int+1+nproc-1)/nproc ))
	info "${FUNCNAME[0]}" "$LINENO" "Split BAM file into ${split_size}G each"
	# split bam files
	samtools view "$bam" \
		| cut -f1,10 \
		| split -C "${split_size}G" -d -a 2 \
			--additional-suffix ".txt" - "$split_bam_dir/" \
		|| die "${FUNCNAME[0]}" "$LINENO" \
			"Failed to split input bam to parts for fishing tags"

	info "${FUNCNAME[0]}" "$LINENO" "Fish reads from each individual split file"
	find "$split_bam_dir" -name "*.txt" \
		| sed 's/\.txt//g' \
		| xargs -P"$nproc" -I{} bash -c "grep -F -f $tag_file {}.txt | cut -f1 >> $out"

	if [ -d "$split_bam_dir" ]; then
		rm -rf "$split_bam_dir"
	fi
  info "${FUNCNAME[0]}" "$LINENO" \
		"Fish TAG reads from alignments using fast mode [DONE]"
}

function fish_from_hla_aln () {
  local bam="$1"
  local hla_bed="$2"
	local out_fished_hla_rids="$3"

  #getting chr6 region
  info "${FUNCNAME[0]}" "$LINENO" "Get reads mapped to HLA regions in BAM" 
  samtools view -ML "$hla_bed" "$bam" \
    | cut -f1 \
    | sort \
    | uniq > "$out_fished_hla_rids" \
    || die "${FUNCNAME[0]}" "$LINENO" \
			"Failed to get reads mapped to HLA regions in BAM"
  info "${FUNCNAME[0]}" "$LINENO" \
		"Get reads mapped to HLA regions in BAM [DONE]" 
}

function fisher_usage () {
	local program
  program="${0##*/}"
	cat << EO
fisher_usage: $program [options]
Options:
EO
	cat << EO | column -s\& -t
	--bam    & Specify the path to the BAM file [Required]
	--tag    & Specify the TAG file, e.g. abc_v14.uniq [Required]
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
	local mode="faster"
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
			--sample) shift; sample="$1";;
			--bam) shift; bam="$1";;
			--tag) shift; tag_file="$1";;
			--bed) shift; hla_bed="$1";;
			--out) shift; out="$1";;
			--mode) shift; mode="$1";;
			-p|--nproc) shift; nproc="$1";;
			--) shift; break;;
			*)
				echo "Invalid option: $1" 1>&2
				fisher_usage; exit 0;
				;;
		esac
		shift
	done

	check_empty_str "$sample" \
		|| die "${FUNCNAME[0]}" "$LINENO" \
			"fisher requires --sample to operate"
	check_empty_str "$bam" \
		|| die "${FUNCNAME[0]}" "$LINENO" "fisher requires a BAM file(--bam)"
	check_empty_str "$tag_file" \
		|| die "${FUNCNAME[0]}" "$LINENO" \
			"fisher requires HLA Kmer tag file (--tag)"
	check_empty_str "$hla_bed" \
		|| die "${FUNCNAME[0]}" "$LINENO" \
			"fisher requires a HLA BED file (--hla_bed)"

	check_file_exists "$bam" \
		|| die "${FUNCNAME[0]}" "$LINENO" \
			"fisher failed to find provided BAM file: $bam"
	check_file_exists "$tag_file" \
		|| die "${FUNCNAME[0]}" "$LINENO" \
			"fisher failed to find provided KMER taga file: $tag_file"
	check_file_exists "$hla_bed" \
		|| die "${FUNCNAME[0]}" "$LINENO" \
			"fisher failed to find provided BED file: $hla_bed"

	bam=$( get_abs_path "$bam" "f" )
	hla_bed=$( get_abs_path "$hla_bed" "f" )
	tag_file=$( get_abs_path "$tag_file" "f" )

	#if [ -z "$bam" ]; then
	#	error "$0" ${LINENO} "polysolver requires a BAM file to work with" 
	#	fisher_usage
	#	exit 1
	#fi

	#if [ -z "$sample" ]; then
	#	error "$0" ${LINENO} "polysolver needs a sample name" 
	#	fisher_usage
	#	exit 1
	#fi

	#if [ -z "$tag_file" ]; then
	#	error "$0" ${LINENO} "polysolver requires a TAG sequence file, such as abc_v14.uniq" 
	#	fisher_usage
	#	exit 1
	#fi

	#if [ -z "$hla_bed" ]; then
	#	error "$0" ${LINENO} "polysolver requires the HLA region file in BED" 
	#	fisher_usage
	#	exit 1
	#fi

	info "${FUNCNAME[0]}" ${LINENO} "Run fisher for HLA-related read candidates" 
	local start_time
	start_time=$(date +%s)

	local outdir="${out%/*}"
	outdir=$( make_dir "$outdir" )
	local logdir=""
	logdir=$( make_dir "$outdir/log")

	# getting matching tag sequences
	local fish_done="${logdir}/${sample}.fish.done"
	if [ -f "$fish_done" ]; then
		info "${FUNCNAME[0]}" ${LINENO} "Fishing has previously done."
		exit 0
	fi

	info "${FUNCNAME[0]}" ${LINENO} "Start fishing" 
	local out_fished_rids="${outdir}/${sample}.fished.ids"
	fish "$bam" "$tag_file" "$hla_bed" "$out_fished_rids" "$mode"
	if [ -z "$out_fished_rids" ] || [ ! -f "$out_fished_rids" ]; then
		die "${FUNCNAME[0]}" "$LINENO" \
			"Failed to find fished read id file $out_fished_rids"
	fi
	local n_fished_reads=0
	n_fished_reads=$( wc -l < "$out_fished_rids" )
	if(( "$n_fished_reads" == 0 )); then
		die "${FUNCNAME[0]}" "$LINENO" "Fisher failed to find any reads. Exit"
	fi

	# cash in actual reads given these ids
	info "${FUNCNAME[0]}" ${LINENO} "Cash in actual reads given fished IDs" 
	local fished_R1="${outdir}/${sample}.fished.R1.fastq"
	local fished_R2="${outdir}/${sample}.fished.R2.fastq"
	samtools view -@"$nproc" -bh -N "$out_fished_rids" "$bam" \
		| samtools sort -n \
		| samtools fastq -n -1 "$fished_R1" -2 "$fished_R2" -0 /dev/null -s /dev/null \
		|| die "${FUNCNAME[0]}" "$LINENO" "Failed to cash in actual reads given ids"

	if [ ! -f "$fished_R1" ] || [ ! -f "$fished_R2" ]; then
		die "${FUNCNAME[0]}" "$LINENO" \
			"Failed to find either $fished_R1 or $fished_R2 fastq file"
	fi
	echo -e "$fished_R1\t$fished_R2" > "$out"

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
	info "${FUNCNAME[0]}" ${LINENO} "Fish done, time to go home" 
}

fisherman "$@"