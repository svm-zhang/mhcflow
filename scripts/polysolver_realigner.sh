#!/usr/bin/env bash

set -e
set -o pipefail


SRC_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
LIB_DIR="${SRC_DIR%/*}/lib"
for libfile in "$LIB_DIR"/*.sh; do
	source "$libfile"
done

function realn_usage () {
	local program
  program="${0##*/}"
	cat << EO
realn_usage: $program [options]
Options:
EO
	cat << EO | column -s\& -t
	--fqs    & Specify the file with 2 fished fastq files [Required]
	--hla_ref    & Specify the Novoalign HLA indexed file, e.g. abc_complete.nix [Required]
	--sample    & Specify the sample name [Required]
	--out    & Specify the path to the reaglined BAM file [Required]
	--mdup    & Specify the path to the reaglined BAM file (false)
	--mdup_ram    & Specify the max amount of RAM for mdup in GB (8)
	--nproc    & Specify the number of CPUs used (8)
EO
}

function realigner_main () {
	local fqs=
	local sample=
	local hla_ref=
	local out=
	local nproc=8
	local mdup=false
	local mdup_ram=8

	if [ "$#" -le 1 ];
	then
		realn_usage
		exit 1
	fi

	while [ $# -gt 0 ]; do
		case $1 in
			-h|--help)
				realn_usage
				exit 0;;
			--fqs) shift; fqs="$1";;
			--hla_ref) shift; hla_ref="$1";;
			--sample) shift; sample="$1";;
			--out) shift; out="$1";;
			--mdup) mdup=true;;
			--mdup_ram) shift; mdup_ram="$1";;
			-p|--nproc) shift; nproc="$1";;
			--) shift; break;;
			*)
				echo "Invalid option: $1" 1>&2
				realn_usage; exit 1;
				;;
		esac
		shift
	done

	check_empty_str "$sample" \
		|| die "$0" "$LINENO" \
			"realigner requires --sample to operate"
	check_empty_str "$fqs" \
		|| die "$0" "$LINENO" "realigner requires file with fastqs (--fqs)"
	check_empty_str "$hla_ref" \
		|| die "$0" "$LINENO" "realigner requires a HLA ref file (--hla_ref)"

	check_file_exists "$fqs" \
		|| die "$0" "$LINENO" "realigner failed to find provided --fqs: $fqs"
	check_file_exists "$hla_ref" \
		|| die "$0" "$LINENO" "realigner failed to find HLA ref: $hla_ref"

	fqs=$( get_abs_path "$fqs" "f" )
	hla_ref=$( get_abs_path "$hla_ref" "f" )

	info "$0" ${LINENO} "Run Polysolver HLA realigner"
	local start_time
	start_time=$(date +%s)

	local outdir="${out%/*}"
	outdir=$( make_dir "$outdir" )
	local logdir=""
	logdir=$( make_dir "$outdir/log")
	local done="${logdir}/${sample}.hla.realn.done"
	if [ -f "${done}" ]; then
		info "$0" ${LINENO} "Previous HLA realignment result exists. Skip realignment..."
		info "$0" ${LINENO} "Run Polysolver HLA realigner [DONE]"
		exit 0
	fi

	local hla_ref_nix="${hla_ref%.*}.nix"
	check_file_exists "$hla_ref_nix" \
		|| die "$0" "$LINENO" \
			"realigner failed to find indexed HLA ref: $hla_ref_nix"

	local realn_input_list=""
	local split_dir="${outdir}/splits"
	local fished_R1=""
	local fished_R2=""
	fished_R1=$( cut -f1 "$fqs" )
	fished_R2=$( cut -f2 "$fqs" )

	if [ ! -f "$fished_R1" ] || [ ! -f "$fished_R2" ]; then
		die "$0" "$LINENO" "Failed to find either fished R1 or R2 fastq file"
	fi

	local split_log="${logdir}/${sample}.seqkit.split2.log"
	realn_input_list="${split_dir}/.${sample}.realn.input.list.txt"
	seqkit split2 -j"$nproc" -p"$nproc" -O "$split_dir" \
		-1 "$fished_R1" -2 "$fished_R2" > "$split_log" 2>&1 \
		|| die "$0" "$LINENO" "Failed to run seqkit split2 on fished reads"

	if [ -f "$realn_input_list" ]; then
		rm "$realn_input_list"
	fi
	for i in $(seq 1 "$nproc"); do
		find "$split_dir" -name "*part_00$i*fastq" \
			| sort \
			| paste -sd " " - >> "$realn_input_list"
	done
	awk 'BEGIN{err=0} {if(NF > 2) {err=1}} END{exit err}' \
		"$realn_input_list" \
		|| die "$0" "$LINENO" \
			"Found more than 2 fq files in (some) line(s) in $realn_input_list"

	info "$0" ${LINENO} "Realign reads to HLA using Novoalign"
	awk '{print}' "$realn_input_list" \
		| xargs -P"$nproc" -I{} bash -c 'run_novoalign "$@"' \
			"run_novoalign" "{}" "$sample" "$hla_ref_nix" "$split_dir"

	# concat individual realigned BAM
	raw_bam="${outdir}/${sample}.hla.realn.bam"
	local bam_list="$split_dir/.${sample}.bams_to_cat.list.txt"
	find "$split_dir" -name "*.bam" > "$bam_list" \
		|| die "$0" "$LINENO" "Failed to find inidividual realigned BAM files"
	local nbams=0
	nbams=$( wc -l <"$bam_list" )
	if (( "$nbams" == 0 )); then
		die "$0" "$LINENO" "Found zero BAM files under $split_dir"
	fi
	samtools_cat "$bam_list" "$raw_bam"

	# sort concated BAM
	local so_bam="${raw_bam%.bam}.so.bam"
	if [ "$mdup" = false ]; then
		so_bam="$out"
	fi
	samtools_sort "$raw_bam" "$so_bam" "$nproc"

	rm -f "$raw_bam"
	if [ "$mdup" = true ]; then
		# if --mdup specified, mark duplicates as well
		picard_mdup "$so_bam" "$out" "$mdup_ram"
		rm -f "${so_bam}" "${so_bam%.bam}.bai" "${so_bam}.bai"
	fi

	local end_time
	local runtime
	end_time=$(date +%s)
	runtime=$(( end_time - start_time ))
	runtime=$( date -u -d @"$runtime" +'%M.%S')
	runtime_file="${outdir}/${sample}.realigner.runtime.tsv"
	printf "%s\t%s\n" "SampleID" "RealignerTime" > "$runtime_file"
	printf "%s\t%s\n" \
		"$sample" \
		"${runtime}m" \
		>> "$runtime_file"

	info "$0" ${LINENO} "Clean out intermediate files"
	rm -rf "$split_dir"

	touch "$done"
	info "$0" ${LINENO} "Realign reads to HLA using Novoalign [DONE]"
}

realigner_main "$@"

exit 0
