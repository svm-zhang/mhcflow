#!/usr/bin/env bash

# script to finalize HLA typing by doing thw following two things
# 1. generate a subject/sample-level HLA reference based on typing result
# 2. realign all reads used for typing again the subject/sample-level reference

set -e

SRC_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &>/dev/null && pwd)
LIBCOMMON="${SRC_DIR%/*}/lib/libcommon.sh"
LIBREALIGN="${SRC_DIR%/*}/lib/librealign.sh"
source "$LIBCOMMON"
source "$LIBREALIGN"

function usage () {
	local program
	program=$(basename "$0")
	cat <<EO
Usage: $program [options]
Options:
EO
	cat <<EO | column -s\& -t
	--sample    & Specify the sample name (Required)
	--realn_dir    & Specify the path to the realigner directory (Required)
	--typeres    & Specify the path to the HLA typing result (Required)
	--hla_ref    & Specify the HLA reference sequences in Fasta (Required)
	--realigner    & Specify realigner [polysolver] (polysolver, hlareforged) 
	--outdir    & Specify the path to the output directory (Required)
	--nproc    & Specify the number of CPUs used [8]
	--mdup_ram    & Specify the max amount of RAM for mdup in GB [8]
EO
}

sample=
realn_dir=
typing_res=
hla_ref=
realigner="polysolver"
outdir=
nproc=8
mdup_ram=8

if [ "$#" -le 1 ]; then
	usage
	exit 1
fi

while [ $# -gt 0 ]; do
	case $1 in
	-h | --help)
		usage
		exit 0
		;;
	--sample)
		shift; sample="$1";;
  --realn_dir)
    shift; realn_dir="$1";;
  --typeres)
    shift; typing_res=$(parse_path "$1");;
	--hla_ref)
		shift; hla_ref=$(parse_path "$1");;
  --realigner)
    shift; realigner="$1";;
	--outdir)
		shift; outdir="$1";;
	-j | --nproc)
		shift; nproc="$1";;
	--mdup_ram)
		shift; mdup_ram="$1";;
	--)
		shift; break;;
	*)
		echo "Invalid option: $1" 1>&2
		usage
		exit 1
		;;
	esac
	shift
done

case "$realigner" in
	polysolver|hlareforged) true;;
	*) die "$0" ${LINENO} "realigner can only be either polysolver or hlareforged"
esac

info "$0" "$LINENO" "Finalize HLA realigner and typer results"
check_dir_exists "$realn_dir"
outdir=$( make_dir "$outdir" )
logdir="$outdir/log"
ready_bam="$outdir/${sample}.hla.realn.ready.bam"
donefile="$logdir/${sample}.finalizer.done"
logdir=$( make_dir "$logdir" )

if [ -f "$ready_bam" ] && [ -f "$donefile" ]; then
	info "$0" "$LINENO" "Found previous finalizer results. Skip"
	info "$0" "$LINENO" "Finalize HLA realigner and typer results [DONE]"
	exit 0
fi

# 1.1 get sample-level HLA reference based on hlatyping result
info "$0" "$LINENO" "Get sample-level HLA reference sequence"
hlatypelist_file="$outdir/.$sample.hlatypelist.txt"
#tail -n +3 "$typing_res" | cut -d $'\t' -f 2- | sed 's/\t/\n/g' > "$hlatypelist_file"
awk '(NR > 1) {print $1}' "$typing_res" > "$hlatypelist_file"
nexpect=$( sort "$hlatypelist_file" | uniq | wc -l )
sample_hla_ref="$outdir/$sample.hla.fasta"
if [ ! -f "$sample_hla_ref" ] || [ ! -s "$sample_hla_ref" ]; then
	seqkit_grep_log="$logdir/$sample.seqkit_grep.log"
	seqkit grep -f "$hlatypelist_file" -o "$sample_hla_ref" "$hla_ref" >"$seqkit_grep_log" 2>&1 \
		|| die "$0" "$LINENO" "Failed to grep sequences of typed HLA alleles"
else
	info "$0" "$LINENO" "Found sample-level HLA sequence file: $sample_hla_ref. Continue"
fi

# check if Fasta is empty
if [ ! -s "$sample_hla_ref" ]; then
	die "$0" "$LINENO" "Failed to extract any sequences in ${sample_hla_ref##*/}"
fi

# check if sequences from 6 alleles were generated
nseq=$( grep -c ">" "$sample_hla_ref" )
if [ "$nseq" -ne "$nexpect" ]; then
	rm -f "$sample_hla_ref"
	die "$0" "$LINENO" "Expect to extract $nexpect alleles, but got $nseq"
fi
info "$0" "$LINENO" "Get sample-level HLA reference sequence [DONE]"

# 1.2 make novoindex off the reference
info "$0" "$LINENO" "Index HLA reference using novoindex"
sample_hla_ref_nix="${sample_hla_ref%.fasta}.nix"
novoidx_log="$logdir/$sample.novoindex.log"
if [ ! -f "$sample_hla_ref_nix" ]; then
	novoindex "$sample_hla_ref_nix" "$sample_hla_ref" >"$novoidx_log" 2>&1 \
		|| die "$0" "$LINENO" "Failed to index the sample-level HLA sequences"
else
	info "$0" "$LINENO" "Found NIX built previously: $sample_hla_ref_nix. Continue"
fi
info "$0" "$LINENO" "Index HLA reference using novoindex [DONE]"

# 1.3 run_novoalign_batch function
info "$0" "$LINENO" "Realign paired reads to sample-level HLA reference using Novoalign"
novodonefile="$logdir/novoalign.done"
realn_bam="$outdir/$sample.hla.realn.bam"
if [ ! -f "$novodonefile" ] || [ ! -f "$realn_bam" ]; then
	fq_dir=""
	fq_search_regex=""
	# this is not pretty. need to think
	if [ "$realigner" = "polysolver" ]; then
		fq_dir="$realn_dir/splits"
		fq_search_regex="*.fastq"
	else
		fq_dir="$realn_dir/pair"
		fq_search_regex=".part_*.fastq"
	fi
	check_dir_exists "$fq_dir"
	novo_dir="$outdir/novoalign"
	novo_dir=$( make_dir "$novo_dir" )
	info "$0" "$LINENO" "Check if intermediate paired reads from realigner exist"
	run_novoalign_batch "$fq_dir" "$novo_dir" "$fq_search_regex" "$sample" "$sample_hla_ref_nix" "$nproc"

	# 1.4 concat bam file
	info "$0" ${LINENO} "Concatenate individually novoaligned BAM files"
	bam_search_regex="*.bam"
	run_samtools_cat "$novo_dir" "$bam_search_regex" "$realn_bam"
	info "$0" ${LINENO} "Concatenate individually novoaligned BAM files [DONE]"

	touch "$novodonefile"
	rm -rf "$novo_dir" "$fq_dir"
	info "$0" "$LINENO" "Realign paired reads to sample-level HLA reference using Novoalign [DONE]"
else
	info "$0" "$LINENO" "Found aligned BAM file: $realn_bam"
fi

# 1.5 sort bam file
info "$0" ${LINENO} "Post-process realigned BAM file" 
sort_donefile="$logdir/samtools_sort.done"
realn_so_bam="${realn_bam%.bam}.so.bam"
if [ ! -f "$sort_donefile" ] || [ ! -f "$realn_so_bam" ]; then
	run_samtools_sort "$realn_bam" "$realn_so_bam" "$nproc"
else
	info "$0" "$LINENO" "Previous sorted BAM file found. Skip"
fi
info "$0" ${LINENO} "Post-process realigned BAM file [DONE]" 

# 1.5 mdup bam file
info "$0" ${LINENO} "Mark PCR duplicates using picard markduplicates"
mdup_donefile="$logdir/mdup.done"
if [ ! -f "$mdup_donefile" ] || [ ! -f "$ready_bam" ]; then
	run_picard_mdup "$realn_so_bam" "$ready_bam" "$mdup_ram"
else
	info "$0" "$LINENO" "Previous dup-marked BAM file found. Skip"
fi
info "$0" ${LINENO} "Mark PCR duplicates using picard markduplicates [DONE]"

rm -f "$realn_bam" "$realn_so_bam" "${realn_so_bam%.bam}.bam.bai" 

info "$0" "$LINENO" "Finalize HLA realigner and typer results [DONE]"

touch "$donefile"

exit 0
