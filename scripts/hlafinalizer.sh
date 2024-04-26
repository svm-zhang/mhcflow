#!/usr/bin/env bash

# script to finalize HLA typing by doing thw following two things
# 1. generate a subject/sample-level HLA reference based on typing result
# 2. realign all reads used for typing again the subject/sample-level reference

set -e

SRC_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &>/dev/null && pwd)
#COMMON_FUNC_LIB_DIR="${SRC_DIR%/*}/lib"
COMMON_FUNC_LIB="${SRC_DIR}/libcommon.sh"
LIBREALIGN="${SRC_DIR}/librealign.sh"
source "$COMMON_FUNC_LIB"
source "$LIBREALIGN"

# FIXME: need too update the help message
function usage () {
	local program
	program=$(basename "$0")
	cat <<EO
Usage: $program [options]
Options:
EO
	cat <<EO | column -s\& -t
	--r1    & Specify the path to R1 fastq (Required)	
	--r2    & Specify the path to R2 fastq (Required)	
	--hla_ref    & Specify the HLA reference sequences in Fasta (Required)
	--sample    & Specify the sample name (Required)
	-o or --out    & Specify the path to the output BAM file (Required)
	--nproc    & Specify the number of CPUs used [8]
EO
}

sample=
realn_dir=
typing_res=
hla_ref=
outdir=
nproc=8

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
		shift
		sample="$1"
		;;
  --realn_dir)
    shift
    realn_dir="$1"
    ;;
  --typeres)
    shift
    typing_res=$(parse_path "$1")
    ;;
	--hla_ref)
		shift
		hla_ref=$(parse_path "$1")
		;;
	--outdir)
		shift
		outdir="$1"
		;;
	-j | --nproc)
		shift
		nproc="$1"
		;;
	--)
		shift
		break
		;;
	*)
		echo "Invalid option: $1" 1>&2
		usage
		exit 1
		;;
	esac
	shift
done

# FIXME: make a final folder
check_dir_exists "$realn_dir"
outdir=$( make_dir "$outdir" )
logdir="$outdir/log"
logdir=$( make_dir "$logdir" )

# 1.1 get sample-level HLA reference based on hlatyping result
info "$0" "$LINENO" "Get sample-level HLA reference sequence"
hlatypelist_file="$outdir/.$sample.hlatypelist.txt"
tail -n +2 "$typing_res" | cut -d $'\t' -f 2- | sed 's/\t/\n/g' > "$hlatypelist_file" 
nexpect=$( wc -l <"$hlatypelist_file" )

sample_hla_ref="$outdir/$sample.hla.fasta"
cmd=("seqkit" "grep" "-f" "$hlatypelist_file" "-o" "$sample_hla_ref" "$hla_ref" )
if ! "${cmd[@]}" >/dev/null 2>&1;
then
	error "$0" "$LINENO" "Failed to get sample-level HLA reference sequence in a Fasta file"
	exit 1
fi

# check if Fasta is empty
if [ ! -s "$sample_hla_ref" ];
then
	error "$0" "$LINENO" "Failed to extract any sequences in ${sample_hla_ref##*/}"
	exit 1
fi

# check if sequences from 6 alleles were generated
nseq=$( grep -c ">" "$sample_hla_ref" )
if [ "$nseq" -ne "$nexpect" ];
then
	error "$0" "$LINENO" "Expect to extract $nexpect alleles, but got $nseq"
	exit 1
fi
info "$0" "$LINENO" "Get sample-level HLA reference sequence [DONE]"

# 1.2 make novoindex off the reference
info "$0" "$LINENO" "Index HLA reference using novoindex"
sample_hla_ref_nix="${sample_hla_ref%.fasta}.nix"
cmd=("novoindex" "$sample_hla_ref_nix" "$sample_hla_ref")
if ! "${cmd[@]}" >/dev/null 2>&1;
then
	error "$0" "$LINENO" "Failed to index the sample-level HLA reference: ${sample_hla_ref##*/}"
	exit 1
fi
info "$0" "$LINENO" "Index HLA reference using novoindex [DONE]"

# 1.3 run_novoalign_batch function
info "$0" "$LINENO" "Realign paired reads to sample-level HLA reference using Novoalign"
pair_dir="$realn_dir/pair"
novo_dir="$outdir/novoalign"
novo_dir=$( make_dir "$novo_dir" )
check_dir_exists "$pair_dir"
fq_search_regex=".part_*.fastq"
run_novoalign_batch "$pair_dir" "$novo_dir" "$fq_search_regex" "$sample_hla_ref_nix" "$nproc"
info "$0" "$LINENO" "Realign paired reads to sample-level HLA reference using Novoalign [DONE]"

# 1.4 concat bam file
info "$0" ${LINENO} "Concatenate individually novoaligned BAM files" 
bam_search_regex="*.bam"
realn_bam="$outdir/$sample.hla.realn.bam"
run_samtools_cat "$novo_dir" "$bam_search_regex" "$realn_bam"
info "$0" ${LINENO} "Concatenate individually novoaligned BAM files [DONE]" 

# 1.4 concat bam file
info "$0" ${LINENO} "Post-process realigned BAM file" 
realn_so_bam="${realn_bam%.bam}.so.bam"
run_samtools_sort "$realn_bam" "$realn_so_bam" "$nproc"
info "$0" ${LINENO} "Post-process realigned BAM file [DONE]" 
info "$0" ${LINENO} "Run Polysolver HLA realigner [DONE]" 

exit 0
