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
bam=
typing_res=
hla_ref=
out_bam=
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
  --bam)
    shift
    bam=$(parse_path "$1")
    ;;
  --typeres)
    shift
    typing_res=$(parse_path "$1")
    ;;
	--hla_ref)
		shift
		hla_ref=$(parse_path "$1")
		;;
	-o | --out)
		shift
		out_bam="$1"
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

# 1.1 get sample-level HLA reference based on hlatyping result

# 1.2 make novoindex off the reference

# 1.3 bam2fq (if I decide to do this way)
# or we could have the file ready from razer_realigner.sh

# 1.4 split fastq files using run_seqkit_split2 function

# 1.5 run_novoalign_batch function

# 1.6 concat bam file
