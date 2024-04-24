#!/usr/bin/env bash

# Given a raw (unsorted, non-deduped) BAM, generate a sorted and
# duplciates-marked BAM

set -e

SRC_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &>/dev/null && pwd)
#COMMON_FUNC_LIB_DIR="${SRC_DIR%/*}/lib"
COMMON_FUNC_LIB="${SRC_DIR}/hla_bash_util_funcs"
source "${COMMON_FUNC_LIB}"

function usage () {
	local program
	program=$(basename "$0")
	cat << EO
Usage: $program [options]
Options:
EO
	cat << EO | column -s\& -t
	-b or --bam    & Specify path to the BAM file (Required)
	-o or --out    & Specify path to the output BAM file (Required)
	--mdup    & Specify to run picard markduplicates [false]
	--rmdup    & Specify to run samtools rmdup [false]
	--nproc    & Specify the number of CPUs used [8]
EO
}

bam=
out=
mdup=false
rmdup=false
nproc=8

if [ "$#" -le 1 ];
then
	usage
	exit 1
fi

while [ $# -gt 0 ]; do
	case $1 in
		-h|--help)
			usage
			exit 0;;
		--bam)
			shift; bam=$(parse_path "$1");;
		--out)
			shift; out="$1";;
		--mdup)
			mdup=true;;
		--rmdup)
			rmdup=true;;
		-p|--nproc)
			shift; nproc="$1";;
		--)
			shift; break;;
		*)
			echo "Invalid option: $1" 1>&2
			usage; exit 1;
			;;
	esac
	shift
done

if [ "${mdup}" = true ] && [ "${rmdup}" = true ]; then
	error "$0" ${LINENO} "--rmdup and --mdup should not be specified at the same time"
	exit 1
fi

if [ -f "${out}" ];
then
	info "$0" ${LINENO} "Found output file from previous run: ${out}"
	exit 0
fi

check_file_exists "${bam}"

info "$0" ${LINENO} "Sort BAM file: ${bam}"
so_bam="${bam%.bam}.so.bam"
if [ "${mdup}" = false ] && [ "${rmdup}" = false ];
then
	so_bam="${out}"
fi
if [ ! -f "${so_bam}" ]; then
	so_tmp="${bam%.bam}.tmp"
	so_bam_bai="$so_bam.bai"
	# reference: https://github.com/samtools/samtools/issues/1196
	cmd="samtools sort -T$so_tmp -@$nproc --write-index -o $so_bam##idx##$so_bam_bai $bam 2>/dev/null"
	if ! eval "${cmd}";
	then
		error "$0" ${LINENO} "Failed to sort BAM file: ${bam}"
		exit 1 
	fi
	info "$0" ${LINENO} "Sorting realignment [DONE]"
else
	info "$0" ${LINENO} "Found previous sorted BAM ${so_bam}. Skip"
fi
check_file_exists "$so_bam"

if [ "${rmdup}" = true ]; then
	# run samtools rmdup
	info "$0" ${LINENO} "Removing PCR duplicates using samtools rmdup"
	cmd="samtools rmdup ${so_bam} ${out} >/dev/null 2>&1 && samtools index ${out}"
	if ! eval "$cmd";
	then
		error "$0" ${LINENO} "Failed to run command: ${cmd}"
		exit 1
	fi
	info "$0" ${LINENO} "Removing PCR duplicates using samtools rmdup [DONE]"
fi

if [ "${mdup}" = true ];
then
	# run picard markduplicates
	info "$0" ${LINENO} "Removing PCR duplicates using picard markduplicates"
	metric="${so_bam%.bam}.mdup.metrics.txt"
	cmd="picard -XX:ParallelGCThreads=1 -Xmx8g MarkDuplicates --INPUT ${so_bam} --OUTPUT ${out} --METRICS_FILE ${metric} --CREATE_INDEX true"
	if ! eval "${cmd}";
	then
		error "$0" ${LINENO} "Failed to run command: ${cmd}"
		exit 1
	fi
	info "$0" ${LINENO} "Removing PCR duplicates using picard markduplicates [DONE]"
fi

exit 0
