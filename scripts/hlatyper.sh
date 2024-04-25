#!/usr/bin/env bash

set -e


function usage () {
	local program
	program=$(basename "$0")
	cat << EO
Usage: $program [options]
Options:
EO
	cat << EO | column -s\& -t
	--sample    & Specify the sample name [Required]
	--bam    & Specify realigned HLA BAM [Required]
	--freq    & Specify the HLA allele population frequency file [Required]
	-r or --race    & Specify the race (Caucasian, Black, Asian, Unknown)
	-o or --outdir    & Specify the path to the output directory [Required]
	--excludeFreq    & Specify to not include population frequency as prior
	--nproc    & Specify the number of CPUs used [8]
EO
}

# common function script
SRC_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
#COMMON_FUNC_LIB_DIR="${SRC_DIR%/*}/lib"
COMMON_FUNC_LIB="${SRC_DIR}/hla_bash_util_funcs"
source "${COMMON_FUNC_LIB}"

sample=
bam=
freq_file=
race=
outdir=
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
		--sample)
			shift; sample="$1";;
		--bam)
			shift; bam=$(parse_path "$1");;
		--freq)
			shift; freq_file=$(parse_path "$1");;
		-r|--race)
			shift; race="$1";;
		-o|--outdir)
			shift; outdir="$1";;
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

check_file_exists "${bam}"
if [ -z "$sample" ];
then
	error "$0" ${LINENO} "polysolver needs a sample name" 
	usage
	exit 1
fi

if [ -z "$freq_file" ];
then
	error "$0" ${LINENO} "polysolver requires the HLA frequency file" 
	usage
	exit 1
else
  check_file_exists "${freq_file}"
fi

if [ -z "$outdir" ];
then
	error "$0" ${LINENO} "polysolver requires an output directory" 
	usage
	exit 1
fi

if [ -z "$race" ];
then
	race="Unknown"
fi

case "$race" in
	Caucasian|Black|Asian|Unknown) true;;
	*) error "$0" ${LINENO} "Race can be either Caucasian, or Black, or Asian, or Unknown";
		exit 1;;
esac

outdir=$(make_dir "$outdir")

a1_dir=$(make_dir "${outdir}/a1")
a1_log_dir=$(make_dir "${a1_dir}/log")
a1_loglik_all="${outdir}/${sample}.a1.loglik.tsv"

info "$0" "$LINENO" "Get a list of HLA alleles to calculate log-likelihood score"
hla_ids_file="${outdir}/hla_ids.txt"
#grep "^>" "$hla_ref" | sed 's/^>//g' > "$hla_ids_file"
cmd="samtools view -H $bam | grep '^@SQ' | cut -f2 | sed 's/^SN://g' > $hla_ids_file"
run_cmd "$cmd" "$LINENO" "Failed to get HLA alleles"
info "$0" "$LINENO" "Get a list of HLA alleles to calculate log-likelihood score [DONE]"

if [ ! -f "$a1_loglik_all" ]; then
	info "$0" ${LINENO} "Calculating likelihood score for the first HLA allele"
	< "$hla_ids_file" tee | \
		xargs -P"$nproc" -I{} bash -c "samtools view -bh $bam {} | \
										samtools sort -n | samtools view | \
										first_allele_calculations $race {} $freq_file $a1_dir \
										>$a1_log_dir/$sample.{}.log"

	n_a1_lik_files=$(find "$a1_dir" -name "*.lik1" | wc -l)
	if [ "$n_a1_lik_files" -eq 0 ];
  then
		error "$0" ${LINENO} "Found no likelihood score for any A1 alleles"
		exit 1
	fi

	find "$a1_dir" -name "*.lik1" | \
		awk '{n=split($1, a, "/"); gsub(/\.lik1/, "", a[n]); print a[n];}' | \
		xargs -I{} bash -c "tail -1 ${a1_dir}/{}.lik1 | cut -f2 | paste <(echo {}) - >> $a1_loglik_all"

	info "$0" ${LINENO} "Calculating likelihood score for the first HLA allele [DONE]" 
else
	info "$0" ${LINENO} "Previous calculation for the first HLA allele exists. Skipping..." 
fi

a2_dir=$(make_dir "${outdir}/a2")
a2_log_dir=$(make_dir "${a2_dir}/log")
a2_log_file="${a2_log_dir}/${sample}.a2.log"
a2_loglik_all="${outdir}/${sample}.a2.loglik.tsv"
if [ ! -f "$a2_loglik_all" ]; then

	info "$0" ${LINENO} "Calculating likelihood score for the second HLA allele" 
	# check if a1_dir exists with contents
	# check if a1 loglik file exists
	content_filt_a1=$(head "${a1_loglik_all}")	
	if [ -z "${content_filt_a1}" ];
	then
		error "$0" ${LINENO} "All A1 alelle calculation return score of zero"
	fi
	filt_a1_loglik_all="${outdir}/${sample}.a1.loglik.filt.tsv"
	awk '($2 > 0)' "$a1_loglik_all" > "$filt_a1_loglik_all"
	second_allele_calculations "$race" "$filt_a1_loglik_all" "$hla_ids_file" \
		"$freq_file" "$a1_dir" "$a2_dir" > "$a2_log_file" 

	n_a2_lik_files=$(find "$a2_dir"/ -name "*.lik2" | wc -l)
	if [ "$n_a2_lik_files" -eq 0 ]; then
		error "$0" ${LINENO} "Found no likelihood score for any A2 alleles"
		exit 1
	fi

	find "$a2_dir" -name "*.lik2" | \
		awk '{n=split($1, a, "/"); gsub(/\.lik2/, "", a[n]); print a[n];}' | \
		xargs -I{} bash -c "tail -1 ${a2_dir}/{}.lik2 | cut -f2 | paste <(echo {}) - >> $a2_loglik_all"
	info "$0" ${LINENO} "Calculating likelihood score for the second HLA allele [DONE]" 
else
	info "$0" ${LINENO} "Previous calculation for the second HLA allele exists. Skipping..." 
fi

hla_typing_res="${outdir}/${sample}.hla_typing.tsv"
if [ ! -f "$hla_typing_res" ]; then
	winner1_a=$(grep "^hla_a" "$a1_loglik_all" | sort -k2,2nr | awk '(NR==1) {print $1}')
	winner1_b=$(grep "^hla_b" "$a1_loglik_all" | sort -k2,2nr | awk '(NR==1) {print $1}')
	winner1_c=$(grep "^hla_c" "$a1_loglik_all" | sort -k2,2nr | awk '(NR==1) {print $1}')
	winner2_a=$(grep "^hla_a" "$a2_loglik_all" | sort -k2,2nr | awk '(NR==1) {print $1}')
	winner2_b=$(grep "^hla_b" "$a2_loglik_all" | sort -k2,2nr | awk '(NR==1) {print $1}')
	winner2_c=$(grep "^hla_c" "$a2_loglik_all" | sort -k2,2nr | awk '(NR==1) {print $1}')

	info "$0" ${LINENO} "winners1 $winner1_a $winner1_b $winner1_c"
	info "$0" ${LINENO} "winners2 $winner2_a $winner2_b $winner2_c"

	printf "SampleID\tHLA-A_1\tHLA-A_2\tHLA-B_1\tHLA-B_2\tHLA-C_1\tHLA-C_2\n" > "$hla_typing_res"
	printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "$sample" "$winner1_a" "$winner2_a" "$winner1_b" "$winner2_b" "$winner1_c" "$winner2_c" >> "$hla_typing_res"
else
	info "$0" ${LINENO} "Previous HLA typing result exists. Skipping..." 
fi

info "$0" ${LINENO} "JSPolysolver class-I HLA typing [DONE]" 

exit 0