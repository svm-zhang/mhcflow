#!/usr/bin/env bash

set -e

SRC_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &>/dev/null && pwd)
COMMON_FUNC_LIB="${SRC_DIR}/libcommon.sh"
LIBREALIGN="${SRC_DIR}/librealign.sh"
source "$COMMON_FUNC_LIB"
source "$LIBREALIGN"

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
	--outdir    & Specify the path to the output directory (Required)
	--nproc    & Specify the number of CPUs used [8]
	--no_gunzip    & Specify skip decompressing fastq.gz prior to Novoalign
	--overwrite    & Specify to overwrite all previous results and start new
EO
}

r1=
r2=
hla_ref=
hla_ref_nix=
sample=
outdir=
nproc=8
nproc_per_job=2
overwrite=false
no_gunzip=false

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
	--r1)
		shift
		r1=$(parse_path "$1")
		;;
	--r2)
		shift
		r2=$(parse_path "$1")
		;;
	--hla_ref)
		shift
		hla_ref=$(parse_path "$1")
		;;
	--sample)
		shift
		sample="$1"
		;;
	--outdir)
		shift
		outdir="$1"
		;;
	-j | --nproc)
		shift
		nproc="$1"
		;;
	-p | --nproc_per_job)
		shift
		nproc_per_job="$1"
		;;
	--no_gunzip)
		no_gunzip=true
		;;
	--overwrite)
		overwrite=true
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

if [ -z "$r1" ]; then
	error "$0" "$LINENO" "--r1 is required"
	usage
	exit 1
fi

if [ -z "$r2" ]; then
	error "$0" "$LINENO" "--r2 is required"
	usage
	exit 1
fi

if [ -z "$hla_ref" ]; then
	error "$0" "$LINENO" "--hla_ref is required"
	usage
	exit 1
else
	hla_ref_nix="${hla_ref%.*}.nix"
	check_file_exists "$hla_ref_nix"
fi

if [ -z "$sample" ]; then
	error "$0" "$LINENO" "--sample is required"
	usage
	exit 1
fi

if [ -z "$outdir" ]; then
	error "$0" "$LINENO" "--outdir is required"
	usage
	exit 1 
fi

info "$0" ${LINENO} "Run RazerS3 HLA realigner" 
start_time=$(date +%s)

outdir=$( make_dir "$outdir" )
realn_bam="$outdir/$sample.hla.realigner.bam"
donefile="$outdir/$sample.hla.realigner.done"
runtime_file="$outdir/$sample.realigner.runtime.tsv"
if [ "$overwrite" = true ];
then
	rm -f "$realn_bam" "$donefile"
fi

if [ -f "${donefile}" ] && [ -f "${realn_bam}" ];
then
	info "$0" ${LINENO} "Previous HLA realignment result exists: $realn_bam"
	info "$0" ${LINENO} "Run RazerS3 HLA realigner [DONE]" 
	exit 0
fi

# 1. fish HLA-relevant read candidates using razers3
# 1.1 split input fastq files using seqkit split2
info "$0" "$LINENO" "Split Fastq file to prepare for razerS3 realignment"
splits_dir="$outdir/splits"
splits_dir=$( make_dir "$splits_dir" )
run_seqkit_split2 "$r1" "$r2" "$splits_dir" "$nproc"
info "$0" "$LINENO" "Split Fastq file to prepare for razerS3 realignment [DONE]"

# 1.2 run razerS3 realignment on each individually split fastq files
info "$0" "$LINENO" "Fish HLA reads using razerS3 realignment"
razer_dir="$outdir/razer"
fq_search_regex="*.part_*.fastq.gz"
njobs=$(( "$nproc" / "$nproc_per_job" ))
run_razers3_batch "$splits_dir" "$razer_dir" "$hla_ref" "$fq_search_regex" "$njobs" "$nproc_per_job"
info "$0" "$LINENO" "Fish HLA reads using razerS3 realignment [DONE]"

# 1.3 extract fished reads from BAM aligned by razerS3
info "$0" "$LINENO" "Extract fished reads"
bam2fq_dir="$outdir/bam2fq"
bam2fq_dir=$( make_dir "$bam2fq_dir" )
bam_search_regex="*.part_*.razers3.bam"
run_bam2fq_batch "$razer_dir" "$bam2fq_dir" "$bam_search_regex" "$nproc"
info "$0" "$LINENO" "Extract fished reads [DONE]"

# 1.4 pair fished reads per split part
# I am not concatenating and split again as what I did before
# the regex belwow only gives the suffix part
info "$0" "$LINENO" "Pair fished reads"
pair_dir="$outdir/pair"
pair_dir=$( make_dir "$pair_dir" )
fq_search_regex=".part_*razers3.fastq.gz"
run_seqkit_pair_batch "$bam2fq_dir" "$pair_dir" "$fq_search_regex" "$nproc"
info "$0" "$LINENO" "Pair fished reads [DONE]"

# --no_gunzip
# two-layers of if... not pretty
if [ "$no_gunzip" = false ]; then
	info "$0" "$LINENO" "Gunzipping paired reads prior to Novoalign"
	gunzip_donefile="$pair_dir/log/gunzip.done"
	if [ ! -f "$gunzip_donefile" ]; then
		paired_fq_files=$( find_files_on_pattern "$pair_dir" "*.razers3.fastq.gz")
		if [ -z "$paired_fq_files" ]; then
			error "$0" "$LINENO" "Failed to find any paired fastq.gz for gunzip"
			exit 1
		fi
		echo "${paired_fq_files[@]}" | \
			xargs -P"$nproc" -I{} gunzip {}
		info "$0" "$LINENO" "Gunzipping paired reads prior to Novoalign [DONE]"
		touch "$gunzip_donefile"
	else
		info "$0" "$LINENO" "Paired reads were already gunzipped previously"
	fi
fi

# 1.5 realign fished reads using novoalign
info "$0" "$LINENO" "Realign paired reads to HLA using Novoalign"
novo_dir="$outdir/novoalign"
novo_dir=$( make_dir "$novo_dir" )
fq_search_regex=
# when no need for gunzipping, we use the fastq.gz as pattern
if [ "$no_gunzip" = false ]; then
	fq_search_regex=".part_*razers3.fastq"
else
	fq_search_regex=".part_*razers3.fastq.gz"
fi
run_novoalign_batch "$pair_dir" "$novo_dir" "$fq_search_regex" "$hla_ref_nix" "$nproc"
info "$0" "$LINENO" "Realign paired reads to HLA using Novoalign [DONE]"

# 1.5 concatenate individual bam files to get the merged realigned BAM
info "$0" ${LINENO} "Concatenate individually novoaligned BAM files" 
bam_search_regex="*.novoalign.bam"
cat_bam="$novo_dir/$sample.novoalign.cat.bam"
run_samtools_cat "$novo_dir" "$bam_search_regex" "$cat_bam"
info "$0" ${LINENO} "Concatenate individually novoaligned BAM files [DONE]" 

# 1.6 sort the merged realigned BAM
info "$0" ${LINENO} "Post-process realigned BAM file" 
run_samtools_sort "$cat_bam" "$realn_bam" "$nproc"
info "$0" ${LINENO} "Post-process realigned BAM file [DONE]" 
info "$0" ${LINENO} "Run Polysolver HLA realigner [DONE]" 

end_time=$(date +%s)
#runtime=$( echo "${end_time} - ${start_time}" | bc -l )
runtime=$(( "$end_time" - "$start_time" ))
echo -e "${sample}\t$(date -u -d @"${runtime}" +'%M.%S')m" > "$runtime_file" 

# the paired fastq files are still useful for hlafinalizer.sh
rm -rf "$splits_dir" "$razer_dir" "$bam2fq_dir" "$novo_dir"

touch "$donefile"
