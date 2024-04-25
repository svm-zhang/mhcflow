#!/usr/bin/env bash

set -e

SRC_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &>/dev/null && pwd)
#COMMON_FUNC_LIB_DIR="${SRC_DIR%/*}/lib"
COMMON_FUNC_LIB="${SRC_DIR}/hla_bash_util_funcs"
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
	-o or --out    & Specify the path to the output BAM file (Required)
	--nproc    & Specify the number of CPUs used [8]
EO
}

r1=
r2=
hla_ref=
hla_ref_nix=
sample=
out=
skip_fish=false
nproc=8
nproc_per_job=2

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
	-o | --out)
		shift
		out="$1"
		;;
	--skip-fish)
		skip_fish=false
		;;
	-j | --nproc)
		shift
		nproc="$1"
		;;
	-p | --nproc_per_job)
		shift
		nproc_per_job="$1"
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

if [ -z "$out" ]; then
	error "$0" "$LINENO" "--out is required"
	usage
	exit 1 
fi

info "$0" ${LINENO} "Run RazerS3 HLA realigner" 
start_time=$(date +%s)

wkdir=${out%/*}
wkdir="${wkdir}/${sample}_hla_realn"
wkdir=$(make_dir "$wkdir")
done="${wkdir}/${sample}.hla.realn.done"
runtime_file="$wkdir/$sample.realn.runtime.tsv"
if [ -f "${done}" ] && [ -f "${out}" ];
then
	info "$0" ${LINENO} "Previous HLA realignment result exists. Skip realignment..."
	info "$0" ${LINENO} "Run RazerS3 HLA realigner [DONE]" 
	exit 0
fi

# 1. fish HLA-relevant read candidates using razers3
# 1.1 split input fastq files using seqkit split2
info "$0" "$LINENO" "Split Fastq file to prepare for razerS3 realignment"
wip_dir="$wkdir/wip"
wip_dir=$( make_dir "$wip_dir" )
run_seqkit_split2 "$r1" "$r2" "$wip_dir" "$nproc"

info "$0" "$LINENO" "Split Fastq file to prepare for razerS3 realignment [DONE]"
# FIXME: check if nproc * 2 files generated, if not, exit
# FIXME: the number of output files should always be an even number

# 1.2 run razerS3 realignment on each individually split fastq files
info "$0" "$LINENO" "Fish HLA reads using razerS3 realignment"
fq_search_regex="*.part_*.fastq.gz"
run_razers3_batch "$wip_dir" "$hla_ref" "$fq_search_regex" "$nproc" "$nproc_per_job"

# 1.3 extract fished reads from BAM aligned by razerS3
info "$0" "$LINENO" "Extract fished reads"
bam_search_regex="*.part_*.razers3.bam"
run_bam2fq_batch "$wip_dir" "$bam_search_regex" "$nproc"
info "$0" "$LINENO" "Extract fished reads [DONE]"

# 1.4 pair fished reads per split part
# I am not concatenating and split again as what I did before
# the regex belwow only gives the suffix part
info "$0" "$LINENO" "Pair fished reads"
fq_search_regex=".part_*razers3.fastq.gz"
gunzip=true
run_seqkit_pair_batch "$wip_dir" "$fq_search_regex" "$gunzip" "$nproc"
info "$0" "$LINENO" "Pair fished reads [DONE]"

info "$0" "$LINENO" "Realign paired reads to HLA using Novoalign"
fq_search_regex=".part_*razers3.paired.fastq"
run_novoalign_batch "$wip_dir" "$fq_search_regex" "$hla_ref_nix" "$nproc"
info "$0" "$LINENO" "Realign paired reads to HLA using Novoalign [DONE]"

info "$0" ${LINENO} "Concatenate individually novoaligned BAM files" 
realn_bam="${wkdir}/${sample}.hla.realn.bam"
bam_search_regex="*.razers3.paired.bam"
bam_input_str=$( find_files_on_pattern "$wip_dir" "$bam_search_regex")
bam_list_file="$wip_dir/bams.list.txt"
echo "${bam_input_str[@]}" > "$bam_list_file"
cmd=(
	"samtools"
	"cat"
	"-o"
	"$realn_bam"
	"-b"
	"$bam_list_file"
)
if ! "${cmd[@]}" ;
then
	error "$0" "$LINENO" "Failed to concatenate novoalign-aligned BAM files"
	exit 1
fi
info "$0" ${LINENO} "Concatenate individually novoaligned BAM files [DONE]" 

info "$0" ${LINENO} "Post-process realigned BAM file" 
cmd=(
	"bash"
	"${SRC_DIR}/bamer.sh"
	"--bam"
	"$realn_bam"
	"--out"
	"$out"
)
if ! "${cmd[@]}" ;
then
	error "$0" "$LINENO" "Failed to post-process realigned BAM file"
fi
info "$0" ${LINENO} "Post-process realigned BAM file [DONE]" 

info "$0" ${LINENO} "Run Polysolver HLA realigner [DONE]" 

end_time=$(date +%s)
runtime=$( echo "${end_time} - ${start_time}" | bc -l )
echo -e "${sample}\t$(date -u -d @"${runtime}" +'%M.%S')m" > "$runtime_file" 

rm -rf "${wip_dir}"

touch "$done"

exit 0

# 1.3 cat individual fished reads
# FIXME: this does not return gzipped fastq
# FIXME: lets not cat and split again
fished_r1_fq="$wip_dir/$sample.fished.R1.fastq.gz"
fished_r2_fq="$wip_dir/$sample.fished.R2.fastq.gz"
cmd="find $wip_dir -name '*R1*.fished.fastq.gz' -exec cat {} \+ > $fished_r1_fq "
run_cmd "$cmd" "$LINENO" "Failed to concatenate fished R1 reads"
cmd="find $wip_dir -name '*R2*.fished.fastq.gz' -exec cat {} \+ > $fished_r2_fq "
run_cmd "$cmd" "$LINENO" "Failed to concatenate fished R2 reads"

# 1.4 pair fished reads
info "$0" "$LINENO" "Pair fished reads"
cmd="seqkit pair -1 $fished_r1_fq -2 $fished_r2_fq >/dev/null 2>&1"
run_cmd "$cmd" "$LINENO" "Failed to pair fished reads"
info "$0" "$LINENO" "Pair fished reads [DONE]"

# 1.5 realign fished reads using novoalign
realn_bam="${wkdir}/${sample}.hla.realn.bam"
if [ ! -f "${out}" ]; then
	info "$0" "$LINENO" "Realign paired reads to HLA using Novoalign"
	paired_r1_fq="${fished_r1_fq%.*.*}.paired.fastq.gz"
	paired_r2_fq="${fished_r2_fq%.*.*}.paired.fastq.gz"
	check_file_exists "$paired_r1_fq"
	check_file_exists "$paired_r2_fq"
	fqhead=$( zcat "$paired_r1_fq" | head )
	if [ -n "$fqhead" ];
	then
		# split paired 
		cmd="seqkit split2 -p $nproc -1 $paired_r1_fq -2 $paired_r2_fq -O $wip_dir -j $nproc 2>/dev/null"
		run_cmd "$cmd" "$LINENO" "Failed to split Fastq files" 
		exit 0

		cmd="find $wip_dir -name '*R1*paired.part_*.fastq.gz' | \
			xargs -n1 -P$nproc bash -c 'gunzip -f \"\$@\"' _ $1"
		run_cmd "$cmd" "$LINENO" "Failed to decomparess R1 paired gzipped Fastq files" 
		cmd="find $wip_dir -name '*R2*paired.part_*.fastq.gz' | \
			xargs -n1 -P$nproc bash -c 'gunzip -f \"\$@\"' _ $1"
		run_cmd "$cmd" "$LINENO" "Failed to decomparess R2 paired gzipped Fastq files" 

		cmd="paste <(find $wip_dir -name '*R1*paired.part_*.fastq' | sort ) <(find $wip_dir -name '*R2*paired.part_*.fastq' | sort ) | \
			awk '{print \$1\"\n\"\$2}' | \
			xargs -n2 -P$nproc bash -c 'run_novoalign \"\$@\"' $LIBREALIGN $1 $2 $hla_ref_nix"
		run_cmd "$cmd" "$LINENO" "Failed to realign paired reads to HLA using Novoalign"

		cmd="samtools cat -o $realn_bam $wip_dir/*.paired.part_*.bam"
		run_cmd "$cmd" "$LINENO" "Failed to concatenate inidividual bams"
		info "$0" "$LINENO" "Realign paired reads to HLA using Novoalign [DONE]"

	else
		error "$0" "$LINENO" "Found no fished reads that are in pairs."
		exit 1
	fi
else
	info "$0" "$LINENO" "Found previous realignment results: $out. Skip..."
fi

info "$0" ${LINENO} "Post-process realigned BAM file" 
cmd="bash ${SRC_DIR}/bamer.sh --bam $realn_bam --out ${out}"
run_cmd "$cmd" "$LINENO" "Failed to post-process realigned BAM file"
info "$0" ${LINENO} "Post-process realigned BAM file [DONE]" 

info "$0" ${LINENO} "Run Polysolver HLA realigner [DONE]" 

end_time=$(date +%s)
runtime=$( echo "${end_time} - ${start_time}" | bc -l )
runtime_file="$wkdir/$sample.realn.runtime.tsv"
echo -e "${sample}\t$(date -u -d @${runtime} +'%M.%S')m" > "$runtime_file" 

rm -rf "${wip_dir}"

touch "$done"

exit 0
