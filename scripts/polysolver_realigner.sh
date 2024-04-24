#!/usr/bin/env bash

set -e


SRC_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
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
	-b or --bam    & Specify the path to the BAM file [Required]	
	-t or --tag    & Specify the TAG file, e.g. abc_v14.uniq [Required]
	--nv_idx    & Specify the Novoalign HLA indexed file, e.g. abc_complete.nix [Required]
	-g or --genome    & Specify the path to the genome (through which .fai will be used) [Required]	
	--sample    & Specify the sample name [Required]
	--bed    & Specify the path to the HLA region defined in BED [Required]
	-o or --out    & Specify the path to the reaglined BAM file [Required]
	--nproc    & Specify the number of CPUs used [8]
EO
}

bam=
genome=
fai=
sample=
tag_file=
nv_idx=
hla_bed=
out=
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
		-b|--bam)
			shift; bam=$(parse_path "$1");;
		-t|--tag)
			shift; tag_file=$(parse_path "$1");;
		-g|--genome)
			shift; genome=$(parse_path "$1");;
		--nv_idx)
			shift; nv_idx=$(parse_path "$1");;
		--sample)
			shift; sample="$1";;
		--bed)
			shift; hla_bed=$(parse_path "$1");;
		-o|--out)
			shift; out="$1";;
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

if [ -z "$bam" ];
then
	error "$0" ${LINENO} "polysolver requires a BAM file to work with" 
	usage
	exit 1
fi

if [ -z "$genome" ];
then
	error "$0" ${LINENO} "polysolver requires a genome" 
	usage
	exit 1
else
	fai="$genome.fai"
	check_file_exists "$fai"
fi

if [ -z "$sample" ];
then
	error "$0" ${LINENO} "polysolver needs a sample name" 
	usage
	exit 1
fi

if [ -z "$tag_file" ];
then
	error "$0" ${LINENO} "polysolver requires a TAG sequence file, such as abc_v14.uniq" 
	usage
	exit 1
fi

if [ -z "$nv_idx" ];
then
	error "$0" ${LINENO} "polysolver requires the Novoalign HLA index, such as abc_complete.nix" 
	usage
	exit 1
fi

if [ -z "$hla_bed" ];
then
	error "$0" ${LINENO} "polysolver requires the HLA region file in BED" 
	usage
	exit 1
fi

info "$0" ${LINENO} "Run Polysolver HLA realigner" 
start_time=$(date +%s)

wkdir=${out%/*}
wkdir="${wkdir}/${sample}_hla_realn"
wkdir=$(make_dir "$wkdir")
done="${wkdir}/${sample}.hla.realn.done"
if [ -f "${done}" ] && [ -f "${out}" ];
then
	info "$0" ${LINENO} "Previous HLA realignment result exists. Skip realignment..."
	info "$0" ${LINENO} "Run Polysolver HLA realigner [DONE]" 
	exit 0
fi

# getting matching tag sequences
info "$0" ${LINENO} "Fish reads with exact TAG sequences in BAM" 
tag_read_ids="${wkdir}/${sample}.tag.ids"
cmd="samtools view $bam | grep -F -f $tag_file | cut -f1 | sort | uniq > $tag_read_ids"
if ! eval "$cmd";
then
	error "$0" ${LINENO} "Failed to extract reads with exact TAG sequences: ${cmd}"
fi
info "$0" ${LINENO} "Fish reads with exact TAG sequences in BAM [DONE]" 

#getting chr6 region
info "$0" ${LINENO} "Get reads mapped to HLA regions in BAM" 

hla_read_ids="${wkdir}/${sample}.hla.ids"
cmd="samtools view -ML $hla_bed $bam | cut -f1 | sort | uniq > $hla_read_ids"
if ! eval "$cmd";
then
	error "$0" ${LINENO} "Failed to extract reads mapped to HLA regions: ${cmd}"
fi
info "$0" ${LINENO} "Get reads mapped to HLA regions in BAM [DONE]" 

info "$0" ${LINENO} "Collect both TAG and HLA reads" 
hla_read_merged_ids="${wkdir}/${sample}.merged.ids"
cat "$tag_read_ids" "$hla_read_ids" | sort | uniq > "$hla_read_merged_ids"

merged_R1="${wkdir}/${sample}.merged.R1.fastq"
merged_R2="${wkdir}/${sample}.merged.R2.fastq"
samtools view "$bam" | grep -F -f "$hla_read_merged_ids" | samtools view -bt "$fai" | \
	samtools sort -n | \
	samtools fastq -1 "$merged_R1" -2 "$merged_R2" -0 /dev/null -s /dev/null -n

info "$0" ${LINENO} "Collect both TAG and HLA reads [DONE]" 

info "$0" ${LINENO} "Realign reads to HLA using Novoalign" 
# below gives the number of lines per file to proc.
nreads_per_proc=$(wc -l "$merged_R1" | awk -v nproc="$nproc" '{print (int(($1/4)/nproc)+1)*4}')
split_dir="${wkdir}/splits"
split_fq_dir=$(make_dir "${split_dir}/fqs")
split_bam_dir=$(make_dir "${split_dir}/bams")
split_log_dir=$(make_dir "${split_bam_dir}/log")

# split R1 into individual files
split -l "$nreads_per_proc" -d -a 2 --additional-suffix ".R1.fastq" "$merged_R1" "${split_fq_dir}/"
split -l "$nreads_per_proc" -d -a 2 --additional-suffix ".R2.fastq" "$merged_R2" "${split_fq_dir}/"
find "$split_fq_dir" -name "*.R1.fastq" | \
	awk '{n=split($1, a, "/"); print a[n]}' | \
	sed 's/\.R1\.fastq//g' | \
	xargs -P"$nproc" -I{} bash -c "novoalign -d ${nv_idx} -F STDFQ -R 0 -r All -o SAM -o FullNW -f ${split_fq_dir}/{}.R1.fastq ${split_fq_dir}/{}.R2.fastq 2>${split_log_dir}/${sample}.nv.{}.log | samtools view -bh -o ${split_bam_dir}/${sample}.{}.bam  && touch ${split_log_dir}/${sample}.{}.done"
info "$0" ${LINENO} "Realign reads to HLA using Novoalign [DONE]" 

hla_realign_bam="${wkdir}/${sample}.hla.realn.bam"
info "$0" ${LINENO} "Concatenate individual split BAM files" 
cmd="samtools cat -o $hla_realign_bam $split_bam_dir/*.bam"
if ! eval "$cmd";
then
	error "$0" ${LINENO} "Failed to concatenate split BAM filesi: ${cmd}"
	exit 1
fi
info "$0" ${LINENO} "Concatenate individual split BAM files [DONE]" 

info "$0" ${LINENO} "Post-process realigned BAM file" 
cmd="bash ${SRC_DIR}/bamer.sh --rmdup --bam ${hla_realign_bam} --out ${out}"
if ! eval "$cmd";
then
	error "$0" ${LINENO} "Failed to post-process realigned BAM file"
	exit 1
fi
info "$0" ${LINENO} "Post-process realigned BAM file [DONE]" 

info "$0" ${LINENO} "Run Polysolver HLA realigner [DONE]" 
end_time=$(date +%s)
runtime=$( echo "${end_time} - ${start_time}" | bc -l )
runtime_file="${wkdir}/${sample}.realn.runtime.tsv"
echo -e "${sample}\t$(date -u -d @${runtime} +'%M.%S')m" > "${runtime_file}" 

rm -r "${split_dir}"

touch "$done"

exit 0