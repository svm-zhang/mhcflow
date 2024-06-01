#!/usr/bin/env bash

set -e
set -o pipefail


SRC_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
LIBCOMMON="${SRC_DIR%/*}/lib/libcommon.sh"
source "$LIBCOMMON"

function run_novoalign () {
  local reads="$1"
  local sample="$2"
  local nix="$3"
  local outdir="$4"

  IFS=" " read -r r1 r2 <<< "$reads"

  local logdir="$outdir/log"
  logdir=$( make_dir "$logdir" )
  local prefix="${r1##*/}"
  local logfile="$logdir/${prefix%.fastq}.novoalign.log"

  local bam="$outdir/${prefix%.fastq}.novoalign.bam"
  info "${FUNCNAME[0]}" "$LINENO" "Run Novoalign on ${r1##*/} ${r2##*/}"
  rg_str="@RG\tID:${sample}\tSM:${sample}"
	novoalign -d ${nix} -F STDFQ -R 0 -r All -o SAM "$rg_str" -o FullNW \
		-f "$r1" "$r2" 2>"$logfile" \
		| samtools view -bh -o "$bam" \
		|| die "$0" "$LINENO" "Failed to novoalign on $r1 and $r2"
}
export -f run_novoalign

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
	-o or --out    & Specify the path to the reaglined BAM file [Required]
	--nproc    & Specify the number of CPUs used [8]
EO
}

#bam=
fqs=
sample=
#tag_file=
hla_ref=
#hla_bed=
out_bam=
outdir=
nproc=8

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
		--fqs)
			shift; fqs=$(parse_path "$1");;
		#-b|--bam)
		#	shift; bam=$(parse_path "$1");;
		#-t|--tag)
		#	shift; tag_file=$(parse_path "$1");;
		--hla_ref)
			shift; hla_ref=$(parse_path "$1");;
		--sample)
			shift; sample="$1";;
		#--bed)
		#	shift; hla_bed=$(parse_path "$1");;
		#--outdir)
		#	shift; outdir="$1";;
		--out)
			shift; out="$1";;
		-p|--nproc)
			shift; nproc="$1";;
		--)
			shift; break;;
		*)
			echo "Invalid option: $1" 1>&2
			realn_usage; exit 1;
			;;
	esac
	shift
done

#if [ -z "$bam" ]; then
#	error "$0" ${LINENO} "polysolver requires a BAM file to work with" 
#	realn_usage
#	exit 1
#fi

if [ -z "$fqs" ]; then
	error "$0" ${LINENO} "polysolver needs a sample name" 
	realn_usage
	exit 1
fi

if [ -z "$sample" ]; then
	error "$0" ${LINENO} "polysolver needs a sample name" 
	realn_usage
	exit 1
fi

#if [ -z "$tag_file" ]; then
#	error "$0" ${LINENO} "polysolver requires a TAG sequence file, such as abc_v14.uniq" 
#	realn_usage
#	exit 1
#fi

if [ -z "$hla_ref" ]; then
	error "$0" "$LINENO" "polysolver requires the Novoalign HLA fasta, such as abc_complete.fasta"
	realn_usage
	exit 1
fi

#if [ -z "$hla_bed" ]; then
#	error "$0" ${LINENO} "polysolver requires the HLA region file in BED" 
#	realn_usage
#	exit 1
#fi

info "$0" ${LINENO} "Run Polysolver HLA realigner" 
start_time=$(date +%s)

outdir="${out%/*}"
outdir=$( make_dir "$outdir" )
logdir=$( make_dir "$outdir/log")
done="${logdir}/${sample}.hla.realn.done"
if [ -f "${done}" ]; then
	info "$0" ${LINENO} "Previous HLA realignment result exists. Skip realignment..."
	info "$0" ${LINENO} "Run Polysolver HLA realigner [DONE]" 
	exit 0
fi
hla_ref_nix="${hla_ref%.*}.nix"
check_file_exists "$hla_ref_nix"

# getting matching tag sequences
#info "$0" ${LINENO} "Fish reads with exact TAG sequences from BAM" 
#fish_tag_done="${logdir}/${sample}.fish.tag.done"
#tag_read_ids="${outdir}/${sample}.tag.ids"
#unmapped_ids="${outdir}/${sample}.unmapped.tag.ids"
#chr6_ids="${outdir}/${sample}.chr6.tag.ids"
#if [ ! -f "$fish_tag_done" ] || [ ! -f "$tag_read_ids" ]; then
#	info "$0" ${LINENO} "Fish TAG reads from alignments on chromosome 6"
#	samtools view -@"$nproc" "$bam" 6 | \
#		cut -f1,10 | \
#		grep -f "$tag_file" | \
#		cut -f1 > "$chr6_ids"
#	info "$0" ${LINENO} "Fish TAG reads from unplaced alignments"
#	awk 'BEGIN{print "4\n8\n"}' | \
#		xargs -I{} bash -c "samtools view -@$nproc -f{} -F2 $bam | cut -f1,10 | grep -f $tag_file | cut -f1 >> $unmapped_ids"
#	cat "$unmapped_ids" "$chr6_ids" > "$tag_read_ids"
#	touch "$fish_tag_done"
#else
#	info "$0" ${LINENO} "Fish TAG reads previously done. Continue"
#fi
#info "$0" ${LINENO} "Fish reads with exact TAG sequences in BAM [DONE]" 

#bam_size=$( du -bs "$bam" | awk '{print $1/2^30}' )
#if [ -z "$bam_size" ]; then
#	die "$0" "$LINENO" "Failed to obtaint input BAM file size"
#fi
#info "$0" ${LINENO} "Obtain input BAM file size: ${bam_size}G" 
#if [ ! -f "$fish_tag_done" ] || [ ! -f "$tag_read_ids" ]; then
#	info "$0" ${LINENO} "Split BAM into parts for faster fishing reads"
#	split_bam_dir="${outdir}/bams"
#	split_bam_dir=$( make_dir "$split_bam_dir" )
#	bam_size_int=$( printf "%.0f" "$bam_size" )
#	split_size=$(( (bam_size_int+1+nproc-1)/nproc ))
#	info "$0" ${LINENO} "Split BAM file into ${split_size}G each"
#	# split bam files
#	samtools view "$bam" \
#		| cut -f1,10 \
#		| split -C "${split_size}G" -d -a 2 --additional-suffix ".txt" - "$split_bam_dir/" \
#		|| die "$0" "$LINENO" "Failed to split input bam to parts for fishing tags"
#
#	info "$0" ${LINENO} "Fish reads from each individual split file"
#	find "$split_bam_dir" -name "*.txt" \
#		| sed 's/\.txt//g' \
#		| xargs -P"$nproc" -I{} bash -c "grep -F -f $tag_file {}.txt | cut -f1 > {}.ids"
#
#	info "$0" ${LINENO} "Merge individually fished reads"
#	find "$split_bam_dir" -name "*.ids" -exec cat {} + \
#		| sort \
#		| uniq > "$tag_read_ids" \
#		|| die "$0" "$LINENO" "Failed to merged individually fished reads "
#
#	if [ -d "$split_bam_dir" ]; then
#		rm -rf "$split_bam_dir"
#	fi
#	touch "$fish_tag_done"
#else
#	info "$0" ${LINENO} "Fish TAG reads previously done. Continue"
#fi
#info "$0" ${LINENO} "Fish reads with exact TAG sequences in BAM [DONE]" 

#getting chr6 region
#info "$0" ${LINENO} "Get reads mapped to HLA regions in BAM" 
#fish_chr6_done="${logdir}/${sample}.fish.chr6.done"
#hla_read_ids="${outdir}/${sample}.hla.ids"
#if [ ! -f "$fish_chr6_done" ] || [ ! -f "$hla_read_ids" ]; then
#	samtools view -ML "$hla_bed" "$bam" \
#		| cut -f1 \
#		| sort \
#		| uniq > "$hla_read_ids" \
#		|| die "$0" "$LINENO" "Failed to get reads mapped to HLA regions in BAM"
#
#	touch "$fish_chr6_done"
#else
#	info "$0" ${LINENO} "Fish reads from HLA region previously done. Continue"
#fi
#info "$0" ${LINENO} "Get reads mapped to HLA regions in BAM [DONE]" 
#
#info "$0" ${LINENO} "Collect both TAG and HLA reads" 
#hla_read_merged_ids="${outdir}/${sample}.merged.ids"
#if [ ! -f "$hla_read_merged_ids" ]; then
#	cat "$tag_read_ids" "$hla_read_ids" \
#		| sort \
#		| uniq > "$hla_read_merged_ids" \
#		|| die "$0" "$LINENO" "Failed to merge TAG and HLA read IDs"
#else
#	info "$0" ${LINENO} "Found combined fished read ID file: $hla_read_merged_ids"
#fi
#n_fished_reads=$( wc -l < "$hla_read_merged_ids" )
#if(( "$n_fished_reads" == 0 )); then
#	die "$0" "$LINENO" "Fished no reads to continue. Exit"
#fi
## collect the number of fished reads
#fisher_stat_file="$outdir/$sample.fisher.stat.tsv"
#printf "%s\t%s\n" "SampleID" "NumFished" > "$fisher_stat_file"
#printf "%s\t%s\n" "$sample" "$n_fished_reads" >> "$fisher_stat_file"

#merged_R1="${outdir}/${sample}.merged.R1.fastq"
#merged_R2="${outdir}/${sample}.merged.R2.fastq"
#if [ ! -f "$merged_R1" ] || [ ! -f "$merged_R2" ]; then
#	samtools view -@"$nproc" -bh -N "$hla_read_merged_ids" "$bam" \
#		| samtools sort -n \
#		| samtools fastq -n -1 "$merged_R1" -2 "$merged_R2" -0 /dev/null -s /dev/null \
#		|| die "$0" "$LINENO" "Failed to collect TAG na HLA reads based on ID in BAM"
#fi
#info "$0" ${LINENO} "Collect both TAG and HLA reads [DONE]" 

if [ ! -f "$fqs" ]; then
	die "$0" "$LINENO" "Failed to find provided fished fastq list file $fqs"
fi

fished_R1=""
fished_R2=""
fished_R1=$( cut -f1 "$fqs" )
fished_R2=$( cut -f2 "$fqs" )

if [ ! -f "$fished_R1" ] || [ ! -f "$fished_R2" ]; then
	die "$0" "$LINENO" "Failed to find either fished R1 or R2 fastq file"
fi

info "$0" ${LINENO} "Realign reads to HLA using Novoalign" 
# below gives the number of lines per file to proc.
novo_done="${logdir}/${sample}.novoalign.done"
raw_bam="${outdir}/${sample}.hla.realn.bam"
if [ -f "$novo_done" ]; then
	info "$0" ${LINENO} "Realignment has previously done."
	exit 0
fi
	#nreads_per_proc=$(
	#	wc -l "$fished_R1" \
	#	| awk -v nproc="$nproc" '{print (int(($1/4)/nproc)+1)*4}' \
	#	|| die "$0" "$LINENO" "Failed to get number of reads per process to run on"
	#)
	# set it 2 if not properly set
	#if [ -z "$nreads_per_proc" ]; then
	#	nreads_per_proc=2
	#fi

split_dir="${outdir}/splits"
split_log="${logdir}/${sample}.seqkit.split2.log"
seqkit split2 -j"$nproc" -p"$nproc" -O "$split_dir" \
	-1 "$fished_R1" -2 "$fished_R2" > "$split_log" 2>&1 \
	|| die "$0" "$LINENO" "Failed to run seqkit split2 on fished reads"

realn_input_list="${split_dir}/${sample}.realn.input.list.txt"
if [ -f "$realn_input_list" ]; then
	rm "$realn_input_list"
fi
for i in $(seq 1 "$nproc"); do
	find "$split_dir" -name "*00$i*fastq" \
		| sort \
		| paste -sd " " - >> "$realn_input_list" 
done
rg_str="@RG\\tID:${sample}\\tSM:${sample}"
awk '{print}' "$realn_input_list" \
	| xargs -P"$nproc" -I{} bash -c 'run_novoalign "$@"' \
		"run_novoalign" "{}" "$sample" "$hla_ref_nix" "$split_dir"

info "$0" ${LINENO} "Concatenate individual split BAM files"
list_bams_to_cat="$split_dir/.bams_to_cat.list.txt"
find "$split_dir" -name "*.bam" > "$list_bams_to_cat" \
	|| die "$0" "$LINENO" "Failed to find inidividual realigned BAM files"
nbams=$( wc -l <"$list_bams_to_cat" )
if (( "$nbams" == 0 )); then
	die "$0" "$LINENO" "Found zero BAM files under $split_dir"
fi
samtools cat -o "$raw_bam" -b "$list_bams_to_cat" \
	|| die "$0" "$LINENO" "Failed to concatenate individual split BAM files"
info "$0" ${LINENO} "Concatenate individual split BAM files [DONE]"

info "$0" ${LINENO} "Sort realigned BAM file" 
so_tmp="${raw_bam%.bam}.tmp"
so_bam="${outdir}/${sample}.hla.realn.so.bam"
out_bai="${so_bam%.bam}.bai"
# reference: https://github.com/samtools/samtools/issues/1196
samtools sort -T"$so_tmp" -@"$nproc" --write-index \
	-o "$so_bam"##idx##"$out_bai" "$raw_bam" 2>/dev/null \
	|| die "$0" "$LINENO" "Failed to sort reaglined BAM file"
info "$0" ${LINENO} "Sort realigned BAM file [DONE]" 

end_time=$(date +%s)
runtime=$(( end_time - start_time ))
runtime=$( date -u -d @"$runtime" +'%M.%S')
runtime_file="${outdir}/${sample}.realigner.runtime.tsv"
printf "%s\t%s\n" "SampleID" "RealignerTime" > "$runtime_file"
printf "%s\t%s\n" \
	"$sample" \
	"${runtime}m" \
	>> "$runtime_file"

awk '{print}' "$list_bams_to_cat" | xargs -P"$nproc" -I{} rm -f {}
rm -f "$list_bams_to_cat"

touch "$done"
info "$0" ${LINENO} "Realign reads to HLA using Novoalign [DONE]"

exit 0

	#split_fq_dir=$( make_dir "${split_dir}/fqs" )
	#split_bam_dir=$( make_dir "${split_dir}/bams" )
	#split_log_dir=$( make_dir "${split_dir}/log" )

	# split R1 into individual files
	#split -l "$nreads_per_proc" -d -a 2 --additional-suffix ".R1.fastq" "$merged_R1" "${split_fq_dir}/"
	#split -l "$nreads_per_proc" -d -a 2 --additional-suffix ".R2.fastq" "$merged_R2" "${split_fq_dir}/"
  #rg_str="@RG\\tID:${sample}\\tSM:${sample}"
	#find "$split_dir" -name "*.R1*fastq" | \
	#	awk '{n=split($1, a, "/"); print a[n]}' | \
	#	sed 's/\.R1\.fastq//g' | \
	#	xargs -P"$nproc" -I{} bash -c "novoalign -d ${hla_ref_nix} -F STDFQ -R 0 -r All -o SAM \"${rg_str}\" -o FullNW -f ${split_fq_dir}/{}.R1.fastq ${split_fq_dir}/{}.R2.fastq 2>${split_log_dir}/${sample}.nv.{}.log | samtools view -bh -o ${split_bam_dir}/${sample}.{}.bam  && touch ${split_log_dir}/${sample}.{}.done"

	#info "$0" ${LINENO} "Concatenate individual split BAM files"
	#list_bams_to_cat="$split_bam_dir/.bams_to_cat.list.txt"
	#find "$split_bam_dir" -name "*.bam" > "$list_bams_to_cat" \
	#	|| die "$0" "$LINENO" "Failed to find BAM files under $split_bam_dir"
	#nbams=$( wc -l <"$list_bams_to_cat" )
	#if (( "$nbams" == 0 )); then
	#	die "$0" "$LINENO" "Found zero BAM files under $split_bam_dir"
	#fi
	#samtools cat -o "$raw_bam" -b "$list_bams_to_cat" \
	#	|| die "$0" "$LINENO" "Failed to concatenate individual split BAM files"
	#info "$0" ${LINENO} "Concatenate individual split BAM files [DONE]"

	#rm -r "$split_bam_dir" "$split_log_dir"

	#touch "$novo_done"
#else
#	info "$0" ${LINENO} "Found realigned BAM from previous run: $raw_bam"
#fi
#info "$0" ${LINENO} "Realign reads to HLA using Novoalign [DONE]"

#info "$0" ${LINENO} "Sort realigned BAM file" 
#so_tmp="${raw_bam%.bam}.tmp"
#so_bam="${outdir}/${sample}.hla.realn.so.bam"
#out_bai="${so_bam%.bam}.bai"
## reference: https://github.com/samtools/samtools/issues/1196
#if [ ! -f "$so_bam" ] || [ ! -f "$out_bai" ]; then
#	samtools sort -T"$so_tmp" -@"$nproc" --write-index \
#		-o "$so_bam"##idx##"$out_bai" "$raw_bam" 2>/dev/null \
#		|| die "$0" "$LINENO" "Failed to sort reaglined BAM file"
#else
#	info "$0" ${LINENO} "Found sorted and realigned BAM from previous run: $so_bam"
#fi
#info "$0" ${LINENO} "Sort realigned BAM file [DONE]" 

#end_time=$(date +%s)
#runtime=$(( end_time - start_time ))
#runtime=$( date -u -d @"$runtime" +'%M.%S')
#runtime_file="${outdir}/${sample}.realigner.runtime.tsv"
#printf "%s\t%s\t%s\n" "SampleID" "InputSize" "RealignerTime" > "$runtime_file"
#printf "%s\t%s\t%s\n" \
#	"$sample" \
#	"$bam_size" \
#	"${runtime}m" \
#	>> "$runtime_file"
#
#touch "$done"
#
#exit 0
