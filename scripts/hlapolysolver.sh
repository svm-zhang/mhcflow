#!/usr/bin/env bash

# script to run hla realigner, typer and finalizer from end to end

SRC_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &>/dev/null && pwd)
LIBCOMMON="${SRC_DIR%/*}/lib/libcommon.sh"
source "$LIBCOMMON"

function usage () {
	local program
  program="${0##*/}"
	cat << EO
Usage: $program [options]
Options:
EO
	cat << EO | column -s\& -t
	-b or --bam    & Specify the path to the BAM file [Required]	
	-t or --tag    & Specify the TAG file, e.g. abc_v14.uniq [Required]
	--nv_idx    & Specify the Novoalign HLA indexed file, e.g. abc_complete.nix [Required]
	--sample    & Specify the sample name [Required]
	--bed    & Specify the path to the HLA region defined in BED [Required]
	--freq    & Specify the HLA allele population frequency file [Required]
	-o or --out    & Specify the path to the reaglined BAM file [Required]
	-r or --race    & Specify the race (Caucasian, Black, Asian, Unknown) [Unknown]
	--nproc    & Specify the number of CPUs used [8]
EO
}

bam=
out_bam=
tag_file=
nv_idx=
hla_bed=
freq_file=
race="Unknown"
sample=
nproc=8

if [ "$#" -le 1 ]; then
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
		-b|--bam)
			shift; bam=$(parse_path "$1");;
		-t|--tag)
			shift; tag_file=$(parse_path "$1");;
		--nv_idx)
			shift; nv_idx=$(parse_path "$1");;
		--bed)
			shift; hla_bed=$(parse_path "$1");;
		--freq)
			shift; freq_file=$(parse_path "$1");;
		--race)
			shift; race="$1";;
		--out_bam)
			shift; out_bam="$1";;
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

# 1.1 run realigner
outdir=${out_bam%/*}
outdir=$( make_dir "$outdir" )
cmd=(
  "polysolver_realigner"
  "--sample"
  "$sample"
  "--bam"
  "$bam"
  "--nv_idx"
  "$nv_idx"
  "--tag"
  "$tag_file"
  "--bed"
  "$hla_bed"
  "--out"
  "$out_bam"
  "--nproc"
  "$nproc"
)
if ! "${cmd[@]}"; then
  exit 1
fi

if [ ! -f "$out_bam" ]; then
  error "$0" "$LINENO" "Failed to find realigned BAM file: $out_bam"
  error "$0" "$LINENO" "Please check realigner result under $outdir"
  exit 1
fi

# 1.2 run typer
cmd=(
  "hlatyper"
  "--sample"
  "$sample"
  "--bam"
  "$out_bam"
  "--freq"
  "$freq_file"
  "--race"
  "$race"
  "--outdir"
  "$outdir"
  "--nproc"
  "$nproc"
)
if ! "${cmd[@]}"; then
  exit 1
fi

exit 0
