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
	--hla_ref    & Specify the Novoalign HLA indexed file, e.g. abc_complete.fasta [Required]
	--sample    & Specify the sample name [Required]
	--bed    & Specify the path to the HLA region defined in BED [Required]
	--freq    & Specify the HLA allele population frequency file [Required]
	--outdir    & Specify the path to output base directory [Required]
	-r or --race    & Specify the race (Caucasian, Black, Asian, Unknown) [Unknown]
	--nproc    & Specify the number of CPUs used [8]
EO
}

bam=
outdir=
tag_file=
hla_ref=
hla_bed=
freq_file=
race="Unknown"
sample=
nproc=8
realigner="polysolver"
realn_only=false

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
		--hla_ref)
			shift; hla_ref=$(parse_path "$1");;
		--bed)
			shift; hla_bed=$(parse_path "$1");;
		--freq)
			shift; freq_file=$(parse_path "$1");;
		--race)
			shift; race="$1";;
		--outdir)
			shift; outdir="$1";;
		-p|--nproc)
			shift; nproc="$1";;
		--realn_only)
			shift; realn_only=true;;
		--)
			shift; break;;
		*)
			echo "Invalid option: $1" 1>&2
			usage; exit 1;
			;;
	esac
	shift
done

echo "i am here"

# 1.1 run realigner
#outdir=${out_bam%/*}
outdir=$( make_dir "$outdir" )
realn_dir="${outdir}/realigner"
cmd=(
  "polysolver_realigner"
  "--sample"
  "$sample"
  "--bam"
  "$bam"
  "--hla_ref"
  "$hla_ref"
  "--tag"
  "$tag_file"
  "--bed"
  "$hla_bed"
  "--outdir"
  "$realn_dir"
  "--nproc"
  "$nproc"
)
if ! "${cmd[@]}"; then
  die "$0" "$LINENO" "Failed to run polysovler_realigner. Exit"
fi

if [ "$realn_only" = true ]; then
	info "$0" "$LINENO" "Realigner-only mode was specified"
	info "$0" "$LINENO" "Will not continue typing"
	exit 0
fi

# 1.2 run typer
realn_bam=$( find "$realn_dir" -name "*.hla.realn.so.bam" )
if [ -z "$realn_bam" ]; then
  die "$0" "$LINENO" "Failed to find realigned BAM file within $realn_dir"
fi
typer_dir="${outdir}/typer"
cmd=(
  "hlatyper"
  "--sample"
  "$sample"
  "--bam"
  "$realn_bam"
  "--freq"
  "$freq_file"
  "--race"
  "$race"
  "--outdir"
  "$typer_dir"
  "--nproc"
  "$nproc"
)
if ! "${cmd[@]}"; then
  die "$0" "$LINENO" "Failed to run hlatyper. Exit"
fi

# 1.3
typeres=$( find "$typer_dir" -name "${sample}.hla_typing.tsv" )
if [ -z "$typeres" ]; then
  error "$0" "$LINENO" "Found no HLA typing result"
  exit 1
fi
cmd=(
  "hlafinalizer"
  "--sample"
  "$sample"
  "--hla_ref"
  "$hla_ref"
  "--realn_dir"
  "$realn_dir"
  "--typeres"
  "$typeres"
  "--outdir"
  "$outdir"
	"--realigner"
	"$realigner"
  "--nproc"
  "$nproc"
)

if ! "${cmd[@]}"; then
  die "$0" "$LINENO" "Failed to run hlafinalizer. Exit"
fi

exit 0
