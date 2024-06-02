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
	--sample    & Specify the sample name [Required]
	-b or --bam    & Specify the path to the BAM file [Required]	
	-t or --tag    & Specify the TAG file, e.g. abc_v14.uniq [Required]
	--hla_ref    & Specify the Novoalign HLA indexed file, e.g. abc_complete.fasta [Required]
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
			realn_only=true;;
		--)
			shift; break;;
		*)
			echo "Invalid option: $1" 1>&2
			usage; exit 1;
			;;
	esac
	shift
done

outdir=$( make_dir "$outdir" )
logdir="$outdir/log"
logdir=$( make_dir "$logdir" )
donefile="$logdir/$sample.hlatyping.done"

if [ -f "$donefile" ]; then
  info "$0" "$LINENO" "Found done file from previous run $donefile"
  info "$0" "$LINENO" "Remove done file if you want to re-run"
	exit 0 
fi

# fishing
fish_dir="$outdir/fisher"
fish_out="$fish_dir/$sample.fished.fqs.txt"
fisher --tag "$tag_file" --bed "$hla_bed" --bam "$bam" \
	--sample "$sample" --out "$fish_out" --nproc "$nproc" \
  || die "$0" "$LINENO" "Failed to run fisher"

if [ ! -f "$fish_out" ]; then
  die "$0" "$LINENO" "Failed to find fisher result $fish_out"
fi

# realigner
realn_dir="$outdir/realigner"
realn_out="$realn_dir/$sample.hla.realn.so.bam"
polysolver_realigner --hla_ref "$hla_ref" --fqs "$fish_out" \
	--sample "$sample" --out "$realn_out" --nproc "$nproc" \
  || die "$0" "$LINENO" "Failed to run realigner"


if [ "$realn_only" = true ]; then
	info "$0" "$LINENO" "Realigner-only mode was specified"
	info "$0" "$LINENO" "Will not continue typing"
	exit 0
fi

# typer
if [ ! -f "$realn_out" ]; then
  die "$0" "$LINENO" "Failed to find realigned BAM file $realn_out"
fi

typer_dir="$outdir/typer"
typer_out="$typer_dir/$sample.hlatyping.res.tsv"
pyhlatyper --freq "$freq_file" --race "$race" --out "$typer_out" \
	--bam "$realn_out" --nproc "$nproc" \
  || die "$0" "$LINENO" "Failed to run typer."
if [ ! -f "$typer_out" ]; then
  error "$0" "$LINENO" "Failed to find typing result $typer_out"
  exit 1
fi

# extractor
final_dir="$outdir/finalizer"
extract_out="$final_dir/$sample.hla.fasta"
extract_sample_hlaref --hla_ref "$hla_ref" \
	--sample "$sample" --typeres "$typer_out" --out "$extract_out" \
  || die "$0" "$LINENO" "Failed to extract sample HLA ref"

# realigner against the sample hla ref
final_realn_out="$final_dir/$sample.hla.realn.ready.bam"
polysolver_realigner --hla_ref "$extract_out" --fqs "$fish_out" \
	--sample "$sample" --out "$final_realn_out" --nproc "$nproc" \
  || die "$0" "$LINENO" \
		"Failed to run realigner on sample HLA ref"

touch "$donefile"

exit 0
