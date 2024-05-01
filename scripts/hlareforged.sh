#!/usr/bin/env bash

# script to run hla realigner, typer and finalizer from end to end

SRC_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &>/dev/null && pwd)
LIBCOMMON="${SRC_DIR%/*}/lib/libcommon.sh"
LIBREALIGN="${SRC_DIR%/*}/lib/librealign.sh"
source "$LIBCOMMON"
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
EO
}

r1=
r2=
hla_ref=
freq_file=
race="Unknown"
sample=
outdir=
mdup_ram=8
nproc=8
nproc_per_job=2
realn_only=false
overwrite=false
no_clean=false

if [ "$#" -le 1 ]; then
	usage
	exit 1
fi

while [ $# -gt 0 ]; do
	case $1 in
	-h | --help)
		usage; exit 0 ;;
	--r1)
		shift; r1=$(parse_path "$1") ;;
	--r2)
		shift; r2=$(parse_path "$1") ;;
	--hla_ref)
		shift; hla_ref=$(parse_path "$1") ;;
	--sample)
		shift; sample="$1" ;;
  --freq)
    shift; freq_file=$(parse_path "$1");;
	--outdir)
		shift; outdir="$1" ;;
	-j | --nproc)
		shift; nproc="$1" ;;
	-p | --nproc_per_job)
		shift; nproc_per_job="$1" ;;
  --mdup_ram)
    shift; mdup_ram="$1" ;;
	--realn_only)
		realn_only="$1";;
	--overwrite)
    overwrite=true ;;
	--no_clean)
    no_clean=true ;;
	--)
		shift; break;;
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

# 1.1 run realigner
outdir=$( make_dir "$outdir" )
realn_dir="$outdir/realigner"
typer_dir="$outdir/typer"
finalizer_dir="$outdir/finalizer"
cmd=(
  "razer_realigner"
  "--sample"
  "$sample"
  "--r1"
  "$r1"
  "--r2"
  "$r2"
  "--hla_ref"
  "$hla_ref"
  "--outdir"
  "$realn_dir"
  "--nproc"
  "$nproc"
  "--nproc_per_job"
  "$nproc_per_job"
)
if ! "${cmd[@]}"; then
  exit 1
fi

if [ "$realn_only" = true ]; then
	info "$0" "$LINENO" "Realigner only mode was specified. Quit now"
	exit 0
fi

realn_bam=$( find_files_on_pattern "$realn_dir" "${sample}*realigner.bam" )
if [ -z "$realn_bam" ]; then
  error "$0" "$LINENO" "Found no realignment result"
  exit 1
fi
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
  exit 1
fi

typeres=$( find_files_on_pattern "$typer_dir" "${sample}.hla_typing.tsv" )
if [ -z "$typeres" ]; then
  error "$0" "$LINENO" "Found no HLA typing result"
  exit 1
fi

cmd=(
  "hlafinalizer"
  "--sample"
  "$sample"
  "--realn_dir"
  "$realn_dir"
  "--typeres"
  "$typeres"
  "--hla_ref"
  "$hla_ref"
  "--outdir"
  "$finalizer_dir"
  "--nproc"
  "$nproc"
  "--mdup_ram"
  "$mdup_ram"
)
if ! "${cmd[@]}"; then
  exit 1
fi

exit 0
