#!/usr/bin/env bash

# script to run hla realigner, typer and finalizer from end to end

SRC_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &>/dev/null && pwd)
LIBCOMMON="${SRC_DIR%/*}/lib/libcommon.sh"
source "$LIBCOMMON"

function usage() {
  local program
  program="${0##*/}"
  cat <<EO
Usage: $program [options]
Options:
EO
  cat <<EO | column -s\& -t
	--sample    & Specify the sample name [Required]
	--bam    & Specify the path to the BAM file [Required]	
	--tag    & Specify the TAG file, e.g. abc_v14.uniq [Required]
	--hla_ref    & Specify the Novoalign HLA indexed file, e.g. abc_complete.fasta [Required]
	--bed    & Specify the path to the HLA region defined in BED [Required]
	--freq    & Specify the HLA allele population frequency file [Required]
	--fish_mode    & Specify the mode fisher to use (faster, fast) [faster]
	--min_ecnt    & Specify the min number of mismatch event allowed [999]
	--outdir    & Specify the path to output base directory [Required]
	--nproc    & Specify the number of CPUs used [8]
EO
}

bam=
outdir=
tag_file=
hla_ref=
hla_bed=
freq_file=
fish_mode="faster"
min_ecnt=999
sample=
nproc=8
realn_only=false

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
  --sample)
    shift
    sample="$1"
    ;;
  --bam)
    shift
    bam="$1"
    ;;
  --tag)
    shift
    tag_file="$1"
    ;;
  --hla_ref)
    shift
    hla_ref="$1"
    ;;
  --bed)
    shift
    hla_bed="$1"
    ;;
  --fish_mode)
    shift
    fish_mode="$1"
    ;;
  --freq)
    shift
    freq_file="$1"
    ;;
  --min_ecnt)
    shift
    min_ecnt="$1"
    ;;
  --outdir)
    shift
    outdir="$1"
    ;;
  -p | --nproc)
    shift
    nproc="$1"
    ;;
  --realn_only) realn_only=true ;;
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

outdir=$(make_dir "$outdir")
logdir="$outdir/log"
logdir=$(make_dir "$logdir")
donefile="$logdir/$sample.hlatyping.done"

if [ -f "$donefile" ]; then
  info "$0" "$LINENO" "Found done file from previous run $donefile"
  info "$0" "$LINENO" "Remove done file if you want to re-run"
  exit 0
fi

# fishing
fish_dir="$outdir/fisher"
fish_out="$fish_dir/$sample.fished.fqs.txt"
fisher --mode "$fish_mode" --tag "$tag_file" --bed "$hla_bed" --bam "$bam" \
  --sample "$sample" --out "$fish_out" --nproc "$nproc" ||
  die "$0" "$LINENO" "Failed to run fisher"

check_file_exists "$fish_out" ||
  die "$0" "$LINENO" "Failed to find fisher result $fish_out"

# realigner
if [ "$realn_only" = false ]; then
  realn_dir="$outdir/realigner"
  realn_out="$realn_dir/$sample.hla.realn.so.bam"
  realigner --hla_ref "$hla_ref" --fqs "$fish_out" \
    --sample "$sample" --out "$realn_out" --nproc "$nproc" ||
    die "$0" "$LINENO" "Failed to run realigner"
else
  realn_dir="$outdir/finalizer"
  realn_out="$realn_dir/$sample.hla.realn.ready.bam"
  realigner --hla_ref "$hla_ref" --fqs "$fish_out" \
    --sample "$sample" --out "$realn_out" --nproc "$nproc" --mdup ||
    die "$0" "$LINENO" "Failed to run realigner"

  info "$0" "$LINENO" "Realigner-only mode was specified"
  info "$0" "$LINENO" "Will not continue typing"
  exit 0
fi

# typer
check_file_exists "$realn_out" ||
  die "$0" "$LINENO" "Failed to find realigned BAM file $realn_out"
typer_dir="$outdir/typer"
typer_out="$typer_dir/$sample.hlatyping.res.tsv"
mhctyper --freq "$freq_file" --outdir "$typer_dir" \
  --bam "$realn_out" --min_ecnt "$min_ecnt" --nproc "$nproc" ||
  die "$0" "$LINENO" "Failed to run typer."
check_file_exists "$typer_out" ||
  die "$0" "$LINENO" "Failed to find typing result $typer_out"

# extractor
final_dir="$outdir/finalizer"
extract_out="$final_dir/$sample.hla.fasta"
extractor --hla_ref "$hla_ref" \
  --sample "$sample" --typeres "$typer_out" --out "$extract_out" ||
  die "$0" "$LINENO" "Failed to extract sample HLA ref"
check_file_exists "$extract_out" ||
  die "$0" "$LINENO" "Failed to find extraced sample HLA ref: $extract_out"

# realigner against the sample hla ref
final_realn_out="$final_dir/$sample.hla.realn.ready.bam"
realigner --hla_ref "$extract_out" --fqs "$fish_out" \
  --sample "$sample" --out "$final_realn_out" \
  --mdup --nproc "$nproc" ||
  die "$0" "$LINENO" \
    "Failed to run realigner on sample HLA ref"

touch "$donefile"

exit 0
