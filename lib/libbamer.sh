#!/usr/bin/env bash

set -e
set -o pipefail

SRC_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
LIBCOMMON="$SRC_DIR/libcommon.sh"
source "$LIBCOMMON"

function samtools_sort () {
  local bam="$1"
  local so_bam="$2"
  local nproc="$3"

  info "$0" ${LINENO} "Sort Bam file $bam"
  check_file_exists "$bam" \
    || die "$0" "$LINENO" "Failed to find $bam to sort"

  if [ -z "$nproc" ];
  then
    nproc=1
  fi

  local outdir="${so_bam%/*}"
  local logdir="$outdir/log"
  logdir=$( make_dir "$logdir" )
  local logfile="$logdir/samtools_sort.log"

  local so_tmp="${bam%.bam}.tmp."
  local so_bai="${so_bam%.bam}.bam.bai"
  # reference: https://github.com/samtools/samtools/issues/1196
  samtools sort -@"$nproc" -T"$so_tmp" \
    --write-index -o "$so_bam##idx##$so_bai" "$bam" >"$logfile" 2>&1 \
    || die "$0" "$LINENO" "Failed to run samtools sort on $bam"
  info "$0" ${LINENO} "Sort Bam file [DONE]"
}

function samtools_cat () {
  local bam_list="$1"
  local out="$2"

  info "$0" ${LINENO} "Concatenate individual BAM files provided in $bam_list"
  if [ ! -f "$bam_list" ];
  then
    die "$0" "$LINENO" "Failed to find the file with a list of BAMs"
  fi

  local outdir="${out%/*}"
  local logdir="$outdir/log"
  logdir=$( make_dir "$logdir" )
  local logfile="$logdir/samtools_cat.log"

  samtools cat -o "$out" -b "$bam_list" >"$logfile" 2>&1 \
    || die "$0" "$LINENO" "Failed to run samtools cat on $bam_list"
  info "$0" ${LINENO} "Concatenate individual BAM files [DONE]"
}

function picard_mdup () {
  local bam="$1"
  local mdup_bam="$2"
  local ram="$3"

  info "$0" ${LINENO} "Mark PCR duplicates in $bam"
  check_file_exists "$bam" \
    || die "$0" "$LINENO" "Failed to find $bam to mark duplicates"

  if [ -z "$ram" ];
  then
    ram=8   # in GB
  fi

  local metric="${mdup_bam%.bam}.mdup.metrics.txt"
  local outdir="${bam%/*}"
  local logdir="$outdir/log"
  logdir=$( make_dir "$logdir" )
  local logfile="$logdir/mdup.log"

  picard -XX:ParallelGCThreads=1 -Xmx"$ram"g MarkDuplicates \
    --INPUT "$bam" --OUTPUT "$mdup_bam" \
    --METRICS_FILE "$metric" --CREATE_INDEX >"$logfile" 2>&1 \
    || die "$0" "$LINENO" "Failed to run picard mdup on $bam"
  info "$0" ${LINENO} "Mark PCR duplicates in [DONE]"
}