#!/usr/bin/env bash

set -e
set -o pipefail

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
	novoalign -d "$nix" -F STDFQ -R 0 -r All -o SAM "$rg_str" -o FullNW \
		-f "$r1" "$r2" 2>"$logfile" \
		| samtools view -bh -o "$bam" \
		|| die "$0" "$LINENO" "Failed to novoalign on $r1 and $r2"
}

export -f run_novoalign
