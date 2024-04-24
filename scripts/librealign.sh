#!/usr/bin/env bash


function run_razer () {
	fq="$1"
	hla_ref="$2"
  nproc="$3"

  if [ -z "$nproc" ];
  then
    nproc=1
  fi

	bam="${fq%.fastq.gz}.fished.bam"
  info "$0" "$LINENO" "Run razers3 realignment on Fastq: $fq"
	cmd="razers3 -i 95 -m 1 -dr 0 -tc $nproc -o $bam $hla_ref $fq"
	run_cmd "$cmd" "$LINENO" "Failed to fish reads using razerS3"
  info "$0" "$LINENO" "Run razers3 realignment on Fastq: $fq [DONE]"

}

function run_bam2fq () {
  bam="$1"

  if [ -z "$bam" ];
  then
    error "$0" "$LINENO" "Missing BAM file to run \'run_bam2fq\' function"
    exit 1
  fi
  check_file_exists "$bam"

  fq="${bam%.bam}.fastq.gz"

  info "$0" "$LINENO" "Extract fished reads from BAM: $bam"
  cmd="samtools bam2fq $bam > $fq 2>/dev/null"
  run_cmd "$cmd" "$LINENO" "Failed to extract fished reads from BAM file"
  info "$0" "$LINENO" "Extract fished reads from BAM: $bam [DONE]"
}

function run_seqkit_pair () {
  r1="$1"
  r2="$2"

  info "$0" "$LINENO" "Pair $r1 and $2"
  cmd="seqkit pair -1 $r1 -2 $r2 >/dev/null 2>&1"
  run_cmd "$cmd" "$LINENO" "Failed to pair fished reads"
  info "$0" "$LINENO" "Pair $r1 and $2 [DONE]"

}


function run_novoalign () {
  nix="$1"
  r1="$2"
  r2="$3"

  bam="${r1%.fastq.gz}.bam"
  done="${bam%.bam}.done"
  info "$0" "$LINENO" "Run Novoalign on ${r1##*/} ${r2##*/}"
  cmd="novoalign -d $nix -F STDFQ -R 0 -r All -o SAM -o FullNW -f $r1 $r2 2>/dev/null | \
    samtools view -bh -o $bam - && \
    touch $done"
  run_cmd "$cmd" "$LINENO" "Failed to realign reads to HLA using Novoalign"
  info "$0" "$LINENO" "Run Novoalign on ${r1##*/} ${r2##*/} [DONE]"
}

export -f run_razer
export -f run_bam2fq
export -f run_novoalign
