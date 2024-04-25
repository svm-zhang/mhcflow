#!/usr/bin/env bash


function find_files_on_pattern () {
  local wkdir="$1"
  local pattern="$2"
  local delimiter="$3"

  if [ -z "$delimiter" ];
  then
    delimiter="\n"
  fi

  if [ -z "$pattern" ];
  then
    error "${FUNCNAME[0]}" "$LINENO" "pattern variable cannot be empty to do a search"
    exit 1
  fi

  if [ -z "$wkdir" ];
  then
    error "${FUNCNAME[0]}" "$LINENO" "wkdir variable cannot be empty to find files within"
    exit 1
  fi
  check_dir_exists "$wkdir"

  local array=()
  while IFS= read -r -d $'\0'; do
    array+=("$REPLY")
  done < <(find "$wkdir" -type f -name "$pattern" -print0)

  if [ "${#array[@]}" -eq 0 ];
  then
    error "${FUNCNAME[0]}" "$LINENO" "Found no files in $wkdir given pattern \"$pattern\""
    exit 1
  fi

  printf "%s$delimiter" "${array[@]}"
}

function run_seqkit_split2 () {
  local r1="$1"
  local r2="$2"
  local wkdir="$3"
  local nproc="$4"

  if [ -z "$r1" ] || [ -z "$r2" ];
  then
    error "${FUNCNAME[0]}" "$LINENO" "r1 and r2 variables cannot be empty to run seqkit split2 command"
    exit 1
  fi

  if [ -z "$wkdir" ];
  then
    error "${FUNCNAME[0]}" "$LINENO" "wkdir variable cannot be empty to run seqkit split2 command"
    exit 1
  fi

  local logdir="${wkdir%/*}/log"
  logdir=$( make_dir "$logdir" )
  local logfile="$logdir/seqkit_split2.log"
  local donefile="$logdir/seqkit_split2.done"

  if [ -f "$donefile" ]; then
    info "${FUNCNAME[0]}" "$LINENO" "Found previous seqkit split2 done file. Skip"
    return 0
  fi

  info "${FUNCNAME[0]}" "$LINENO" "Split read files into $nproc parts"
  local cmdArgs=( "-p" "$nproc" "-e" ".gz" "-j" "$nproc" )
  local cmd=(
    "seqkit"
    "split2"
    "-1"
    "$r1"
    "-2"
    "$r2"
    "-O"
    "$wkdir"
  )
  cmd+=( "${cmdArgs[@]}" )

  if ! "${cmd[@]}" >"$logfile" 2>&1 ;
  then
    error "${FUNCNAME[0]}" "$LINENO" "Failed to run seqkit split2 command"
    error "${FUNCNAME[0]}" "$LINENO" "Please check command: ${cmd[*]}"
    return 255
  fi

  info "${FUNCNAME[0]}" "$LINENO" "Check if expected number of Fastq files generated"
  local n_expect_files
  local n_found_files
  n_expect_files=$(( "$nproc" * 2 ))
  local out_files
  # for some reason, I cannot figure out using -regex and -regextype
  # to make find command return the correct files
  # thats why I borrow the power of grep...
  out_files=$( find_files_on_pattern "$wkdir" "*part_[0-9]*.fastq.gz" | \
    grep -E "*part_[0-9]+.fastq.gz"
  )
  n_found_files=$( echo "$out_files" | wc -l )
  if [ "$n_found_files" -ne "$n_expect_files" ];
  then
    error "${FUNCNAME[0]}" "$LINENO" "Seqkit split2 expected to generate $n_expect_files"
    error "${FUNCNAME[0]}" "$LINENO" "Seqkit split2 generated $n_found_files"
    exit 1
  fi

  # check if there are even number of files generated
  check_even=$(( "$n_found_files" % 2 ))
  if [ "$check_even" -ne 0 ];
  then
    error "${FUNCNAME[0]}" "$LINENO" "Seqkit split2 on paired-end data generates odd number of files"
    error "${FUNCNAME[0]}" "$LINENO" "Please check command: ${cmd[*]}"
    exit 1
  fi

  touch "$donefile"
  
}

function run_razers3_batch () {
  local wkdir="$1"
  local hlaref="$2"
  local pattern="$3"
  local nproc="$4"
  local nproc_per_job="$5"

  if [ -z "$wkdir" ];
  then
    error "${FUNCNAME[0]}" "$LINENO" "wkdir variable cannot be empty to run razers3 in batch"
    exit 1
  fi

  if [ -z "$hlaref" ];
  then
    error "${FUNCNAME[0]}" "$LINENO" "razers3 reaglinment requires hlaref variable to be not empty"
    exit 1
  fi

  # make sure nproc and nproc_per_job have a default value of 1
  if [ -z "$nproc" ];
  then
    nproc=1
  fi

  if [ -z "$nproc_per_job" ];
  then
    nproc_per_job=1
  fi

  local logdir="${wkdir%/*}/log"
  logdir=$( make_dir "$logdir" )
  local donefile="$logdir/razers3_fish.done"

  if [ -f "$donefile" ]; then
    info "${FUNCNAME[0]}" "$LINENO" "Found previous razers3 fished alignment done file. Skip"
    return 0
  fi

  local fq_files
  fq_files=$( find_files_on_pattern "$wkdir" "$pattern" | \
    grep -E "*part_[0-9]+.fastq.gz"
  )

  if [ -z "$fq_files" ];
  then
    error "$0" "$LINENO" "Failed to find any Fastq files within ${wkdir} on pattern \"$pattern\""
    exit 1
  fi

  # now run razers3 in batch
  echo "${fq_files[@]}" | \
    xargs -P"$nproc" -I{} bash -c 'run_razer "$@" || exit 255 ' "$0" "{}" "$hlaref" "$nproc_per_job"

  # when one fastq file is empty, razers3 will issue error
  # so the above xargs run will terminate when that happens
  # therefore, if there are 4 files going in, I should expect 4 files going out
  # when razers3 finds no reads alignemnt, it still generates a BAM with only header
  local n_expect_files
  local n_found_files=0
  n_expect_files=$( echo "${fq_files[@]}" | wc -l )
  local out_bams
  out_bams=$( find_files_on_pattern "$wkdir" "*part_[0-9]*.razers3.bam" | \
    grep -E "*part_[0-9]+.razers3.bam"
  )
  n_found_files=$( echo "${out_bams[@]}" | wc -l ) 
  if [ "$n_found_files" -ne "$n_expect_files" ];
  then
    error "${FUNCNAME[0]}" "$LINENO" "razers3 expected to generate $n_expect_files"
    error "${FUNCNAME[0]}" "$LINENO" "razers3 generated $n_found_files"
    exit 1
  fi

  touch "$donefile"
}

function run_razer () {
	fq="$1"
	hlaref="$2"
  nproc="$3"

  if [ -z "$fq" ];
  then
    error "${FUNCNAME[0]}" "$LINENO" "razers3 requires fq variable to be non-empty"
    exit 255
  fi
  check_file_exists "$hlaref"

  if [ -z "$hlaref" ];
  then
    error "${FUNCNAME[0]}" "$LINENO" "HLA reference variable cannot be empty to run razerS3"
    exit 255
  fi
  check_file_exists "$hlaref"

  if [ -z "$nproc" ];
  then
    nproc=1
  fi

  local curdir="${fq%/*}"
  local logdir="${curdir%/*}/log"
  logdir=$( make_dir "$logdir" )
  local logprefix="${fq##*/}"
  local logfile="$logdir/${logprefix%.fastq.gz}.razers3.log"
  local donefile="$logdir/${logprefix%.fastq.gz}.razers3.done"
  local failfile="$logdir/${logprefix%.fastq.gz}.razers3.fail"

  local bam
	bam="${fq%.fastq.gz}.razers3.bam"
  if [ -f "$donefile" ] && [ -f "$bam" ]; then
    info "${FUNCNAME[0]}" "$LINENO" "Found previous razers3 fished alignment done file. Skip"
    return 0
  fi

  info "${FUNCNAME[0]}" "$LINENO" "Run razers3 realignment on Fastq: ${fq##*/}"
  local cmdArgs=("-i" "95" "-m" "1" "-dr" "0" "-tc" "$nproc")
  local program="razers3"
  local cmd=(
    "$program"
    "${cmdArgs[@]}"
    "-o"
    "$bam"
    "$hlaref"
    "$fq"
  )

  if ! "${cmd[@]}" >"$logfile" 2>&1 ;
  then
    error "${FUNCNAME[0]}" "$LINENO" "Failed to run razers3 command"
    error "${FUNCNAME[0]}" "$LINENO" "Please check command: ${cmd[*]}"
    touch "$failfile"
    exit 255
  fi

  touch "$donefile"

}

function run_bam2fq_batch () {

  local wkdir="$1" 
  local pattern="$2"
  local nproc="$3"

  if [ -z "$wkdir" ];
  then
    error "${FUNCNAME[0]}" "$LINENO" "wkdir varialbe cannot be empty to run run_bam2fq_batch function"
    exit 1
  fi
  check_dir_exists "$wkdir" 

  local logdir="${wkdir%/*}/log"
  logdir=$( make_dir "$logdir" )
  local donefile="$logdir/bam2fq.done"

  if [ -f "$donefile" ]; then
    info "${FUNCNAME[0]}" "$LINENO" "Found previous bam2fq done file. Skip"
    return 0
  fi

  local bam_files
  bam_files=$( find_files_on_pattern "$wkdir" "$pattern" )
  if [ -z "$bam_files" ];
  then
    error "${FUNCNAME[0]}" "$LINENO" "Failed to find any BAM files to run bam2fq function"
    exit 1
  fi

  echo "${bam_files[@]}" | \
    xargs -P"$nproc" -I{} bash -c 'run_bam2fq "$@" || exit 255' "_" "{}"

  # I dont need to check if number of inputs equals to the number of outputs
  # becaues some bam2fq might not return any reads from bam file
  touch "$donefile"

}

# this function only dumps all reads in one fq
function run_bam2fq () {
  local bam="$1"

  if [ -z "$bam" ];
  then
    error "${FUNCNAME[0]}" "$LINENO" "Missing BAM file to run run_bam2fq function"
    exit 255
  fi
  check_file_exists "$bam"

  local curdir="${bam%/*}"
  local logdir="${curdir%/*}/log"
  logdir=$( make_dir "$logdir" )
  local logprefix="${bam##*/}"
  local logfile="$logdir/${logprefix%.bam}.bam2fq.log"
  local donefile="$logdir/${logprefix%.bam}.bam2fq.done"
  local failfile="$logdir/${logprefix%.bam}.bam2fq.fail"

  local fq="${bam%.bam}.fastq.gz"
  if [ -f "$donefile" ] && [ -f "$fq" ]; then
    info "${FUNCNAME[0]}" "$LINENO" "Found previous bam2fq done file. Skip"
    return 0
  fi

  info "${FUNCNAME[0]}" "$LINENO" "Extract fished reads from BAM: ${bam##*/}"
  # note: samtools bam2fq bam > fq get non-suppl/secondary
  # in this case, I was losing data
  # Using the command below gets all the reads
  local cmd=("samtools" "bam2fq" "-0" "$fq" "$bam")
  if ! "${cmd[@]}"  >"$logfile" 2>&1;
  then
    error "${FUNCNAME[0]}" "$LINENO" "Failed to run samtools bam2fq command"
    error "${FUNCNAME[0]}" "$LINENO" "Please check command: ${cmd[*]}"
    touch "$failfile"
    exit 255
  fi

  local fq_size
  # first field returns the size of the compressed file
  fq_size=$(zcat "$fq" | head | wc -l)
  if [ "$fq_size" -eq 0 ];
  then
    info "${FUNCNAME[0]}" "$LINENO" "samtools bam2fq returned empty Fastq file. Deleting"
    # -f in case they do not exists...
    rm -f "$fq"
  fi

  touch "$donefile"

}

function run_seqkit_pair_batch () {
  local wkdir="$1"
  local pattern="$2"
  # reason to have gunzip as part of the argument is because
  # novoalign (next step) does not take gzipped fastq files
  # kinda lazy to write separate xargs gunzip on all paired reads
  local gunzip="$3"
  local nproc="$4"

  if [ -z "$wkdir" ];
  then
    error "${FUNCNAME[0]}" "$LINENO" "wkdir varialbe cannot be empty to search for files given pattern"
    exit 1
  fi
  check_dir_exists "$wkdir" 

  if [ -z "$nproc" ];
  then
    nproc=1
  fi

  local logdir="${wkdir%/*}/log"
  logdir=$( make_dir "$logdir" )
  local donefile="$logdir/seqkit_pair.done"

  if [ -f "$donefile" ]; then
    info "${FUNCNAME[0]}" "$LINENO" "Found previous seqkit pair done file. Skip"
    return 0
  fi

  paired_input_str=$( make_paired_r1_r2_str "$wkdir" "$pattern" )
  # r1 and r2 files passed in as one argument, separated by comma
  # thats why in run_seqkit_pair funciton, we first split string to get
  # r1 and r2 respectively
  echo "$paired_input_str" | \
    xargs -P"$nproc" -I{} bash -c 'run_seqkit_pair "$@" || exit 255' "_" "{}" "$gunzip"

  # if no paired reads found, the empty file will be deleted in the subprocess
  # therefore, I dont need to check file counts
  touch "${donefile}"
}

function run_seqkit_pair () {

  local reads="$1"
  local gunzip="$2"
  IFS="," read -r r1 r2 <<< "$reads"

  # check if two input read files are compressed
  if ! gzip -t "$r1" >/dev/null 2>&1;
  then
    error "${FUNCNAME[0]}" "$LINENO" "Input R1 Fastq file is not gzipped: ${r1}"
    exit 255
  fi
  if ! gzip -t "$r2" >/dev/null 2>&1;
  then
    error "${FUNCNAME[0]}" "$LINENO" "Input R2 Fastq file is not gzipped: ${r2}"
    exit 255
  fi

  local curdir="${r1%/*}"
  local logdir="${curdir%/*}/log"
  logdir=$( make_dir "$logdir" )
  local logprefix="${r1##*/}"
  local logfile="$logdir/${logprefix%.fastq.gz}.paired.log"
  local donefile="$logdir/${logprefix%.fastq.gz}.paired.done"
  local failfile="$logdir/${logprefix%.fastq.gz}.paired.fail"

  local paired_r1="${r1%.fastq.gz}.paired.fastq.gz"
  local paired_r2="${r2%.fastq.gz}.paired.fastq.gz"

  # when guznip is true, I should check for suffix with "fastq", not
  # "fastq.gz"
  if [ "$gunzip" = true ];
  then
    paired_r1_unzip="${paired_r1%.fastq.gz}.fastq"
    paired_r2_unzip="${paired_r1%.fastq.gz}.fastq"
    if [ -f "$donefile" ] && [ -f "$paired_r1_unzip" ] && [ -f "$paired_r2_unzip" ];
    then
      info "${FUNCNAME[0]}" "$LINENO" "Found previous seqkit pair done file. Skip"
      return 0
    fi
  else
    if [ -f "$donefile" ] && [ -f "$paired_r1" ] && [ -f "$paired_r2" ];
    then
      info "${FUNCNAME[0]}" "$LINENO" "Found previous seqkit pair done file. Skip"
      return 0
    fi
  fi

  info "${FUNCNAME[0]}" "$LINENO" "Pair fished reads from ${r1##*/} ${r2##*/}"
  local cmd
  cmd=("seqkit" "pair" "-1" "$r1" "-2" "$r2")
  if ! "${cmd[@]}"  >"$logfile" 2>&1;
  then
    error "${FUNCNAME[0]}" "$LINENO" "Failed to run samtools bam2fq command"
    error "${FUNCNAME[0]}" "$LINENO" "Please check command: ${cmd[*]}"
    rm -f "$paired_r1" "$paired_r2"
    touch "$failfile"
    exit 255
  fi

  # check size of the paired read files
  local r1_size
  local r2_size
  r1_size=$(zcat "$paired_r1" | head | wc -l)
  r2_size=$(zcat "$paired_r2" | head | wc -l)

  if [ "$r1_size" -eq 0 ] || [ "$r2_size" -eq 0 ];
  then
    info "${FUNCNAME[0]}" "$LINENO" "Empty paired R1 or R2 Fastq file generated. Deleting"
    # -f in case they do not exists...
    rm -f "$paired_r1" "$paired_r2"
  else
    if [ "$gunzip" = true ];
    then
      gunzip -f "$paired_r1" "$paired_r2"
    fi
  fi

  touch "$donefile"
}

function make_paired_r1_r2_str () {
  local wkdir="$1"
  local pattern="$2"

  # Watch-out alert 
  # need more tests on different input r1 r2 files to make sure
  # this works
  local r1_pattern="*1$pattern"
  local r2_pattern="*2$pattern"
  local r1_files r2_files
  r1_files=$( find_files_on_pattern "$wkdir" "$r1_pattern" )
  r2_files=$( find_files_on_pattern "$wkdir" "$r2_pattern" )
  if [ -z "$r1_files" ] || [ -z "$r2_files" ];
  then
    error "${FUNCNAME[0]}" "$LINENO" "Failed to find any read 1 and 2 files to make a string"
    exit 255
  fi

  paste -d ',' <(echo "${r1_files[@]}" | sort) <(echo "${r2_files[@]}" | sort)

}

# FIXME: now if you look at all batch function, too much duplicated code
function run_novoalign_batch () {
  local wkdir="$1"
  local pattern="$2"
  local hla_ref_nix="$3"
  local nproc="$4"

  if [ -z "$wkdir" ];
  then
    error "${FUNCNAME[0]}" "$LINENO" "wkdir varialbe cannot be empty to search for files given pattern"
    exit 1
  fi
  check_dir_exists "$wkdir" 

  if [ -z "$nproc" ];
  then
    nproc=1
  fi

  local logdir="${wkdir%/*}/log"
  logdir=$( make_dir "$logdir" )
  local donefile="$logdir/novoalign.done"

  if [ -f "$donefile" ]; then
    info "${FUNCNAME[0]}" "$LINENO" "Found previous novoalign alignment done file. Skip"
    return 0
  fi

  paired_input_str=$(make_paired_r1_r2_str "$wkdir" "$pattern")

  echo "$paired_input_str" | \
    xargs -P"$nproc" -I{} bash -c 'run_novoalign "$@" || exit 255' "_" "{}" "$hla_ref_nix"

  touch "$donefile"
}

function run_novoalign () {
  local reads="$1"
  local nix="$2"

  IFS="," read -r r1 r2 <<< "$reads"

  local curdir="${r1%/*}"
  local logdir="${curdir%/*}/log"
  logdir=$( make_dir "$logdir" )
  local logprefix="${r1##*/}"
  local logfile="$logdir/${logprefix%.fastq}.novoalign.log"
  local donefile="$logdir/${logprefix%.fastq}.novoalign.done"
  local failfile="$logdir/${logprefix%.fastq}.novoalign.fail"

  local bam
  bam="${r1%.fastq}.novoalign.bam"
  if [ -f "$donefile" ] && [ -f "$bam" ]; then
    info "${FUNCNAME[0]}" "$LINENO" "Found previous novoalign alignment done file. Skip"
    return 0
  fi

  info "${FUNCNAME[0]}" "$LINENO" "Run Novoalign on ${r1##*/} ${r2##*/}"
  local program="novoalign"
  local cmdArgs=("-F" "STDFQ" "-R" "0" "-r" "All" "-o" "SAM" "-o" "FullNW")
  local novocmd
  novocmd=(
    "$program"  
    "-d"
    "$nix"
    "${cmdArgs[@]}"
    "-f"
    "$r1"
    "$r2"
  )
  local samtools_cmd=(
    "samtools"
    "view"
    "-bh"
    "-o"
    "$bam"
    "-"
  )

  if ! "${novocmd[@]}" 2>"$logfile" | "${samtools_cmd[@]}" >/dev/null 2>&1;
  then
    error "${FUNCNAME[0]}" "$LINENO" "Failed to run novoalign command"
    error "${FUNCNAME[0]}" "$LINENO" "Please check command: ${cmd[*]}"
    touch "$failfile"
    exit 255
  fi

  touch "$donefile"
}

function run_samtools_cat () {
  local wkdir="$1"
  local pattern="$2"
  local out="$3"

  if [ -z "$wkdir" ];
  then
    error "${FUNCNAME[0]}" "$LINENO" "wkdir variable cannot be empty to run razers3 in batch"
    exit 1
  fi

  local logdir="${wkdir%/*}/log"
  logdir=$( make_dir "$logdir" )
  local donefile="$logdir/samtools_cat.done"
  local logfile="$logdir/samtools_cat.log"

  if [ -f "$donefile" ] && [ -f "$out" ];
  then
    info "${FUNCNAME[0]}" "$LINENO" "Found previous samtools cat done file. Skip"
    return 0
  fi

  bam_input_str=$( find_files_on_pattern "$wkdir" "$pattern")
  bam_list_file="$wkdir/bams.list.txt"
  echo "${bam_input_str[@]}" > "$bam_list_file"
  cmd=(
    "samtools"
    "cat"
    "-o"
    "$out"
    "-b"
    "$bam_list_file"
  )
  if ! "${cmd[@]}" >"$logfile" 2>&1 ;
  then
    error "$0" "$LINENO" "Failed to concatenate novoalign-aligned BAM files"
    exit 1
  fi

  touch "$donefile"

}

export -f run_razer
export -f run_bam2fq
export -f run_seqkit_pair
export -f run_novoalign
