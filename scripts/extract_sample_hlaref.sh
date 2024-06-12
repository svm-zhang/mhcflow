#!/usr/bin/env bash

set -e
set -o pipefail

SRC_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
LIB_DIR="${SRC_DIR%/*}/lib"
source "$LIB_DIR/libcommon.sh"

function usage () {
	local program
	program=$(basename "$0")
	cat <<EO
usage: $program [options]
Options:
EO
	cat <<EO | column -s\& -t
	--sample    & Specify the sample name (Required)
	--typeres    & Specify the path to the HLA typing result (Required)
	--hla_ref    & Specify the HLA reference sequences in Fasta (Required)
	--out    & Specify the path to the final analysis-ready bam (Required)
EO
}

function extractor_main() {
	local sample=
	local typing_res=
	local hla_ref=
	local out=

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
		--sample) shift; sample="$1";;
		--typeres) shift; typing_res="$1";;
		--hla_ref) shift; hla_ref="$1";;
		--out) shift; out="$1";;
		--) shift; break;;
		*)
			echo "Invalid option: $1" 1>&2
			usage
			exit 1
			;;
		esac
		shift
	done

	check_empty_str "$sample" \
		|| die "$0" "$LINENO" \
			"extractor requires --sample to operate"
	check_empty_str "$typing_res" \
		|| die "$0" "$LINENO" "extractor requires typing result (--typing_res)"
	check_empty_str "$hla_ref" \
		|| die "$0" "$LINENO" "extractor requires a HLA ref file (--hla_ref)"

	check_file_exists "$typing_res" \
		|| die "$0" "$LINENO" "extractor failed to find typing result: $typing_res"
	check_file_exists "$hla_ref" \
		|| die "$0" "$LINENO" "extractor failed to find HLA ref: $hla_ref"

	typing_res=$( get_abs_path "$typing_res" "f" )
	hla_ref=$( get_abs_path "$hla_ref" "f" )

  local outdir="${out%/*}"
  outdir=$( make_dir "$outdir" )
  local logdir="$outdir/log"
  logdir=$( make_dir "$logdir" )
  donefile="$logdir/extract_sample_hlaref.done"
  if [ -f "$donefile" ]; then
    info "$0" "$LINENO" "Found previously HLA reference extraction result"
    info "$0" "$LINENO" "Move on"
    return
  fi

	# et sample-level HLA reference based on hlatyping result
	info "$0" "$LINENO" "Get sample-level HLA reference sequence"
	local hlatypelist_file="${outdir}/.${sample}.hlatypelist.txt"
	awk '(NR > 1) {print $1}' "$typing_res" > "$hlatypelist_file"
  # get number of unique alleles in the typing result
  # each locus can be homozygous so thats why we need to uniq
	local nexpect
	nexpect=$( sort "$hlatypelist_file" | uniq | wc -l )
	#local sample_hla_ref="$outdir/$sample.hla.fasta"
	local seqkit_grep_log="$logdir/$sample.seqkit_grep.log"
	seqkit grep -f "$hlatypelist_file" \
		-o "$out" "$hla_ref" >"$seqkit_grep_log" 2>&1 \
		|| die "$0" "$LINENO" "Failed to grep sequences of typed HLA alleles"

	# check if Fasta is empty
	if [ ! -s "$out" ]; then
		die "$0" "$LINENO" \
      "Failed to extract any sequences in ${out##*/}"
	fi

	# check if the expected # alleles extracted
	local nseq
	nseq=$( grep -c ">" "$out" )
	if [ "$nseq" -ne "$nexpect" ]; then
    # if not, remove the fasta before die
		rm -f "$out"
		die "$0" "$LINENO" "Expect to extract $nexpect alleles, but got $nseq"
	fi
	info "$0" "$LINENO" "Get sample-level HLA reference sequence [DONE]"

	# novoindex the extracted hlaref as new ref
	info "$0" "$LINENO" "Index HLA reference using novoindex"
	local hlaref_nix="${out%.fasta}.nix"
	local novoidx_log="$logdir/$sample.novoindex.log"
	novoindex "$hlaref_nix" "$out" >"$novoidx_log" 2>&1 \
		|| die "$0" "$LINENO" "Failed to index the sample-level HLA sequences"
	info "$0" "$LINENO" "Index HLA reference using novoindex [DONE]"

  touch "$donefile"
}

extractor_main "$@"

exit 0
