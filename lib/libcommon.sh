#!/usr/bin/env bash

set -e

function error () {
	echo "[${1##*/}:$2] ERROR - $3" 1>&2
}

function info () {
	echo "[${1##*/}:$2] INFO - $3"
}

function die () {
	rc=$?
	# when there is argument passing in "$#"
	# error the message
	(( $# )) && error "$1" "$2" "$3"
	exit "$(( rc == 0 ? 1 : rc ))"
}

function check_empty_str () {
	local var="$1"
	if [ -z "$var" ]; then
		return 1
	fi
}

function get_abs_path () {
	case "$2" in
		d) check_dir_exists "$path" \
			|| die "$0" "$LINENO" \
				"Failed to generate abs path b/c non-existence of provided dir path $1" ;;
		f) check_file_exists "$1" \
			|| die "$0" "$LINENO" \
				"Failed to generate abs path b/c non-existence of provided file path $1" ;;
		*) die "$0" "$LINENO" "Unrecognized type passing in get_abs_path: $2";;
	esac
	local abs_path
	abs_path=$(cd -P -- "$(dirname -- "$1")" && printf '%s\n' "$(pwd -P)/$(basename -- "$1")")
	echo "$abs_path"
}

function check_file_exists () {
	local paths=("$@")
	for f in "${paths[@]}"; do
		if [ ! -f "$f" ]; then
			#error "${FUNCNAME[0]}" ${LINENO} "Failed to find the file path at the given path: ${f}"
			return 1
		fi
	done
	return 0
}

function check_dir_exists () {
	local paths=("$@")
	for d in "${paths[@]}"; do
		if [ ! -d "$d" ]; then
			#error "$0" ${LINENO} "Failed to find directory for the given path: ${d}"
			#exit 1
			return 1
		fi
	done
	return 0
}

function make_dir () {
	local path="$1"
	if [ ! -d "$path" ];
	then
		mkdir -p "$path"	
	fi
	abs_path=$( get_abs_path "$path" "d")
	echo "$abs_path"

	#local dir_path
	#dir_path=$(cd -P -- "$(dirname -- "$1")" && printf '%s\n' "$(pwd -P)/$(basename -- "$1")")
	#echo "$dir_path"
}

export -f info error die
export -f get_abs_path check_file_exists check_dir_exists make_dir
