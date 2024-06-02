#!/usr/bin/env bash

set -e

function error () {
	echo "[${1##*/}:$2] ERROR: $3" 1>&2
}

function info () {
	echo "[${1##*/}:$2] INFO: $3"
}

function die () {
	rc=$?
	# when there is argument passing in "$#"
	# error the message
	(( $# )) && error "$1" "$2" "$3"
	exit "$(( rc == 0 ? 1 : rc ))"
}

function parse_path () {
	check_file_exists "$1"
	local abs_path
	abs_path=$(cd -P -- "$(dirname -- "$1")" && printf '%s\n' "$(pwd -P)/$(basename -- "$1")")
	echo "$abs_path"
}

function check_file_exists () {
	local paths=("$@")
	for f in "${paths[@]}"; do
		if [ ! -f "$f" ]; then
			error "${FUNCNAME[0]}" ${LINENO} "Failed to find the file path at the given path: ${f}"
			return 1
		fi
	done
	return 0
}

function check_dir_exists () {
	local paths=("$@")
	for d in "${paths[@]}"; do
		if [ ! -d "$d" ]; then
			error "$0" ${LINENO} "Failed to find directory for the given path: ${d}"
			exit 1
		fi
	done
}

function make_dir () {

	if [ ! -d "$1" ];
	then
		mkdir -p "$1"	
	fi

	local dir_path
	dir_path=$(cd -P -- "$(dirname -- "$1")" && printf '%s\n' "$(pwd -P)/$(basename -- "$1")")
	echo "$dir_path"
}

function run_cmd () {

	cmd="$1"
	where="$2"
	err_msg="$3"

	if ! eval "$cmd";
	then
		error "$0" "$where" "$err_msg"
		error "$0" "$where" "Please check command: \"$cmd\""
		exit 255 
	fi

}

export -f info error die parse_path check_file_exists check_dir_exists make_dir run_cmd
