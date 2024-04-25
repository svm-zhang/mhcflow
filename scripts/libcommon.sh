#!/usr/bin/env bash

set -e

function error () {
	echo "[${1##*/}:$2] ERROR: $3" 1>&2
	#echo "[$(basename "$1"):$2:$3] ERROR: $4" 1>&2
}

function info () {
	echo "[${1##*/}:$2] INFO: $3"
}

function send_err_msg () {
	echo "$1"
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
			error "$0" ${LINENO} "the given file path does not exist: ${f}"
			exit 1
		fi
	done
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

export -f info error parse_path check_file_exists check_dir_exists make_dir run_cmd
