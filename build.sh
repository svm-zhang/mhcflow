#!/bin/bash

set -e

mkdir -p "$PREFIX/bin"
mkdir -p "$PREFIX/lib"

# copy scripts to bin folder
cp "$RECIPE_DIR/scripts/mhcflow.sh" "$PREFIX/bin/mhcflow"
cp "$RECIPE_DIR/scripts/fisher.sh" "$PREFIX/bin/fisher"
cp "$RECIPE_DIR/scripts/mhcrealigner.sh" "$PREFIX/bin/realigner"
cp "$RECIPE_DIR/scripts/extract_sample_hlaref.sh" "$PREFIX/bin/extractor"
# cp "$RECIPE_DIR/pyhlatyper/pyhlatyper.py" "$PREFIX/bin/pyhlatyper"

# make it executable
chmod +x "$PREFIX"/bin/*

# copy perl module to lib folder
cp "$RECIPE_DIR/lib/libcommon.sh" "$PREFIX/lib"
cp "$RECIPE_DIR/lib/librealign.sh" "$PREFIX/lib"
cp "$RECIPE_DIR/lib/libbamer.sh" "$PREFIX/lib"

set -x
"$PREFIX/bin/pip3" install "mhctyper>=0.1.4" -vvv
set +x
