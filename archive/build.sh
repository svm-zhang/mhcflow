#!/bin/bash

set -e

mkdir -p "$PREFIX/bin"
mkdir -p "$PREFIX/lib"

# copy scripts to bin folder
cp "$RECIPE_DIR/scripts/hlareforged.sh" "$PREFIX/bin/hlareforged"
cp "$RECIPE_DIR/scripts/razer_realigner.sh" "$PREFIX/bin/razer_realigner"
cp "$RECIPE_DIR/scripts/hlatyper.sh" "$PREFIX/bin/hlatyper"
cp "$RECIPE_DIR/scripts/hlafinalizer.sh" "$PREFIX/bin/hlafinalizer"
cp "$RECIPE_DIR/scripts/hlapolysolver.sh" "$PREFIX/bin/hlapolysolver"
cp "$RECIPE_DIR/scripts/polysolver_realigner.sh" "$PREFIX/bin/polysolver_realigner"
cp "$RECIPE_DIR/scripts/first_allele_calculations.pl" "$PREFIX/bin/first_allele_calculations"
cp "$RECIPE_DIR/scripts/second_allele_calculations.pl" "$PREFIX/bin/second_allele_calculations"

# make it executable
chmod +x "$PREFIX"/bin/*

# copy perl module to lib folder
cp "$RECIPE_DIR/lib/PolysolverUtils.pm" "$PREFIX/lib/perl5/site_perl/5.22.0"
cp "$RECIPE_DIR/lib/libcommon.sh" "$PREFIX/lib"
cp "$RECIPE_DIR/lib/librealign.sh" "$PREFIX/lib"
