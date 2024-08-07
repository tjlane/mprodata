#!/bin/zsh

PDBFILE=$1
MTZFILE=$2

phenix.ready_set $PDBFILE
phenix.ensemble_refinement ${PDBFILE::-4}.updated.pdb ${MTZFILE} ptls=0.7
