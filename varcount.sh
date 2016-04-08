#!/bin/bash

SCRIPT_DIR=$(readlink -f ${0%/*})
REF=$1
SAMPLE=$2

  qsub $SCRIPT_DIR/submit_varcount.sh $REF $SAMPLE
