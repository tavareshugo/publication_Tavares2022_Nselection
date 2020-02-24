#!/bin/bash

#### Input arguments (in order) ####

INPUT="$1"   # input file
OUTPUT="$2"  # output file
REF="$3"     # reference genome file

# Make sure that it doesn't stop if the exit code is != 0
set +e

# Run picard
picard ValidateSamFile \
  INPUT="$INPUT" \
  OUTPUT="$OUTPUT" \
  REFERENCE_SEQUENCE="$REF" \
  MODE=SUMMARY

# Issue exit code
exitcode=$?
if [ $exitcode -eq 1 ]
then
    exit 1
else
    exit 0
fi