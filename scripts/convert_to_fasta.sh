#!/bin/bash
# Convert text to FASTA file
# Assumes no header

# Default parameter
length=74

while getopts f:l:: option
do
  case "${option}"
  in
  f) fasta=${OPTARG};;
  l) length=${OPTARG};;
  exit;;
  esac
done

echo ">${fasta%.*}" > ${fasta%.*}.fa
fold -w $length ${fasta} >> ${fasta%.*}.fa

echo 'Done!'
