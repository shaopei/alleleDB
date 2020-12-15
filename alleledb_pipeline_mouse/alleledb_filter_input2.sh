#!/bin/env sh
# This script should take whatever input is provided and stream a fasta format suitable for bowtie

PL=$1
INFILE=$2

zcat ${INFILE} | python2 ${PL}/fastq2result.py - - |  python2 ${PL}/filter_query.py - - | python2 ${PL}/alleledb_ConvertTags.py
