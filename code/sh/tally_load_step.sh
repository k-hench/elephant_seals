#!/usr/bin/env bash
cat ${1} | \
  SnpSift filter "${2}" | \
  SnpSift extractFields \
  /dev/stdin CHROM POS S | \
  grep -v "CHROM" | \
  wc -l