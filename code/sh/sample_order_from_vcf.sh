#!/usr/bin/env bash
bcftools view -O v ${1} | grep -v "^##" | head -n 1 | cut -f 10- | sed 's/\t/\n/g'