#!/bin/bash

../bin/Example \
  --mode Parse \
  --anc example.anc.gz \
  --mut example.mut.gz  \
  -o output

../bin/Example \
  --mode ParseData \
  --haps data/example.haps.gz \
  --sample data/example.sample.gz \
  --poplabels data/example.poplabels \
  -o output

../bin/Convert \
  --mode ConvertToTreeSequence \
  --anc example.anc.gz \
  --mut example.mut.gz  \
  -o output

../bin/Convert \
  --mode ConvertFromTreeSequence \
  --anc output.anc \
  --mut output.mut  \
  -i ./data/example.trees

