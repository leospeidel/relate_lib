#!/bin/bash

../bin/Convert --mode ConvertToTreeSequence --anc example.anc.gz --mut example.mut.gz -o example

python3 ./load.py
#../bin/Convert --mode ConvertFromTreeSequence --anc example.anc --mut example.mut -i example.trees
