
# relate\_lib 
### C++ library containing functions for parsing/dumping Relate output files 
Leo Speidel <br/> Department of Statistics, University of Oxford, Oxford, UK

## Requirements

C++11 and cmake version 3.1.

## Installation

Clone/fork from this github repository to a local directory.
Then execute

```` bash
cd PATH_TO_RELATE_LIB/
mkdir build
cd build
cmake ..
make
````

## Documentation

### Conversion from tskit format
````bash
relate_lib/bin/Convert --mode ConvertFromTreeSequence \
              --anc example.anc.gz \
              --mut example.mut.gz \
              -i example
````

### Conversion to tskit format
````bash
relate_lib/bin/Convert --mode ConvertToTreeSequence \
              --anc example.anc.gz \
              --mut example.mut.gz \
              -o example
````
Thanks to <b>Nathaniel S. Pope</b>, you can now also specify an argument that <b>compresses</b> these Relate-converted tree sequences by assigning the same age to nodes with identical descendant sets across adjacent trees.
````bash
relate_lib/bin/Convert --mode ConvertToTreeSequence \
              --compress \
              --anc example.anc.gz \
              --mut example.mut.gz \
              -o example
````

### C++ library 
Please see ./include/example/Example.cpp for example code and ./example/example.sh for examples and how to run them.

