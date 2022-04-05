# Ring-RPQ

Repository for the prototype source code of the paper Time- and Space-Efficient Regular Path Queries on Graphs. This repository aims at reproducing the experiments of the paper.

## Queries and graph

The queries are available in data/paths.tsv of this repository.

The data used are available here: [Wikidata (about 1000M triples)](http://compact-leapfrog.tk/files/wikidata-enumerated.dat.gz).

## Instructions for running the ring

To run our code, please install an extended version of the library SDSL. Go to this [this repository](https://github.com/darroyue/sdsl-lite) and follow the instructions.

After the extended version of SDSL is installed, clone this repository and follow these steps:

1. Compile the source code by first moving to the directory that contains the code, and then execute: 
```Bash
bash build.sh
```
This shall create two executable files: build-index and query-index.

2. Download the version of Wikidata used in the paper (note the original triples have been enumerated):
[Wikidata (about 1000M triples)](http://compact-leapfrog.tk/files/wikidata-enumerated.dat.gz) and uncompress it.

3. To build the index run:
```Bash
./build-index <path-to-wikidata-file> 
```
This will create the index on the same directory were the wikidata file is. Keep all these files in the same directory.

4. Move files data/wikidata-enumerated.dat.P and data/wikidata-enumerated.dat.SO to the directory were the index is stored. Please keep these file names, or change them acordingly, keeping the same prefix for all of them.   

5. To run queries, do as follows:
```Bash
./query-index <path-to-index-file> data/paths.tsv 
```

## Running other systems

You can find scripts and instructions on [how to run the benchmark for Blazegraph, Jena, Virtuoso on the following page](http://compact-leapfrog.tk/) in the section "Instructions for using SPARQL Engines". Note that we did not include RDF3x in these experients as it does not support RPQs. You can follow the same instructions, but rather use the RPQ queries from this repository. Please note that RPQs were run under `DISTINCT`/set semantics for comparability (in SPARQL, while queries like `:p*` are defined under set semantics, queries like `:p/:q`, `:p|:q`, etc., are rewritten to joins, unions, etc., that use bag semantics).

For the Datalog engine, you can find code to convert RPQs to Datalog and the optimisations [in the following repository](https://github.com/aidhog/rpqs-to-datalog).

