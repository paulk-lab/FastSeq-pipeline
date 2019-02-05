# FastSeq Pipeline
##### Simple tools to process sequencing data from FastSeq

## Requirements
- At least 2GB of RAM to be dedicated to Java heap, ideally a computer with 
8GB should suffice.
- Docker
- Python 3

## Instructions

If you haven't already, please install docker from here:
https://www.docker.com/get-started

If you have an Apple computer with Mojave you will have python 3.7 installed 
by default. Otherwise install python from here:
https://www.python.org/downloads/

This pipeline uses a simple python script to launch a prebuilt docker container
that contains all the needed sequence processing tools. These tools are run in
sequence according to a different python script.

To run the pipeline simply type in the command line:

`python3 fastseq_pipeline.py <data_directory> <csv_file>`

The data directory should contain all sequencing files, reference sequences
and adapter sequences. The sequencing results should be in zipped fastq, the
reference and adapters should be in fasta. The adapters are trimmed using 
trimmomatic, so refer to documentation for trimmomatic for how to format the
adapter file or mirror the included example.

The CSV file should contain the following columns (order doesn't matter):
- Sample
- Forward Read Path
- Reverse Read Path
- Adapter Path
- Reference Path

Each of the paths should be relative to the data directory, for example if the
adapters are in `data/adapters/adapter1.fa` then the csv file should contain
`adapters/adapter1.fa`.

Finally the script will generate a final stats file which contains a selection
of statistics from picard and bcftools. This file will be output to the 
`Output` folder under the name `final_stats.tsv.txt`.

An example is included, to run it simply run:

`python3 fast_seq_pipeline.py data example.csv`

in this directory.

Feel free to modify the code as needed, the code and tooling is licensed under
the Creative Commons with Attribution license.