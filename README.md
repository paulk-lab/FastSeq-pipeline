# FastSeq Pipeline
##### Simple tools to process sequencing data from FastSeq

## Requirements
- At least 2GB of RAM to be dedicated to Java heap, ideally a computer with 
8GB should suffice.
- Docker
- Python 3
- Pandas

## Instructions

If you haven't already, please install Docker from here:
https://www.docker.com/get-started

If you have an Apple computer with Mojave you will have Python 3.7 installed 
by default. Otherwise install Python from here:
https://www.python.org/downloads/

This pipeline uses a simple Python script to launch a prebuilt Docker container
that contains all the needed sequence processing tools and will run the 
analysis **exactly** as was run for the paper. These tools are run in sequence 
according to a different Python script, `process_seq.py`. 

For the most part, if you wish to understand the exact parameters used to 
analyze the sequences, please look at the top of `process_seq.py`, you will 
find a section where configuration variables are set. You can choose to run 
the pipeline using the docker container *or* you will need to modify the paths 
(also listed at the top of the `process_seq.py` file) to point to a set of 
locally installed tools. **However** should you modify the paths and use your
own set of tools we cannot guarantee that the pipeline will work as expected.

There may be some other configurations of interest, for this look deeper in the
`process_seq.py` script.

To run the pipeline as originally configured, simply type the following in the 
command line:

`python3 run_fastseq_pipeline.py process_seq.py <data_directory> <csv_file>`

The data directory should contain all sequencing files, reference sequences
and adapter sequences. The sequencing results should be in zipped FASTQ, the
reference and adapters should be in FASTA. The adapters are trimmed using 
trimmomatic, so refer to documentation for trimmomatic for how to format the
adapter file or mirror the included example.

The CSV file should contain the following columns (order doesn't matter):
- Sample
- Forward Read Path
- Reverse Read Path
- Adapter Path
- Reference Path

Each of the paths should be relative to the data directory. For example, if the
adapters are in `data/adapters/adapter1.fa` then the csv file should contain
`adapters/adapter1.fa`.

Finally, the script will generate a final stats file which contains a selection
of statistics from picard and bcftools. This file will be output to the 
`Output` folder under the name `final_stats.csv`.

An example is included. To run it, simply execute the following command:

`python3 run_fastseq_pipeline.py process_seq.py data example.csv`

in this directory.

Most of the output the script makes are put into a single `Output` directory
that is newly generated. This folder must be deleted every run of the script.
However, some files generated by the script (mostly non-overwriting 
modifications of the original sequence data) are also output in the same 
location as the input sequencing data.

When run on the demo input data, the `final_stats.csv` file should match the 
`ref_stats.csv` file.

## Assumptions

The scripts make certain assumptions:

- Sample names are unique
- Some stats will be generated for all files
- No intermediate script will fail catastrophically

These are safe assumptions for our purposes as far as we know, but if not in
effect the scripts may behave unexpectedly.

Feel free to modify the code as needed. The code and tooling is licensed under
the Creative Commons with Attribution license.
