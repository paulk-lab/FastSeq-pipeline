#!/usr/local/bin/python3
"""
Wrapper script to launch docker and process sequencing results for packaged
viral genomes.
"""
from argparse import ArgumentParser
from os import getcwd
from pathlib import Path
from subprocess import run

# Argument parsing setup
parser = ArgumentParser(description='Launch a docker container and process '
                                    'sequences with the same environment as '
                                    'the Paulk team.')

parser.add_argument('py_file', type=str,
                    help='The path to the sequence processing pipeline in '
                         'python. In most cases this is process_seq.py')

parser.add_argument('base_dir', type=str,
                    help='Base of where processing takes place. All paths '
                         'in the csv are assumed to be relative to this path '
                         'and results will be placed in "Output" directory '
                         'within this path.')

parser.add_argument('csv_file', type=str,
                    help='CSV file detailing samples, and where the relevant '
                         'files for those samples can be found, all paths '
                         'are relative to the base_dir.')

args = parser.parse_args()

working = Path(getcwd())
base = working / args.base_dir
csv = working / args.csv_file
pyscript = working / args.py_file

run(["docker", "run", "--rm",
     "-v", f"{pyscript}:/home/worker/process_seq.py",
     "-v", f"{base}:/home/worker/host",
     "-v", f"{csv}:/home/worker/samples.csv",
     "paulklab/fastseq-pipeline:version0.3"])
