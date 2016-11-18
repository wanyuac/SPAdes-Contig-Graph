#!/usr/bin/env python

"""
This script submits jobs of spades_contig_graph.py through SLURM.

Usage:
    python spades_contig_graph_slurm.py --fasta *_spades.fasta --fastg *_spades.fastg --paths *_spades.paths \
    --suffix_fasta _spades.fasta --suffix_fastg _spades.fastg --suffix_paths _spades.paths -c --suffix_out c
	
Note: increase the "mem" value if the BLAST database error "No alias or index file found for nucleotide database" occurs when
you run this script under a "-c" mode.
    
Dependency: BLAST
    
Author: Yu Wan (wanyuac@gmail.com, https://github.com/wanyuac)
Python v2.7+
Development: 17/11/2016
License: GNU GPL 2.1
"""

import sys
import os
import subprocess
import collections
from argparse import ArgumentParser

"""
Currently, this script is designed to work on the computational cluster Snowy of VLSCI.
Nonetheless, it is easy to adapt this script to other SLURM-running systems through configuring the following constants.
Set an element an empty string ("") if any module is not running as a module on your computer.
BLAST is necessary when "-c" option is laid.
"""

PARTITION = "sysgen"  # the partition option for SLURM
MODULES = ["Python/2.7.10-vlsci_intel-2015.08.25-SG", "BLAST+/2.2.30-vlsci_intel-2015.08.25-Python-2.7.10"]

def parse_arguments():
    parser = ArgumentParser(description = "Generate contig graphs from assembly graphs")
    parser.add_argument("--fasta", nargs = "+", type = str, required = True, help = "FASTA files produced by SPAdes")
    parser.add_argument("--fastg", nargs = "+", type = str, required = True, help = "FASTG files produced by SPAdes")
    parser.add_argument("--paths", nargs = "+", type = str, required = True, help = "Path files produced by SPAdes")
    parser.add_argument("--suffix_fasta", type = str, required = False, default = "_spades.fasta", \
                        help = "Suffix that will be used to extract sample names from file names")
    parser.add_argument("--suffix_fastg", type = str, required = False, default = "_spades.fastg", \
                        help = "Suffix that will be used to extract sample names from file names")
    parser.add_argument("--suffix_paths", type = str, required = False, default = "_spades.paths", \
                        help = "Suffix that will be used to extract sample names from file names")
    parser.add_argument("-c", action = "store_true", help = "set to prioritise graph connections over segment length")
    parser.add_argument("--suffix_out", type = str, required = False, default = "ctg", \
                        help = "Suffix that will be appended to sample names to form an output file name.")
    parser.add_argument("--script_path", type = str, required = False, default = ".", \
                        help = "Path to spades_contig_graph.py")
    parser.add_argument("--mem", type = str, required = False, default = "4096", help = "Memory assignment (Mb) for every job.")
    return parser.parse_args()

def check_inputs(fasta, fastg, paths):
    n1 = len(fasta)
    n2 = len(fastg)
    n3 = len(paths)
    if (n1 != n2 or n1 != n3 or n2 != n3):
        sys.exit("Error: input contig, graph and paths files must of the same number.")
    return

def get_sample_name(fn, sf):
    """
    Extracts a sample name from a path (fn) using the suffix (sf) of the file name.
    For instance, ERR562360 will be extracted from the string "~/assemblies/ERR562360_spades.fastg".
    """
    return os.path.basename(fn).rstrip(sf)
    
def organise_files(fasta, fastg, paths, sf_fasta, sf_fastg, sf_paths):
    """
    This function returns a two dimensional dictionary with the hierarchical structure: [sample name][fastg/fasta/paths] = file path.
    Herein I assume samples are the same for FASTA, FASTG and paths files. Otherwise, an error will rise.
    """
    files = collections.defaultdict(dict)
    
    # store information about FASTA files first
    for f in fasta:
        try:
            sample = get_sample_name(f, sf_fasta)
            files[sample] = {"contigs" : f, "graph" : "", "paths" : ""}
        except NameError:
            print("Error: there are problems in sample names within names of FASTA files.")
            raise
    
    # process names of FASTG files
    for f in fastg:
        try:
            sample = get_sample_name(f, sf_fastg)
            files[sample]["graph"] = f
        except NameError:
            print("Error: sample names within some FASTG files are not matched to those in FASTA files.")
            raise
        
    # process names of path files
    for f in paths:
        try:
            sample = get_sample_name(f, sf_paths)
            files[sample]["paths"] = f
        except NameError:
            print("Error: sample names within some path files are not matched to those in FASTA files or FASTG files.")
            raise
        
    return files

def submit_jobs(fs, is_c, par, modules, sf, script, mem):
    if is_c:
        style = "-c"  # to make connection-priority contig graphs
    else:
        style = "-l"  # to make length-priority contig graphs
    
    for sample, vals in fs.iteritems():
        outdir = os.path.dirname(os.path.abspath(vals["contigs"]))
        
        # generate SLURM commands
        cmd = "#!/bin/bash"
        if len(par) > 0:
            cmd += "\n#SBATCH -p " + par  # otherwise, use the system"s default partition of jobs
        cmd += "\n#SBATCH --job-name=ctg:" + sample
        cmd += "\n#SBATCH --ntasks=1"
        cmd += "\n#SBATCH --mem-per-cpu=" + mem
        cmd += "\n#SBATCH --time=0-1:30:0\n"
        for m in modules:
            if m != "":
                cmd += "\nmodule load " + m
        cmd += "\n\npython %s %s %s %s %s %s\n" % \
               (script, style, vals["graph"], vals["contigs"], vals["paths"], \
                os.path.join(outdir, sample + "_" + sf + ".fastg"))
        
        slurm_filename = sample + ".slurm"
        with open(slurm_filename, "w") as slurm_script:
            slurm_script.write(cmd)
            
        os.system("sbatch " + slurm_filename)  # submit this job

def main():
    args = parse_arguments()
    check_inputs(args.fasta, args.fastg, args.paths)
    files = organise_files(args.fasta, args.fastg, args.paths, args.suffix_fasta, args.suffix_fastg, args.suffix_paths)
    script = os.path.join(args.script_path, "spades_contig_graph.py")
    submit_jobs(files, args.c, PARTITION, MODULES, args.suffix_out, script, args.mem)
    
if __name__ == "__main__":
    main()
