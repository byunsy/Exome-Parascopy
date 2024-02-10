""" ===========================================================================
TITLE   : run_getBamCounts.py
AUTHOR  : Sang Yoon Byun
PURPOSE : runs ExomeDepth's getBamCounts() concurrently using asyncio
USAGE   : python run_getBamCounts.py -i BAMLIST -o OUTDIR -x FILE -t INT
=========================================================================== """

import os
import sys
import glob
import pyreadr
import asyncio
import argparse
import subprocess
import pandas as pd

# -----------------------------------------------------------------------------
# Procedures
# -----------------------------------------------------------------------------

def parse_args():
    
    # Parse arguments from command-line
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", required=True, 
                        help="input BAM file list")
    parser.add_argument("-o", "--output", required=True, 
                        help="output directory")
    parser.add_argument("-x", "--exons", required=False, 
                        help="custom exons file")
    parser.add_argument("-t", "--threads", default=8, 
                        help="Number of threads")
    args_parsed = parser.parse_args()
    
    return args_parsed


async def proc_count(bam_fp, sem, exons, outdir):

    path = os.path.dirname(os.path.abspath(__file__))
    sample_id = bam_fp.split("/")[-1].split(".")[0]

    async with sem:
        proc = await asyncio.create_subprocess_exec(
                'Rscript', path + '/run_ExomeDepthCount.r', 
                '-s', bam_fp,                       # input sample BAM
                '-o', outdir,                       # output directory
                '-p', sample_id,                    # prefix for output file
                '-x', exons,                        # path to exons file
                '-v',                               # verbosity
                stdout=asyncio.subprocess.PIPE, 
                stderr=asyncio.subprocess.PIPE)
            
        stdout, stderr = await proc.communicate()
        
        if stdout:
            print(stdout.decode())
        if stderr:
            print(stderr.decode())


def proc_merge(outdir):

    rds_files = glob.glob(f"{outdir}/counts_df_*.rds")
    rds_files.sort()
    rds_list = [pyreadr.read_r(f)[None] for f in rds_files]

    all_counts = pd.concat(rds_list, axis=1)
    all_counts.astype(int).to_csv(f"{outdir}/all.counts.tsv", sep="\t", index=None)


# -----------------------------------------------------------------------------
# Main Function
# -----------------------------------------------------------------------------
async def main():
    
    # Parse args
    args = parse_args()

    # Get a list of input BAM filepaths 
    with open(args.input, "r") as listfile:
        data = listfile.read()
        f_list = data.splitlines()

    # Asynchronous subprocesses
    MAX_PROCESSES = int(args.threads)
    RESULTS_DIR = args.output
    EXONS_FP = args.exons

    print(f"Asynchronous getBamCounts() started with {MAX_PROCESSES} processes.\n")

    sem = asyncio.Semaphore(MAX_PROCESSES)
    await asyncio.gather(*[proc_count(bam_fp, sem, EXONS_FP, RESULTS_DIR) for bam_fp in f_list])
    
    proc_merge(RESULTS_DIR)

if __name__ == "__main__":
    asyncio.run(main())

