""" ===========================================================================
TITLE   : run_getBamCounts.py
AUTHOR  : Sang Yoon Byun
PURPOSE : runs ExomeDepth's getBamCounts() concurrently using asyncio
USAGE   : python run_getBamCounts.py -i BAMLIST -o OUTDIR -x FILE -t INT
=========================================================================== """

import asyncio
import argparse
import subprocess

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
args = parser.parse_args()

# Get a list of input BAM filepaths 
with open(args.input, "r") as listfile:
    data = listfile.read()
    f_list = data.splitlines()

# Asynchronous subprocesses
MAX_PROCESSES = int(args.threads)
RESULTS_DIR = args.output

print(f"Asynchronous getBamCounts() started with {MAX_PROCESSES} processes.\n")

# -----------------------------------------------------------------------------
# Procedures
# -----------------------------------------------------------------------------
async def proc_count(bam_fp, sem):
  
    sample_id = bam_fp.split("/")[-1].split(".")[0]
    async with sem:
        proc = await asyncio.create_subprocess_exec(
                'Rscript', 'run_ExomeDepthCount.r', 
                '-s', bam_fp,                       # input sample BAM
                '-o', RESULTS_DIR,                  # output directory
                '-p', sample_id,                    # prefix for output file
                '-x', args.exons,                   # path to exons file
                '-v',                               # verbosity
                stdout=asyncio.subprocess.PIPE, 
                stderr=asyncio.subprocess.PIPE)
            
        stdout, stderr = await proc.communicate()
        
        if stdout:
            print(stdout.decode())
        if stderr:
            print(stderr.decode())


def proc_merge(outdir):
    subprocess.call([f'Rscript run_ExomeDepthCount_merge.r -o {outdir}'], shell=True)

# -----------------------------------------------------------------------------
# Main Function
# -----------------------------------------------------------------------------
async def main():
    sem = asyncio.Semaphore(MAX_PROCESSES)
    await asyncio.gather(*[proc_count(bam_fp, sem) for bam_fp in f_list])
    
    proc_merge(RESULTS_DIR)

if __name__ == "__main__":
    asyncio.run(main())

