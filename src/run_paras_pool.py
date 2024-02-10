""" ===========================================================================
TITLE   : run_paras_pool.py
AUTHOR  : Sang Yoon Byun
PURPOSE : aggregate reads of duplicated gene and extract read counts
=========================================================================== """

import os
import asyncio
import argparse
import subprocess

import run_getBamCounts as _bamcounts

# -----------------------------------------------------------------------------
# Procedures
# -----------------------------------------------------------------------------

def parse_args():
    
    # Parse arguments from command-line
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", required=True, 
                        help="path to input BAM file list")
    parser.add_argument("-o", "--output", required=True, 
                        help="path to output directory")
    parser.add_argument("-f", "--reference", required=True, 
                        help="path to reference FASTA")
    parser.add_argument("-t", "--hom-table", required=True,
                        help="path to homology table")
    parser.add_argument("-x", "--exons", required=True,
                        help="path to full list of exons")
    parser.add_argument("-l", "--loci-name", required=True,
                        help="path to background depth regions")
    parser.add_argument("-r", "--loci-regions", required=True,
                        help="path to file with loci of interest")
    parser.add_argument("-g", "--hg-version", required=True,
                        help="hg19 or hg38")
    parser.add_argument("-@", "--threads", default=8, 
                        help="Number of threads")
    args_parsed = parser.parse_args()

    return args_parsed


def create_custom_exon_file(loci, hg_ver, exons_fp):
    
    exons_dir = os.path.dirname(exons_fp)
    if not exons_dir:
        exons_dir = os.path.dirname(os.path.abspath(__file__))

    # Obtain exons of loci of interest
    exons_bed = f'{exons_dir}/exons.{hg_ver}.{loci.upper()}.bed'

    if not os.path.isfile(exons_bed):
        subprocess.call([f'head -n 1 {exons_fp} > {exons_bed}'], shell=True)    
        subprocess.call([f'grep -P "\t{loci.upper()}_" {exons_fp} >> {exons_bed}'], shell=True)    
        print(f"Created custom exon file for {loci}.")
    else:
        print(f"Found an existing exon file: {exons_bed}.")

    return exons_bed


def create_output_dirs(output_dir, loci_name):
    
    # results output directory
    if not os.path.isdir(output_dir):
        subprocess.call([f'mkdir {output_dir}'], shell=True)    
    
    # pooled_reads directory
    pooled_reads_dir = f'pooled-reads-bam-{loci_name.upper()}'
    if not os.path.isdir(pooled_reads_dir):
        subprocess.call([f'mkdir {pooled_reads_dir}'], shell=True)

    print("Created output directory.")
    return pooled_reads_dir


def create_input_lists(input_list_fp, outdir):
    """ Read input filepath list and get input sample names """

    with open(input_list_fp, 'r') as inp:
        input_list = inp.read().splitlines()
    
    sample_list = [line.split("::")[1] for line in input_list]

    with open(f"{outdir}/input.sample.list", "w") as sample_out:
        sample_out.writelines("\n".join(sample_list))

    return f"{outdir}/input.sample.list"
        

async def run_paras_pool(inp_bam, args, sem):
    """ Run Parascopy pool on the specified loci """
    
    path = os.path.dirname(os.path.abspath(__file__))

    async with sem:
        proc = await asyncio.create_subprocess_exec(
                f'{path}/scripts/paras_pool.sh', 
                inp_bam,
                args.loci_name.upper(),
                args.loci_regions,
                args.hom_table,
                args.reference,
                stdout=asyncio.subprocess.PIPE,
                stderr=asyncio.subprocess.PIPE)

        stdout, stderr = await proc.communicate()

        if stdout:
            print(stdout.decode())
        if stderr:
            print(stderr.decode())
    

async def run_sort_fix_sort(sample_id, pool_dir, sem):
    
    path = os.path.dirname(os.path.abspath(__file__))

    async with sem:
        proc = await asyncio.create_subprocess_exec(
                f'{path}/scripts/sort_fix_sort.sh', sample_id, pool_dir,
                stdout=asyncio.subprocess.PIPE,
                stderr=asyncio.subprocess.PIPE)

        stdout, stderr = await proc.communicate()

        if stdout:
            print(stdout.decode())
        if stderr:
            print(stderr.decode())

    subprocess.call([f'ls -d $PWD/{pool_dir}/*.final.bam > {pool_dir}/pooled.fp.list'], shell=True)


# -----------------------------------------------------------------------------
# Main Function
# -----------------------------------------------------------------------------

async def main():
    
    # Parse arguments
    args = parse_args()
    
    # Prepare necessary files and directories
    exons_fp = create_custom_exon_file(args.loci_name, args.hg_version, args.exons)
    pool_dir = create_output_dirs(args.output, args.loci_name)
    inp_samples_fp = create_input_lists(args.input, args.output)
    
    # list of input fp
    with open(args.input, "r") as inpf1:
        input_fp = inpf1.read().splitlines()

    # list of input samples
    with open(inp_samples_fp, "r") as inpf2:
        inp_samples = inpf2.read().splitlines()
       
    MAX_PROCESSES = int(args.threads)
    sem = asyncio.Semaphore(MAX_PROCESSES)

    # Run Parascopy pool
    await asyncio.gather(*[run_paras_pool(f, args, sem) for f in input_fp])

    # Run sort_fix_sort
    await asyncio.gather(*[run_sort_fix_sort(sample_id, pool_dir, sem) for sample_id in inp_samples])
    
    # list of pooled reads BAM files
    with open(f"{pool_dir}/pooled.fp.list", "r") as listfile:
        pooled_f_list = listfile.read().splitlines()
    
    # run getBamCounts
    print(f"Asynchronous getBamCounts() started with {MAX_PROCESSES} processes.\n")
    await asyncio.gather(*[_bamcounts.proc_count(bam_fp, sem, exons_fp, args.output) for bam_fp in pooled_f_list]) 
    _bamcounts.proc_merge(args.output)

    print(f"Completed extracting pooled read counts for {args.loci_name}.")
    

if __name__ == "__main__":
    asyncio.run(main())

