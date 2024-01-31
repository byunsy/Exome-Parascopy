import os
import argparse
import subprocess


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


def create_custom_exon_file(loci, hg_ver):
    
    # Obtain exons of loci of interest
    exons_bed = f'exons/exons.{hg_ver}.{loci.upper()}.bed'

    if not os.path.isfile(exons_bed):
        subprocess.call([f'head -n 1 exons/exons.{hg_ver}.noalt.bed > {exons_bed}'], shell=True)    
        subprocess.call([f'grep -P "\t{loci.upper()}_" exons/exons.{hg_ver}.noalt.bed >> {exons_bed}'], shell=True)    
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


def create_input_lists(input_list_fp):
    with open(input_list_fp, 'r') as inp:
        input_list = inp.read().splitlines()
    
    sample_list = [line.split("::")[1] for line in input_list]

    with open("input.sample.list", "w") as sample_out:
        sample_out.writelines("\n".join(sample_list))

    return "input.sample.list"
        

def run_paras_pool(args):
    
    command = (
        f'cat {args.input} | '
        f'parallel -j {args.threads} '
        f'./scripts/paras_pool.sh {{}} {args.loci_name.upper()} {args.loci_regions} {args.hom_table} {args.reference}'
    )
    subprocess.call([command], shell=True)    
    #print(command)


def run_sort_fix_sort(sample_list, pool_dir, num_thread):
    subprocess.call([f'cat {sample_list} | parallel -j {num_thread} ./scripts/sort_fix_sort.sh {{}} {pool_dir}'], shell=True)
    subprocess.call([f'ls -d $PWD/{pool_dir}/*.final.bam > {pool_dir}/pooled.fp.list'], shell=True)


def run_getBamCounts(pool_dir, args, exons):
    
    command = (
        f'time python run_getBamCounts.py '
        f'-i {pool_dir}/pooled.fp.list '
        f'-o {args.output} '
        f'-x {exons} '
        f'-t {args.threads}'
    )
    subprocess.call([command], shell=True)


def main():
    
    # Parse arguments
    args = parse_args()
    
    # Prepare necessary files and directories
    exons_fp = create_custom_exon_file(args.loci_name, args.hg_version)
    pool_dir = create_output_dirs(args.output, args.loci_name)
    inp_samples_fp = create_input_lists(args.input)

    # Run Parascopy pool
    run_paras_pool(args)

    # Run sort_fix_sort
    run_sort_fix_sort(inp_samples_fp, pool_dir, args.threads)

    # run_getBamCounts.py
    run_getBamCounts(pool_dir, args, exons_fp)


if __name__ == "__main__":
    main()
