import sys
import argparse
import pandas as pd


parser = argparse.ArgumentParser()
group_inp = parser.add_mutually_exclusive_group()
group_inp.add_argument("-i", "--input", help="path to final output file")
group_inp.add_argument("-I", "--input-list", nargs='+', default=[], help="a list of paths to final output files")
parser.add_argument("-c", "--refcn", required=True, help="reference copy number")
args = parser.parse_args()

# Parse arguments from command-line
REF_CN = int(args.refcn)
fps = [args.input]
if args.input_list:
    fps = args.input_list

# Get data from command-line
#fp, REF_CN = sys.argv[1], int(sys.argv[2])

total_dup_tp, total_dup_fp, total_dup_fn = 0, 0, 0
total_del_tp, total_del_fp, total_del_fn = 0, 0, 0

for fp in fps:

    # Prepare dataframe
    df = pd.read_csv(fp, sep=" ")
    df = df.reset_index()
    df.columns = ['final-CN', 'component', 'component_size', 'sample_name', 
                  'integerCN', 'penalizedCN', 'fractionalCN', 'trueCN']
    df = df.loc[df.integerCN != 'likelihood']
    df = df.loc[df.trueCN != 0]
    df = df.reset_index()
    df["integerCN"] = pd.to_numeric(df["integerCN"])
    df["trueCN"] = pd.to_numeric(df["trueCN"])
    
    df['dup_tp'] = (df.integerCN == df.trueCN) & (df.integerCN > REF_CN)
    df['del_tp'] = (df.integerCN == df.trueCN) & (df.integerCN < REF_CN)
    
    df['dup_fp'] = (df.integerCN > df.trueCN) & (df.trueCN == REF_CN)
    df['del_fp'] = (df.integerCN < df.trueCN) & (df.trueCN == REF_CN)
    
    df['dup_fn'] = (REF_CN < df.trueCN) & (df.integerCN == REF_CN)
    df['del_fn'] = (REF_CN > df.trueCN) & (df.integerCN == REF_CN)

    # Summarize stats
    print(fp)
    print(f"DUP: TP {df['dup_tp'].sum()}, FP {df['dup_fp'].sum()}, FN {df['dup_fn'].sum()}")
    print(f"DEL: TP {df['del_tp'].sum()}, FP {df['del_fp'].sum()}, FN {df['del_fn'].sum()}")
    print()
    
    total_dup_tp += df['dup_tp'].sum()
    total_dup_fp += df['dup_fp'].sum()
    total_dup_fn += df['dup_fn'].sum()
    
    total_del_tp += df['del_tp'].sum()
    total_del_fp += df['del_fp'].sum()
    total_del_fn += df['del_fn'].sum()

print("Duplication Events")
print("- Precis:", total_dup_tp / (total_dup_tp + total_dup_fp))
print("- Recall:", total_dup_tp / (total_dup_tp + total_dup_fn))
print()
print("Deletion Events")
print("- Precis:", total_del_tp / (total_del_tp + total_del_fp))
print("- Recall:", total_del_tp / (total_del_tp + total_del_fn))
