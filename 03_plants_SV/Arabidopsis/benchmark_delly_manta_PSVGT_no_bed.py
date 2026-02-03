import pandas as pd
import numpy as np
import sys
import argparse
import os

def get_header_line_index(filepath):
    """Returns the zero-based index of the header line (first line not starting with ##)."""
    with open(filepath, 'r') as f:
        for i, line in enumerate(f):
            if not line.startswith('##'):
                return i
    return 0

def load_delly_vcf(vcf_path):
    print(f"Loading Delly VCF: {vcf_path}")
    header = get_header_line_index(vcf_path)
    # Use low_memory=False to prevent dtype warnings on large files
    vcf_df = pd.read_csv(vcf_path, header=header, sep="\t", dtype=str)
    
    df_fil = vcf_df.copy()
    # Extract Genotype (GT)
    df_fil['GT'] = df_fil[df_fil.columns[-1]].str.split(":", expand=True)[0]
    df_fil = df_fil[df_fil['GT'] != "0/0"]
    
    df_fil["POS"] = pd.to_numeric(df_fil["POS"], errors='coerce').astype(int)
    # Parse INFO field for SVTYPE
    df_fil['SVTYPE'] = df_fil['INFO'].str.extract(r'SVTYPE=([^;]+)', expand=False)
    if df_fil['SVTYPE'].isnull().all():
         df_fil['SVTYPE'] = df_fil['ID'].astype(str).str[:3]

    # Process Deletions
    dels = df_fil[df_fil['SVTYPE'] == 'DEL'].copy()
    dels["END"] = pd.to_numeric(dels['INFO'].str.extract(r'END=([^;]+)', expand=False), errors='coerce')
    dels['SVLEN'] = (dels["END"] - dels['POS']).abs()
    dels = dels[dels['SVLEN'] >= 45]

    # Process Insertions
    ins = df_fil[df_fil['SVTYPE'] == 'INS'].copy()
    ins['SVLEN_STR'] = ins['INFO'].str.extract(r'SVLEN=([^;]+)', expand=False)
    ins['SVLEN'] = pd.to_numeric(ins['SVLEN_STR'], errors='coerce').abs()
    
    # Fallback if SVLEN is missing in INFO
    if ins['SVLEN'].isnull().any():
         ins["END"] = pd.to_numeric(ins['INFO'].str.extract(r'END=([^;]+)', expand=False), errors='coerce')
         
    ins = ins[ins['SVLEN'] >= 45]
    
    indel = pd.concat([dels, ins], axis=0, ignore_index=True)
    indel.sort_values(by=["#CHROM", "POS"], inplace=True)
    
    # Create Unique ID
    indel['SVID'] = (indel["#CHROM"].astype(str) + ":" + 
                     indel['POS'].astype(str) + "_" + 
                     indel["SVTYPE"] + "=" + 
                     indel['SVLEN'].astype(str))
    
    return indel

def load_manta_vcf(vcf_path):
    print(f"Loading Manta VCF: {vcf_path}")
    header = get_header_line_index(vcf_path)
    vcf = pd.read_csv(vcf_path, header=header, sep="\t", dtype=str)
    
    vcf['GT'] = vcf[vcf.columns[-1]].str.split(":", expand=True)[0]
    vcf = vcf[vcf['GT'] != "0/0"]
    vcf['POS'] = pd.to_numeric(vcf['POS']).astype(int)
    
    # Extract SVLEN and SVTYPE
    vcf['SVLEN'] = pd.to_numeric(vcf['INFO'].str.extract(r'SVLEN=([^;]+)', expand=False), errors='coerce').abs()
    vcf = vcf.dropna(subset=['SVLEN'])
    vcf = vcf[vcf['SVLEN'] >= 45]
    
    vcf['SVTYPE'] = vcf['INFO'].str.extract(r'SVTYPE=([^;]+)', expand=False)
    
    vcf['SVID'] = (vcf["#CHROM"].astype(str) + ":" + 
                   vcf['POS'].astype(str) + "_" + 
                   vcf["SVTYPE"] + "=" + 
                   vcf['SVLEN'].astype(int).astype(str))
    return vcf

def load_generic_vcf(vcf_path):
    print(f"Loading Generic VCF: {vcf_path}")
    header = get_header_line_index(vcf_path)
    vcf = pd.read_csv(vcf_path, header=header, sep="\t", dtype=str)
    
    vcf['POS'] = pd.to_numeric(vcf['POS']).astype(int)
    vcf.sort_values(by=['#CHROM', 'POS'], inplace=True)
    
    # GT Processing
    vcf['GT'] = vcf[vcf.columns[-1]].str.split(":", expand=True)[0]
    vcf = vcf[~vcf['GT'].isin(["0/0", "./."])]
    vcf["GT"] = vcf["GT"].str.replace("|", "/")
    
    def process_gt(gt):
        if pd.isna(gt): return gt
        alleles = gt.split("/")
        return "0/1" if len(set(alleles)) == 2 else gt
        
    vcf["GT"] = vcf["GT"].apply(process_gt)

    # INFO Processing
    vcf['SVLEN'] = pd.to_numeric(vcf['INFO'].str.extract(r'SVLEN=([^;]+)', expand=False), errors='coerce').abs()
    vcf = vcf.dropna(subset=['SVLEN'])
    vcf = vcf[vcf['SVLEN'] >= 45]
    
    vcf['SVTYPE'] = vcf['INFO'].str.extract(r'SVTYPE=([^;]+)', expand=False)
    
    vcf['SVID'] = (vcf["#CHROM"].astype(str) + ":" + 
                   vcf['POS'].astype(str) + "_" + 
                   vcf["SVTYPE"] + "=" + 
                   vcf['SVLEN'].astype(int).astype(str))
    return vcf

def find_closest_row(df, target_number, col_name='POS'):
    """Finds the row in df with the value in col_name closest to target_number."""
    if df.empty:
        return None
    # Calculate absolute difference
    diff = (df[col_name] - target_number).abs()
    closest_index = diff.idxmin()
    return df.loc[closest_index]

def main():
    parser = argparse.ArgumentParser(description="Benchmark SV VCF against a truth set VCF (All variants).")
    
    parser.add_argument("-b", "--bench", required=True, help="Path to Benchmark/Truth VCF file")
    parser.add_argument("-i", "--input", required=True, help="Path to Input/Test VCF file")
    parser.add_argument("-t", "--type", choices=['delly', 'manta', 'generic'], default='generic', 
                        help="Type of input VCF (delly, manta, or generic default)")
    parser.add_argument("-o", "--output", default="comparison_result", help="Prefix for output files")

    args = parser.parse_args()

    # 1. Load Data
    print("--- Loading Files ---")
    bench_df = load_generic_vcf(args.bench)
    
    if args.type == 'delly':
        vcf2_df = load_delly_vcf(args.input)
    elif args.type == 'manta':
        vcf2_df = load_manta_vcf(args.input)
    else:
        vcf2_df = load_generic_vcf(args.input)

    print(f"Benchmark SVs (Total): INS={len(bench_df[bench_df['SVTYPE']=='INS'])}, DEL={len(bench_df[bench_df['SVTYPE']=='DEL'])}")
    print(f"Test SVs (Total): INS={len(vcf2_df[vcf2_df['SVTYPE']=='INS'])}, DEL={len(vcf2_df[vcf2_df['SVTYPE']=='DEL'])}")

    # 2. Comparison Logic
    print("--- Comparing SVs ---")
    
    # Initialize Counters
    stats = {k: {'true': 0, 'false': 0, 'bad_gt': 0, 'total_bench': 0} for k in ['DEL', 'INS', 'DUP', 'INV']}
    
    for k in stats:
        stats[k]['total_bench'] = len(bench_df[bench_df['SVTYPE'] == k])

    true_ids = set()

    # Iterate through Test VCF
    for idx, row in vcf2_df.iterrows():
        chrom = row['#CHROM']
        pos = int(row['POS'])
        svlen = row['SVLEN']
        svtype = row['SVTYPE']
        gt = row['GT']
        
        # Define window for matching
        window = 1000
        
        # Filter benchmark for candidate matches
        if svtype == 'DEL':
            vcf_end = pos + svlen
            # Overlap logic: Bench Start < Test End AND Bench End > Test Start (simplified)
            candidates = bench_df[
                (bench_df['#CHROM'] == chrom) &
                (bench_df['SVTYPE'] == 'DEL') &
                (bench_df['POS'] < vcf_end + window) & 
                ((bench_df['POS'] + bench_df['SVLEN']) > pos - window)
            ]
        else:
            # For INS/INV/DUP, usually look near the insertion point
            candidates = bench_df[
                (bench_df['#CHROM'] == chrom) &
                (bench_df['SVTYPE'] == svtype) &
                (bench_df['POS'] > pos - window) &
                (bench_df['POS'] < pos + window)
            ]

        if candidates.empty:
            if gt != "0/1": 
                 if svtype in stats: stats[svtype]['false'] += 1
        else:
            # Find best match
            nearest = find_closest_row(candidates, pos, 'POS')
            
            bench_len = nearest['SVLEN']
            bench_gt = nearest['GT']
            bench_id = nearest['SVID']
            
            # Size Ratio Check (0.5 < ratio < 2.0)
            ratio = svlen / bench_len if bench_len > 0 else 0
            if ratio < 0.5 or ratio > 2.0:
                 if svtype in stats: stats[svtype]['false'] += 1
            else:
                # It is a True Positive
                true_ids.add(bench_id)
                if svtype in stats: stats[svtype]['true'] += 1
                
                # Check Genotype
                if bench_gt != gt:
                    if svtype in stats: stats[svtype]['bad_gt'] += 1

    # 3. Calculate Metrics
    print(f"\n--- Results: {args.input} vs {args.bench} (Whole Genome) ---")
    print(f"SV\tRecall\tPrecision\tF1\tGT_Concordance")

    total_tp = 0
    total_fp = 0
    total_bench_count = 0
    total_bad_gt = 0

    for svtype, data in stats.items():
        tp = data['true']
        fp = data['false']
        # Note: In whole genome mode, all benchmark SVs are considered "Targets"
        fn = data['total_bench'] - tp 
        bad_gt = data['bad_gt']
        
        recall = tp / data['total_bench'] if data['total_bench'] > 0 else 0
        precision = tp / (tp + fp) if (tp + fp) > 0 else 0
        f1 = 2 * (precision * recall) / (precision + recall) if (precision + recall) > 0 else 0
        gt_conc = (tp - bad_gt) / tp if tp > 0 else 0
        
        print(f"{svtype}\t{recall:.2%}\t{precision:.2%}\t{f1:.2%}\t{gt_conc:.2%}")

        total_tp += tp
        total_fp += fp
        total_bench_count += data['total_bench']
        total_bad_gt += bad_gt

    # Overall
    ov_recall = total_tp / total_bench_count if total_bench_count > 0 else 0
    ov_prec = total_tp / (total_tp + total_fp) if (total_tp + total_fp) > 0 else 0
    ov_f1 = 2 * (ov_prec * ov_recall) / (ov_prec + ov_recall) if (ov_prec + ov_recall) > 0 else 0
    ov_gt_conc = (total_tp - total_bad_gt) / total_tp if total_tp > 0 else 0
    
    print(f"Total\t{ov_recall:.2%}\t{ov_prec:.2%}\t{ov_f1:.2%}\t{ov_gt_conc:.2%}")

    # Output False Negatives list
    fn_df = bench_df[~bench_df['SVID'].isin(true_ids)]
    fn_outfile = f"{args.output}_fn.list"
    fn_df.to_csv(fn_outfile, sep="\t", index=False)
    print(f"\nFalse Negatives saved to: {fn_outfile}")

if __name__ == "__main__":
    main()
