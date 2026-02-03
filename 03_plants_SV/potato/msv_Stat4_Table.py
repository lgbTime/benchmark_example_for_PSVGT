import re
import argparse
import pandas as pd


def parse_input_line(line):
    """Parse input line and return extended interval, SV information, and allele count"""
    parts = line.strip().split('\t')[1:]  # Skip first column (likely ID)
    starts, ends, sv_types, sv_sizes = [], [], [], []
    for part in parts:
        # Regex pattern to extract SV details (chrom:start-end_type=size)
        match = re.match(r'chr\d+:(?P<start>\d+)-(?P<end>\d+)_(?P<type>DEL|INS)=(?P<size>\d+)', part)
        if match:
            starts.append(int(match.group('start')))
            ends.append(int(match.group('end')))
            sv_types.append(match.group('type'))
            sv_sizes.append(int(match.group('size')))
    
    if not starts:  # Handle empty lines
        return 0, 0, [], [], 0  # Allele count = 0
    
    # Extend interval with boundary protection
    min_start = max(1, min(starts) - 300)
    max_end = max(ends) + 300
    if max_end <= min_start:
        max_end = min_start + 1  # Avoid invalid interval
    
    allele_num = len(sv_types)  # Allele count equals number of SVs
    return min_start, max_end, sv_types, sv_sizes, allele_num


def parse_vcf_file(vcf_file):
    """Parse VCF file into DataFrame (consistent with original logic)"""
    data = []
    with open(vcf_file, 'r') as f:
        for line in f:
            if line.startswith('#'):  # Skip header lines
                continue
            fields = line.strip().split('\t')
            info = fields[7]
            start = int(fields[1])
            sv_type = sv_size = end = None
            # Extract SV attributes from INFO field
            for item in info.split(';'):
                if item.startswith('SVTYPE='):
                    sv_type = item.split('=')[1]
                elif item.startswith('END='):
                    end = int(item.split('=')[1])
                elif item.startswith('SVLEN='):
                    sv_size = abs(int(item.split('=')[1]))  # Use absolute value for length
            if sv_type and end and sv_size:  # Only keep valid entries
                data.append([start, end, sv_type, sv_size])
    return pd.DataFrame(data, columns=['start', 'end', 'sv_type', 'sv_size'])


def find_matching_sv(input_file, vcf_file):
    """Find matching SVs and return detailed results (preserves original statistics)"""
    results = []
    vcf_df = parse_vcf_file(vcf_file)
    
    # Original statistics variables
    total_line_count = 0     # Total lines in input
    matched_line_count = 0   # Lines with matches
    # Additional: Statistics for INS and DEL types
    ins_total = 0
    ins_matched = 0
    del_total = 0
    del_matched = 0

    with open(input_file, 'r') as in_f:
        for input_line in in_f:
            # Parse input line (includes allele count)
            min_start, max_end, sv_types, sv_sizes, allele_num = parse_input_line(input_line)
            
            if not sv_types:  # Skip lines with no SVs
                continue
                
            total_line_count += 1
            input_line_clean = input_line.strip()  # Original input line (without newline)
            # Assume single SV type per line, take first type
            sv_type = sv_types[0]
            if sv_type == 'INS':
                ins_total += 1
            else:
                del_total += 1
            
            # Filter VCF SVs within extended interval (original logic)
            filtered_df = vcf_df[
                (vcf_df['start'] >= min_start) & (vcf_df['start'] <= max_end) &
                (vcf_df['end'] >= min_start) & (vcf_df['end'] <= max_end)
            ]
            
            # Matching logic (original criteria)
            input_sv_matched = [False] * len(sv_types)
            matching_sv = []
            for i in range(len(sv_types)):
                input_type = sv_types[i]
                input_size = sv_sizes[i]
                # Filter VCF by SV type
                type_filtered = filtered_df[filtered_df['sv_type'] == input_type]
                for _, row in type_filtered.iterrows():
                    # Size similarity check (70-130% of input size)
                    size_rate = input_size / row['sv_size']
                    if 0.7 <= size_rate <= 1.3 and not input_sv_matched[i]:
                        if (row['start'], row['end'], row['sv_type'], row['sv_size']) not in matching_sv:
                            matching_sv.append((row['start'], row['end'], row['sv_type'], row['sv_size']))
                            input_sv_matched[i] = True
            
            # Determine if matched (original logic: ≥2 SV matches)
            is_matched = sum(input_sv_matched) >= 2
            if is_matched:
                matched_line_count += 1
                if sv_type == 'INS':
                    ins_matched += 1
                else:
                    del_matched += 1
            
            # Save results (includes allele count and match status)
            results.append({
                'input_line': input_line_clean,
                'matching_sv': matching_sv,
                'allele_num': allele_num,
                'is_matched': is_matched,
                'original_line': input_line  # Preserve original line for output file
            })
    
    return results, total_line_count, matched_line_count, ins_total, ins_matched, del_total, del_matched


def main():
    parser = argparse.ArgumentParser(description='Multi-allele SV matching analysis (preserves original stats + categorized output)')
    parser.add_argument('input_file', type=str, help='Path to input file')
    parser.add_argument('vcf_file', type=str, help='Path to VCF file')
    args = parser.parse_args()

    # Run core matching function (preserves original statistics)
    results, total_line_count, matched_line_count, ins_total, ins_matched, del_total, del_matched = find_matching_sv(args.input_file, args.vcf_file)
    
    # ========== Original statistics printing ==========
    print(f"===== Matching Statistics =====")
    print(f"Total lines in input file: {total_line_count}")
    print(f"Lines with matches in VCF: {matched_line_count}")
    print(f"Matching rate: {matched_line_count/total_line_count:.2%}")
    print("===================\n")

    # Print detailed matching results (consistent with original logic)
    for result in results:
        print(f"Input line: {result['input_line']}")
        if result['matching_sv']:
            print("Matching SVs in VCF:")
            for sv in result['matching_sv']:
                print(f"  Start: {sv[0]}, End: {sv[1]}, Type: {sv[2]}, Size: {sv[3]}")
        else:
            print("  No matching SVs found.")
        print()

    # ========== Additional: Create output file with 6th column "yes"/"No" ==========
    output_file = "sv_matching_results.txt"
    with open(output_file, 'w') as out_f:
        with open(args.input_file, 'r') as in_f:  # Read original input to preserve all columns
            for line in in_f:
                line_clean = line.strip()
                # Find corresponding result to check match status
                match_result = next((r for r in results if r['original_line'].strip() == line_clean), None)
                if match_result:
                    # Split original line into columns
                    columns = line_clean.split('\t')
                    # Ensure at least 5 columns exist (pad with empty if necessary)
                    while len(columns) < 5:
                        columns.append('')
                    # Add 6th column: "yes" if matched, else "No"
                    columns.append("YES" if match_result['is_matched'] else "NO")
                    # Write modified line to output file
                    out_f.write('\t'.join(columns) + '\n')
                else:
                    # Write original line if no match found (add "No" in 6th column)
                    columns = line_clean.split('\t')
                    while len(columns) < 5:
                        columns.append('')
                    columns.append("No")
                    out_f.write('\t'.join(columns) + '\n')
    print(f"Additional output file created: {output_file}")

    # ========== Additional: 3-allele/4-allele categorized output ==========
    # Initialize storage for categorized results
    output_tables = {
        '3_allele_matched': [],    # 3-allele with matches
        '3_allele_unmatched': [],  # 3-allele without matches
        '4_allele_matched': [],    # 4-allele with matches
        '4_allele_unmatched': []   # 4-allele without matches
    }

    for result in results:
        allele_num = result['allele_num']
        is_matched = result['is_matched']
        input_line = result['input_line']
        matching_sv = result['matching_sv']

        # Only process 3 or 4 allele cases
        if allele_num not in (3, 4):
            continue

        # Build row data for output tables
        row = {
            'Input_Line': input_line,
            'Matched_SVs': '; '.join([f'Start:{s[0]},End:{s[1]},Type:{s[2]},Size:{s[3]}' for s in matching_sv]) 
                          if matching_sv else 'No Matches'
        }

        if allele_num == 3:
            if is_matched:
                output_tables['3_allele_matched'].append(row)
            else:
                output_tables['3_allele_unmatched'].append(row)
        elif allele_num == 4:
            if is_matched:
                output_tables['4_allele_matched'].append(row)
            else:
                output_tables['4_allele_unmatched'].append(row)

    # Save categorized results to CSV files
    for key, data in output_tables.items():
        df = pd.DataFrame(data)
        df.to_csv(f'{key}.csv', index=False)
    print("Categorized CSV files created for 3-allele and 4-allele cases")

    # ========== Additional: INS and DEL matching statistics ==========
    ins_match_rate = ins_matched / ins_total if ins_total > 0 else 0
    del_match_rate = del_matched / del_total if del_total > 0 else 0
    print("\n===== INS and DEL Matching Statistics =====")
    print(f"Total INS lines: {ins_total}")
    print(f"Matched INS lines: {ins_matched}")
    print(f"INS matching rate: {ins_match_rate:.2%}")
    print(f"Total DEL lines: {del_total}")
    print(f"Matched DEL lines: {del_matched}")
    print(f"DEL matching rate: {del_match_rate:.2%}")
    print("==========================================")


if __name__ == "__main__":
    main()
