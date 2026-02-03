# benchmark example of PSVGT
## 01 simulation data benckmark
```
├── 01_simualted_data
│   ├── bench2vcf.py 
│   ├── benchmarking4suppTable.py 
│   ├── filelist 
│   ├── sim_find_pos.py
│   ├── simulate_sv_genome.py 
│   └── TRA
│       ├── 0.get_genome.sh 
│       ├── tra_bench_mark_GT_all_hetero.py 
│       ├── tra_bench_mark_GT_all_homo.py
│       └── tra_noY.bed 
```
### genome simulated and variants table (allsv-random-sorted-rev) generated
```shell
python3 simulate_sv_genome.py 
```
***diploid-resovled genome will be generated with SV(10000 DEL, 10000 INS, 1000 DUP and 1000 INV)***

### translocation variats genome
```shell
bash 0.get_genome.sh 
```
***tra_noY.bed is translocation variants table collected from cuteSV methods***

### To benchmark the performance of each method in simulated datasets
```
python bench2vcf.py allsv-random-sorted-rev allsv-random-sorted-rev.vcf
python benchmarking4suppTable.py  script.py allsv-random-sorted-rev.vcf method_result.vcf
```

## 02 real data benchmark 
```
├── 02_real_data
│   ├── 20bp_breakpoint_plot4truvari_output_v2.py 
│   └── f1_for_size_range.py
```
### To assess the breakpoint accuracy
```
python 20bp_breakpoint_plot4truvari_output_v2.py truvari_bench/tp-comp.vcf breakpoint.pdf
```
***truvari_bench folder is from truvari program, the tp-com.vcf file can be use to assess the breakpoint accuracy.***

### To assess different size SV accuracy
```
python f1_for_size_range.py -h
usage: f1_for_size_range.py [-h] -t TP -f FP -n FN [-o OUTPUT]

Calculate F1 scores for SV size ranges, the input file is from truvari output

options:
  -h, --help            show this help message and exit
  -t TP, --tp TP        True positive VCF file
  -f FP, --fp FP        False positive VCF file
  -n FN, --fn FN        False negative VCF file
  -o OUTPUT, --output OUTPUT
                        Output CSV file

python f1_for_size_range.py -t truvari_bench/tp-base.vcf -f truvari_bench/fp.vcf -n truvari_bench/fn.vcf -o size_range_f1.score.txt
```
***the reqiured input can also be found from travari output folder***

## benchmark example of plant data
```
├── 03_plants_SV
│   ├── Arabidopsis
│   │   ├── Arabidopsis.bed
│   │   ├── benchmark_delly_manta_PSVGT_no_bed.py
│   │   ├── benchmark_delly_manta_PSVGT_with_bed.py
│   │   ├── chr1_Arabidopsis4genotype_concrodance_hifi.vcf
│   │   ├── chr1_Arabidopsis4reacall_precision_f1.cr.vcf
│   │   ├── chr1_Arabidopsis_short_read_assemblies.PSVGT_lrGT.vcf
│   │   ├── chr1_Arabidopsis_sorted.bam_bpgenotype.PSVGT_srGT.vcf
│   │   ├── chr1_Arabidopsis_sorted.bam.delly.vcf
│   │   ├── chr1_Arabidopsis_sorted.bam.manta.vcf
│   │   ├── comparison_result_fn.list
│   │   └── example_bench.sh
│   ├── potato
│   │   ├── chr1_debreak.vcf
│   │   ├── chr1_msv_no_suppp
│   │   ├── chr1_PSVGT_signal_falseGT.vcf
│   │   ├── chr1_sniffles.vcf
│   │   ├── genome_msv_no_supp
│   │   ├── msv.sh
│   │   └── msv_Stat4_Table.py
│   └── rice
│       ├── benchmarkSVInDEL4rice.py
│       ├── chr1_GP58_cr_benchmark.vcf
│       ├── chr1_GP58_sr.bam_delly.vcf.recode.vcf
│       ├── chr1_GP58_sr.bam_manta.vcf.recode.vcf
│       ├── chr1_GP58_sr.bam_PSVGT.txt
│       ├── example.benchmark.sh
│       └── GP58_5k.bed 
```
***The usage and benchmark command can be found in example_xx.sh***

