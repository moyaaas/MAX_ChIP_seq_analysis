# ChIP-Seq Analysis pipeline (MAX)

End-to-End Bash pipeline for ChIP-Seq peak calling and characterization of the **MAX** transcription factor, K562 cell line and GRCh38. 
Starting from unfiltered BAM files retrieved from the ENCODE database, two replicates and a control, it filters multi-mapping read, calls peaks with MACS2 compares replicates, benchmarks the ENCODE gold standard, removes blacklisted regions and annotates the final summits with ChromHMM chromatin states.

## Requirements
- `samtools`
- `MACS2`
- `bedtools`
- `R` (used by MACS2 model scripts; boxplots of q-values are produced downstream from `qval_*.txt`)

## Input files
| Variable | File | Description |
|---|---|---|
| `REP1_RAW` | `rep1_unfiltered.bam` | Replicate 1 alignment |
| `REP2_RAW` | `rep2_unfiltered.bam` | Replicate 2 alignment |
| `CONTROL_RAW` | `control_unfiltered.bam` | Input/Control alignment |
| `GOLD_STANDARD` | `starpeaks.bed` | ENCODE gold standard peaks |
| `BLACKLIST` | `GRCh38_unified:blacklist.bed` | Blacklist Regions |
| `CHROMHMM` | `K562_ChromHMM_15states.bed` | ChromHMM 15-state annotation for K562

## Usage
1. Make the file executable by running this command in your terminal:
    `chmod +x pipeline_ChIPSeq_MAX.sh`
2. Run the script by typing:
    `./pipeline_ChIPSeq_MAX.sh`

If Linux cannot run the file (Windows line endings), run `dos2unix pipeline_ChIPSeq_MAX.sh` first. The script uses `set -e` (stops at the first error) and `set -x` (prints each command).
Place all input BAM files in the same directory as the script.

## Pipeline Steps
0. **Mapping QC and filtering**
1. **Peak Calling (MACS2)**
2. **Sorting**
3. **Replicate Comparison**
4. **ENCODE Comparison**
5. **Blacklist Filtering**
6. **Shared vs Unique Peaks**
7. **Summits and ChromHMM annotation**

## Output files
| File | Content |
|---|---|
| `qc_summary.txt` | Summary of all QC metrics and results |
| `*_flagstat.txt` | Mapping statistics before/after filtering |
| `rep1/rep2/pooled_peaks.narrowPeak`, `*_summits.bed`, `*.bdg`, `*_model.pdf` | MACS2 outputs |
| `overlaps_sorted.bed` | Rep1 peaks overlapping rep2 |
| `rep*_summits_near_rep*.bed` | Summits within 100 bp between replicates |
| `jaccard_*.txt`, `jaccard_table.txt` | Jaccard indices vs ENCODE |
| `final_peaks_clean.bed` | Final peaks after blacklist removal |
| `shared_peaks.bed`, `unique_peaks.bed` | Peaks shared with / unique vs ENCODE |
| `qval_shared.txt`, `qval_unique.txt` | -log10(q-value) for boxplots |
| `final_summits.bed`, `top5000_summits.bed` | Summits (all / top 5,000, e.g. for GREAT) |
| `summits_chromhmm.txt`, `chromhmm_summary.txt`, `chromhmm_major_categories.txt` | ChromHMM annotation and counts |

## Notes
The q-value cutoff of 0.01 was chosen for this CRISP-edited dataset, that is stricter than the MACS2 default of 0.05.
The macro-category of chromatin annotation must be checked against the state labels of ChromHMM file.

## Project Presentation
For an overview of the biological context, quality control steps, and final biological results (including ChromHMM state distributions), please check the project presentation: 
**[View the Presentation PDF](./MAX_ChIP_Seq_Presentation.pdf)**

## License
MIT License © 2026 Martina Montonati