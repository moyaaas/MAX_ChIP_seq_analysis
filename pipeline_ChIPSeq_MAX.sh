#!/bin/bash

# ==============================================================================
# HOW TO USE THIS SCRIPT:
# 1. Make the file executable by running this command in your terminal:
#    chmod +x pipeline_ChIPSeq_MAX.sh
#  
# 2. Run the script by typing:
#    ./pipeline_ChIPSeq_MAX.sh
# if Linux doesn't find the file, try
#   dos2unix pipeline_ChIPSeq_MAX.sh
# ==============================================================================

# Exit immediately if any command fails, and print commands as they are executed
set -e 
set -x 

echo "Starting end-to-end MAX pipeline..."

# Define input raw file names
REP1_RAW="rep1_unfiltered.bam"
REP2_RAW="rep2_unfiltered.bam"
CONTROL_RAW="control_unfiltered.bam"
GOLD_STANDARD="starpeaks.bed"
BLACKLIST="GRCh38_unified_blacklist.bed"
CHROMHMM="K562_ChromHMM_15states.bed"

# ==========================================
# 0. Filtering Multi-Mapping Reads (samtools)
# ==========================================
echo "### STEP 0: mapping QC + filtering multi-mappingreads ###"

samtools flagstat $REP1_RAW    > rep1_unfiltered_flagstat.txt
samtools flagstat $REP2_RAW    > rep2_unfiltered_flagstat.txt
samtools flagstat $CONTROL_RAW > control_unfiltered_flagstat.txt

# Filter out multi-mapping reads (MAPQ = 0) to generate unique BAM files
samtools view -b -q 1 $REP1_RAW > rep1_filtered_unique.bam
samtools view -b -q 1 $REP2_RAW > rep2_filtered_unique.bam
samtools view -b -q 1 $CONTROL_RAW > control_filtered_unique.bam

samtools flagstat rep1_filtered_unique.bam    > rep1_filtered_flagstat.txt
samtools flagstat rep2_filtered_unique.bam    > rep2_filtered_flagstat.txt
samtools flagstat control_filtered_unique.bam > control_filtered_flagstat.txt

# Report QC metrics: % mapped e % multi-mapping (mapped_unfiltered - mapped_unique)
report_qc () {
  local name=$1 raw_stat=$2 filt_stat=$3
  local total mapped_raw mapped_uniq multi pct_mapped pct_multi
  total=$(grep "in total" "$raw_stat" | cut -d' ' -f1)
  mapped_raw=$(grep " mapped (" "$raw_stat" | head -1 | cut -d' ' -f1)
  mapped_uniq=$(grep " mapped (" "$filt_stat" | head -1 | cut -d' ' -f1)
  multi=$((mapped_raw - mapped_uniq))
  pct_mapped=$(awk -v a="$mapped_raw" -v b="$total" 'BEGIN{printf "%.2f", a*100/b}')
  pct_multi=$(awk -v a="$multi" -v b="$total" 'BEGIN{printf "%.2f", a*100/b}')
  echo "${name}: total=${total} mapped=${mapped_raw} (${pct_mapped}%) unique=${mapped_uniq} multimapping=${multi} (${pct_multi}%)"
}

{
  report_qc REP1    rep1_unfiltered_flagstat.txt    rep1_filtered_flagstat.txt
  report_qc REP2    rep2_unfiltered_flagstat.txt    rep2_filtered_flagstat.txt
  report_qc CONTROL control_unfiltered_flagstat.txt control_filtered_flagstat.txt
} | tee -a qc_summary.txt

# ==========================================
# 1. Peak Calling (MACS2)
# ==========================================
echo "### STEP 1: peak calling ###"

# Merge filtered BAM files to simulate a higher-depth pooled experiment
samtools merge -f pooled_filtered_unique.bam rep1_filtered_unique.bam rep2_filtered_unique.bam

# Call narrow peaks and generate bedGraph files (-B --SPMR) for UCSC visualization
# This is a CRISPR-edited dataset, so we use a q-value cutoff of 0.01 for peak calling, instead of the default 0.05
macs2 callpeak -t rep1_filtered_unique.bam -c control_filtered_unique.bam -g hs -n rep1 -q 0.01 -B --SPMR
macs2 callpeak -t rep2_filtered_unique.bam -c control_filtered_unique.bam -g hs -n rep2 -q 0.01 -B --SPMR
macs2 callpeak -t pooled_filtered_unique.bam -c control_filtered_unique.bam -g hs -n pooled -q 0.01 -B --SPMR

# Generate PDF model plots using the R scripts created by MACS2
R < rep1_model.r --vanilla
R < rep2_model.r --vanilla
R < pooled_model.r --vanilla

# QC: redundant rate and fragment size estimation 
{
  for n in rep1 rep2 pooled; do
    echo "--- $n ---"
    grep -E "Redundant rate|# d =|tags after filtering" ${n}_peaks.xls
  done
} | tee -a qc_summary.txt

# plot estimation fragment size (positive/negative distance + cross-correlation)
for n in rep1 rep2 pooled; do
  if [ -f ${n}_model.r ]; then
    R < ${n}_model.r --vanilla
  fi
done

# ==========================================
# 2. Intersection and Sorting
# ==========================================
echo "### STEP 2: sorting ###"

# bedtools jaccard requires all files to be strictly sorted by chromosome and start position
sort -k1,1 -k2,2n rep1_peaks.narrowPeak > rep1_sorted.bed
sort -k1,1 -k2,2n rep2_peaks.narrowPeak > rep2_sorted.bed
sort -k1,1 -k2,2n pooled_peaks.narrowPeak > pooled_sorted.bed
sort -k1,1 -k2,2n $GOLD_STANDARD > starpeaks_sorted.bed

# Comparison REP vs REP (Challenge 1)
#    -> two criteria: overlap peaks + summit proximity (100bp)
# ==========================================
echo "### STEP 3: comparison rep1 vs rep2 ###"
# Extract peaks that physically overlap between the two replicates (-u no duplicates)
bedtools intersect -a rep1_peaks.narrowPeak -b rep2_peaks.narrowPeak -u > overlaps.bed
sort -k1,1 -k2,2n overlaps.bed > overlaps_sorted.bed

n_rep1=$(wc -l < rep1_sorted.bed)
n_rep2=$(wc -l < rep2_sorted.bed)
n_overlap=$(wc -l < overlaps_sorted.bed)
echo "Overlap regions: ${n_overlap} peaks rep1 (out of ${n_rep1}) overlap with a peak in rep2 (total rep2=${n_rep2})" | tee -a qc_summary.txt

# 3b. summit 100bp (file *_summits.bed generated directly by MACS2)
bedtools window -w 100 -a rep1_summits.bed -b rep2_summits.bed -u > rep1_summits_near_rep2.bed
bedtools window -w 100 -a rep2_summits.bed -b rep1_summits.bed -u > rep2_summits_near_rep1.bed

n_s1=$(wc -l < rep1_summits.bed)
n_s2=$(wc -l < rep2_summits.bed)
n_close1=$(wc -l < rep1_summits_near_rep2.bed)
n_close2=$(wc -l < rep2_summits_near_rep1.bed)
echo "Summit proximity (<=100bp): ${n_close1}/${n_s1} summit rep1 close to summit rep2" | tee -a qc_summary.txt
echo "Summit proximity (<=100bp): ${n_close2}/${n_s2} summit rep2 close to summit rep1" | tee -a qc_summary.txt


# ==========================================
# 3. ENCODE comparison: overlap % + Jaccard Index Calculation
# ==========================================
echo "### STEP 4: comparison with ENCODE (gold standard) ###"

count_overlap () { bedtools intersect -a "$1" -b starpeaks_sorted.bed -u | wc -l; }

for f in rep1_sorted.bed rep2_sorted.bed overlaps_sorted.bed pooled_sorted.bed; do
  n=$(wc -l < "$f")
  ov=$(count_overlap "$f")
  echo "$f: ${ov}/${n} peaks overlap with a peak in the ENCODE Gold Standard" | tee -a qc_summary.txt
done

# Compare your peak sets against the ENCODE Gold Standard
bedtools jaccard -a rep1_sorted.bed -b starpeaks_sorted.bed > jaccard_rep1.txt
bedtools jaccard -a rep2_sorted.bed -b starpeaks_sorted.bed > jaccard_rep2.txt
bedtools jaccard -a overlaps_sorted.bed -b starpeaks_sorted.bed > jaccard_intersect.txt
bedtools jaccard -a pooled_sorted.bed -b starpeaks_sorted.bed > jaccard_pooled.txt

get_jaccard () { tail -n1 "$1" | cut -f3; }

{
  echo -e "rep1_sorted.bed\t$(get_jaccard jaccard_rep1.txt)"
  echo -e "rep2_sorted.bed\t$(get_jaccard jaccard_rep2.txt)"
  echo -e "overlaps_sorted.bed\t$(get_jaccard jaccard_intersect.txt)"
  echo -e "pooled_sorted.bed\t$(get_jaccard jaccard_pooled.txt)"
} > jaccard_table.txt

cat jaccard_table.txt | tee -a qc_summary.txt

# automatic choice of the set with the highest Jaccard 
FINAL_PEAKS=$(grep -E "overlaps_sorted.bed|pooled_sorted.bed" jaccard_table.txt | sort -k2,2gr | head -n1 | cut -f1)
echo "==> Set of final peaks chosen (highest Jaccard between intersect and pooled): ${FINAL_PEAKS}" | tee -a qc_summary.txt

# ==========================================
# 4. Blacklist Filtering 
# ==========================================
echo "### STEP 5: removal of blacklist regions ###"

n_before=$(wc -l < "$FINAL_PEAKS")
bedtools intersect -a "$FINAL_PEAKS" -b $BLACKLIST -v > final_peaks_clean.bed
n_after=$(wc -l < final_peaks_clean.bed)
echo "Blacklist: removed $((n_before-n_after)) peaks out of ${n_before} (remaining ${n_after})" | tee -a qc_summary.txt


# ==========================================
# 5. Split shared vs unique (vs ENCODE) + q-value profiling
# ==========================================
# Split peaks into two groups: shared with ENCODE vs. unique to your pipeline
bedtools intersect -a final_peaks_clean.bed -b starpeaks_sorted.bed -u > shared_peaks.bed
bedtools intersect -a final_peaks_clean.bed -b starpeaks_sorted.bed -v > unique_peaks.bed

# Extract the -log10(q-value) from column 9 for boxplot generation in R
awk '{print $9}' shared_peaks.bed > qval_shared.txt
awk '{print $9}' unique_peaks.bed > qval_unique.txt

# ==========================================
# 6. Summit Extraction & ChromHMM Annotation
# ==========================================
echo "### STEP 7: final summit extraction and chromatin state annotation ###"

# Calculate absolute summit coordinates (start + offset) and create 1-bp regions
awk 'BEGIN{OFS="\t"} {summit=$2+$10; print $1, summit, summit+1, $4, $9}' final_peaks_clean.bed > final_summits.bed

# Extract the top 5000 most significant summits to avoid GREAT over-saturation
sort -k9,9nr final_peaks_clean.bed | head -n 5000 | \
    awk 'BEGIN{OFS="\t"} {summit=$2+$10; print $1, summit, summit+1, $4}' > top5000_summits.bed

# Overlap summits with the 15-state ChromHMM annotations
bedtools intersect -a final_summits.bed -b $CHROMHMM -wa -wb > summits_chromhmm.txt
# Summarize the count of summits falling into each chromatin state (column 9)
awk '{print $9}' summits_chromhmm.txt | sort | uniq -c | sort -nr > chromhmm_summary.txt

echo "Distribuzione stati chromHMM (grezza, per singolo stato):" | tee -a qc_summary.txt
cat chromhmm_summary.txt | tee -a qc_summary.txt


awk '{
  state=$9;
  if (state ~ /Tss/)                 group="Promoter";
  else if (state ~ /Enh/)             group="Enhancer";
  else if (state ~ /Tx/)              group="Transcribed";
  else if (state ~ /Repr|Het|Quies/)  group="Repressed/Quiescent";
  else                                 group="Other (ZNF/Biv/...)";
  print group
}' summits_chromhmm.txt | sort | uniq -c | sort -nr > chromhmm_major_categories.txt

echo "Distribution per macro categories (verify/tweak the regex above!):" | tee -a qc_summary.txt
cat chromhmm_major_categories.txt | tee -a qc_summary.txt

echo ""
echo "Pipeline completed. See qc_summary.txt for a summary of QC metrics and results."
