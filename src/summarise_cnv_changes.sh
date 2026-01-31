#!/bin/bash

# Usage: ./summarize_cn_changes_with_cytoband.sh input.seg output.csv centromeres.txt

seg_file="$1"
output_file="$2"
cyto_file="$3"

if [[ -z "$seg_file" || -z "$output_file" || -z "$cyto_file" ]]; then
  echo "Usage: $0 <input_seg_file> <output_csv_file> <centromere_file>"
  exit 1
fi

# Write header to output
#echo "CHROMOSOME,ARM,STATUS,NOTES" > "$output_file"

# Preprocess centromere positions: get lowest and highest position of acen per chromosome
awk 'BEGIN { OFS="\t" }
    $5 == "acen" {
        gsub("chr", "", $1)
        chr = $1
        if (chr_cent_start[chr] == "" || $2 < chr_cent_start[chr]) chr_cent_start[chr] = $2
        if (chr_cent_end[chr] == "" || $3 > chr_cent_end[chr]) chr_cent_end[chr] = $3
    }
    END {
        for (chr in chr_cent_start)
            print chr, chr_cent_start[chr], chr_cent_end[chr]
    }
' "$cyto_file" > .centromere_ranges.tsv

# Main processing
awk -F'\t' -v OFS=',' -v CENTROMERE_FILE=".centromere_ranges.tsv" '
BEGIN {
    # Read centromere regions
    while ((getline < CENTROMERE_FILE) > 0) {
        chr = $1
        cent_start[chr] = $2
        cent_end[chr] = $3
    }
}
NR > 1 {
    chrom = $2
    gsub(/^chr/, "", chrom)
    start = $3
    stop = $4
    ratio = $6 + 0

    if (chrom != "X" && chrom != "Y") {

    if (ratio < -0.5 || ratio > 0.5) {
        status = (ratio > 0.5) ? "Gain" : "Loss"

        cent_s = cent_start[chrom]
        cent_e = cent_end[chrom]

        if (cent_s == "" || cent_e == "") {
            arm = "unknown"
            note = "No centromere data"
        } else if (stop < cent_s) {
            arm = "p"
            note = "log2=" ratio
        } else if (start > cent_e) {
            arm = "q"
            note = "log2=" ratio
        } else {
            arm = "p+q"
            note = "Spans centromere log2=" ratio
        }
        print "CNV," chrom " " arm " " status " " note
    }
  }
}
' "$seg_file" >> "$output_file"

# Clean up
rm -f .centromere_ranges.tsv
