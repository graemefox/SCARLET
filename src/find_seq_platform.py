#!/usr/bin/python3

import argparse
import pysam

parser = argparse.ArgumentParser(description='Arguments for find_seq_platform.py')
parser.add_argument('-i','--input_bam', \
                    help='The input bam from which to find the source seq platform', \
                    required=True)
args = parser.parse_args()

input_bam = pysam.AlignmentFile(args.input_bam, "rb")

prom_flowcell_positions = [
    "1A", "2A", "3A", "4A", "5A", "6A", \
    "1B", "2B", "3B", "4B", "5B", "6B", \
    "1C", "2C", "3C", "4C", "5C", "6C", \
    "1D", "2D", "3D", "4D", "5D", "6D", \
    "1E", "2E", "3E", "4E", "5E", "6E", \
    "1F", "2F", "3F", "4F", "5F", "6F", \
    "1G", "2G", "3G", "4G", "5G", "6G", \
    "1H", "2H", "3H", "4H", "5H", "6H"]


# read header
header = input_bam.header.to_dict()
seq = ""
# look for RG tag
for key, value in header.items():
    if key == "RG":
        for item in value:
            PM_code = item['PM'][0:3]
            PL_code = item['PL']
        if PM_code:
            if PM_code == "GXB":
                seq = "GridION"
            if PM_code == "P2S":
                seq = "P2Solo"
            if PM_code in prom_flowcell_positions:
                seq = "PromethION"

for i, alignment in enumerate(input_bam):
    # if it didn't find the RG tag and platform, try the fn and f5 tags
    if seq == "":
        data = str(alignment).split("\t")[11]
        if "\'fn\'" in data:
            if data.split("\'fn\', \'")[1][0] == "F":
                seq = "MinION/GridION"
            if data.split("\'fn\', \'")[1][0] == "P":
                seq = "PromethION"
        if "\'f5\'" in data:
            if data.split("\'f5\', \'")[1][0] == "F":

                # check for contradiction
                if seq == "" or seq == "MinION/GridION":
                    seq = "MinION/GridION"
                elif seq != "" and seq != "MinION/GridION":
                    seq = "Unknown"

            if data.split("\'f5\', \'")[1][0] == "P":
                # TODO - contradcition check below needs work
                if seq == "" or seq == "PromethION":
                    seq = "PromethION"
                elif seq != "" and seq != "PromethION":
                    seq = "Unknown"
        # break after the first line of the file
        if i >= 0:
            break
    else:
        break

if seq == "":
    seq = "Unknown"
print(seq.rstrip("\n"))

input_bam.close()
