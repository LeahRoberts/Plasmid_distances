"""
script to parse nucmer alignment results and determine pairwise
plasmid cluster identity

Requires:
1. list of plasmid lengths
2. plasmid fasta files

Dependencies: 
1. mummer (nucmer)
2. biopython

Script will only accept alignment regions >1000bp and >90% identity

IMPORTANT: put the merge_fasta.sh script somewhere and point to it within this script (line 85)

To run: ebenn_plasmids.py <directory with plasmids in fasta> <list of plasmids and their length>
"""

import glob
import os
import subprocess
import sys
from subprocess import CalledProcessError
from statistics import mean
from Bio import SeqIO

input_dir = sys.argv[1]
lengths = sys.argv[2]

plasmid_lengths = {}
plasmid_clusters = {}
cluster_alignments = {}


def sort_plasmids_into_clusters(plasmid):
    clust = plasmid.split("_")[2].rstrip()
    clust = clust.replace(".fasta", "")
    if clust not in plasmid_clusters.keys():
        plasmid_clusters[clust] = []
    plasmid_clusters[clust].append(plasmid)


def bashcommand(command):
    with open("command_log.txt", "a") as fout:
        fout.write("%s\n" % command)
    try:
        subprocess.check_output(command, stderr=subprocess.STDOUT, shell=True)
    except CalledProcessError as e:
        print(str(e.output))
        exit(1)


def run_nucmer(plasmidA, plasmidB, path):
    A = path + "/" + plasmidA
    B = path + "/" + plasmidB
    bashcommand("nucmer --maxmatch --minmatch=1000 %s %s" % (A, B)) # --minmatch controls min length match
    bashcommand("show-coords -B -I 90 -r out.delta > out.coords") # -I controls min identity filter


# there might be a slight problem with this if there are repeats that appear as multiple separate blocks
# which will inflate the total alignment
def parse_nucmer(nucmer_output_file):
    with open(nucmer_output_file, "r") as align_in:
        total_align = 0
        for block in align_in:
            pos1 = block.split("\t")[6]
            pos2 = block.split("\t")[7]
            block_length = abs(float(pos1) - float(pos2))
            total_align += block_length
    return total_align


all_plasmids = glob.glob(input_dir + "/*.fasta")

# merge multi-contig fasta files
for p_check in all_plasmids:
    with open(p_check) as handle:
        record_count = 0
        for record in SeqIO.parse(handle, "fasta"):
            record_count += 1
        if record_count > 1:
            p = os.path.splitext(p_check)[0]+'.merged.fasta'
            #### CHANGE THE PATH TO MERGE_FASTA.SH HERE ####
            bashcommand("/path/to/merge_fasta.sh %s %s" % (p_check, p)) 
            all_plasmids = list(map(lambda x: x.replace(p_check, p), all_plasmids))
        else:
            p = p_check
    p = os.path.split(p)[1]
    sort_plasmids_into_clusters(p)

with open(lengths, "r") as fin:
    for line in fin:
        plas_name = line.split()[0]
        bps = line.split()[1].rstrip()
        plasmid_lengths[plas_name] = bps

# go through all clusters and do pairwise alignments
with open("plasmid_alignments_summary.tsv", "w") as sum_out:
    sum_out.write("plas1\tplas2\tplas1_align\tplas2_align\n")
for c in plasmid_clusters.keys():
    plasmid_list = plasmid_clusters[c]
    compared = []
    cluster_alignments[c] = []
    for pl1 in plasmid_list:
        if pl1 not in compared:
            compared.append(pl1)
            for pl2 in plasmid_list:
                if pl1 != pl2:
                    run_nucmer(pl1, pl2, input_dir)
                    alignment = parse_nucmer("out.coords")
                    try:
                        pl1_alignment_proportion = float(alignment) / float(plasmid_lengths[pl1])
                    except KeyError:
                        pl1_name = os.path.splitext(pl1)[0].replace("merged", "")
                        pl1_name = pl1_name + "fasta"
                        pl1_alignment_proportion = float(alignment) / float(plasmid_lengths[pl1_name])
                    try:
                        pl2_alignment_proportion = float(alignment) / float(plasmid_lengths[pl2])
                    except KeyError:
                        pl2_name = os.path.splitext(pl2)[0].replace("merged", "")
                        pl2_name = pl2_name + "fasta"
                        pl2_alignment_proportion = float(alignment) / float(plasmid_lengths[pl2_name])
                    cluster_alignments[c].append(pl1_alignment_proportion)
                    cluster_alignments[c].append(pl2_alignment_proportion)
                    with open("plasmid_alignments_summary.tsv", "a") as result_out:
                        result_out.write("%s\t%s\t%s\t%s\n" % (pl1, pl2, pl1_alignment_proportion, pl2_alignment_proportion))
                    os.remove("out.delta")
                    os.remove("out.coords")
    if len(cluster_alignments[c]) > 1:
        avg = mean(cluster_alignments[c])
    else:
        print(c, "has only one representative...skipping")
    with open("cluster_alignment_summary.tsv", "a") as clust_out:
        clust_out.write("%s\t%s\n" % (c, avg))
