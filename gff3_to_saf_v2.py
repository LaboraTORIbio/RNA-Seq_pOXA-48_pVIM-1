#!/usr/bin/env python3

import argparse
import pandas as pd


parser = argparse.ArgumentParser(description="""GFF3 to SAF""")
parser.add_argument("--gff", help="Path to input GFF3 file", type=str, dest="gff", required=True)
parser.add_argument("--saf", help="Path to output SAF file", type=str, dest="saf", required=True)
args = parser.parse_args()
gff_file = args.gff
saf_file = args.saf

gff_colnames = ["Chr", "Source", "Type", "Start", "End", "Score", "Strand", "Phase", "Attributes"]
types_to_keep = ["CDS", "tRNA", "ncRNA", "tmRNA", "antisense_RNA", "RNase_P_RNA", "SRP_RNA"]
cols_for_saf = ["GeneID", "Chr", "Start", "End", "Strand", "Type", "UniRef90", "Gene", "Product"]

gff_df = pd.read_csv(gff_file, sep="\t", names=gff_colnames, usecols=["Chr", "Type", "Start", "End", "Strand", "Attributes"])
gff_df = gff_df[gff_df["Chr"].str.contains("#")==False].reset_index(drop=True)
gff_df = gff_df[gff_df["Type"].isin(types_to_keep)==True].reset_index(drop=True)
gff_df[["Start", "End"]] = gff_df[["Start", "End"]].astype(int)

gff_df["GeneID"] = gff_df["Attributes"].str.extract(r"ID=([^;]+)")
gff_df["Product"] = gff_df["Attributes"].str.extract(r"product=([^;]+)")
gff_df["Gene"] = gff_df["Attributes"].str.extract(r"gene=([^;]+)")
gff_df["UniRef90"] = gff_df["Attributes"].str.extract(r"UniRef90_([^,;]+)")

saf_df = gff_df[cols_for_saf]
saf_df = saf_df.fillna("-")

saf_df.to_csv(saf_file, sep="\t", index=False)
