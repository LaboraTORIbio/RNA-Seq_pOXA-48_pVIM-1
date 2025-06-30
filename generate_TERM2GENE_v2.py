#!/usr/bin/env python3

import argparse
import pandas as pd


parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                 description="""Generates the TERM2GENE table for clusterProfiler"""
                                )

parser.add_argument("-g", "--go_list", help="List of unique GO terms", type=str, dest="go_list", required=True)
parser.add_argument("-a", "--annot", help="Table of GO term annotations retrived from UniProt", type=str, dest="uniref2annot", required=True)
parser.add_argument("-m", "--gene2uniref", help="Mapping of GeneIDs (1st column) and UniRef90 IDs (2nd column)", type=str, dest="gene2uniref", required=True)
parser.add_argument("-o", "--term2gene", help="Path and name for the TERM2GENE file", type=str, dest="term2gene", required=True)

args = parser.parse_args()

go_list = args.go_list
uniref2annot = args.uniref2annot
gene2uniref = args.gene2uniref
term2gene = args.term2gene


### Store the GO terms in a list:

with open(go_list, "r") as file:
    go_terms = [line.strip() for line in file]

### Create a dataframe with GO term and UniRef90 mappings:

uniref2term_df = pd.read_csv(uniref2annot, sep="\t")

df_exploded = uniref2term_df[["Entry", "Gene Ontology (GO)"]].copy()
df_exploded["GO_annot_list"] = df_exploded["Gene Ontology (GO)"].str.split("; ")
df_exploded = df_exploded.explode("GO_annot_list")
df_exploded["GO_ID"] = df_exploded["GO_annot_list"].str.extract(r"\[(GO:\d+)\]")

term2uniref_df = df_exploded[df_exploded["GO_ID"].isin(go_terms)][["GO_ID", "Entry"]].reset_index(drop=True)
term2uniref_df = term2uniref_df.rename(columns={"GO_ID": "GO_terms", "Entry": "UniRef90"})

### Read table with Gene ID and UniRef90 mappings:

gene2uniref_df = pd.read_csv(gene2uniref, sep="\t")

### Create TERM2GENE table (GO ID - Gene ID):

df_merged = pd.merge(term2uniref_df, gene2uniref_df, on="UniRef90")
df_merged = df_merged.drop("UniRef90", axis=1)
df_merged.sort_values("GO_terms")

### Export table:

df_merged.to_csv(term2gene, sep="\t", index=False, header=False)
print(f"TERM2GENE table generated in {term2gene}")
