#!/usr/bin/env python3

import pandas as pd
import os, argparse
from functools import reduce


def table_to_dict(directory, file):
    """
    Reads a tsv file containing differential expression results into a pandas dataframe
    and returns a dictionary with the dataframe as value and the strain name as key.
    """
    strain_name = directory.replace("RNAseq_", "")
    df = pd.read_csv(f"{directory}/{file}", sep="\t")
    return {strain_name : df}


def prepare_dfs_in_dict_for_merging(input_dict):
    """
    Prepares the dataframes stored in a dictionary for posterior merging:
    Performs a filtering step to keep only DEGs that are CDSs, excluding hypothetical proteins,
    and, if present, filters for 'pOXA-VIM' replicons, adding an empty GO column for compatibility.
    Additionally, renames log2FC columns to the strain name for merging.
    """
    out_dict = {}
    for strain, df in input_dict.items():
        if "pOXA-VIM" in df["Chr"].values:
            df = df.loc[df["Chr"] == "pOXA-VIM"].copy()
            df = df.loc[df["Product"] != "IS4 family IS10A transposase"]
            df["GO"] = "-"
        df = df[["Chr", "Type", "UniRef90", "Gene", "Product", "GO", "log2FoldChange"]]
        df = df.loc[(df["Type"] == "CDS") & (df["Product"] != "hypothetical protein")]
        df = df.reset_index(drop=True)
        df = df.rename({"log2FoldChange": strain}, axis=1)
        df = df.drop(columns=["Chr"])
        out_dict.update({strain : df})
    return out_dict


def merge_dfs_in_dict_from_a_comparative_analysis(dictionary):
    """
    Merges the dataframes stored in a dictionary,
    all corresponding to a specific comparative analysis,
    by columns "Type", "UniRef90", "Gene", "Product" and "GO",
    and adds a column with the count of shared DEGs across strains.
    """
    merged_df = reduce(
        lambda left, right: pd.merge(
              left,
              right,
              on=["Type", "UniRef90", "Gene", "Product", "GO"],
               how="outer"
           ),
           dictionary.values()
    ).fillna("0")
    logfc_cols = [col for col in merged_df.columns if col not in ["Type", "UniRef90", "Gene", "Product", "GO"]]
    merged_df[logfc_cols] = merged_df[logfc_cols].astype(float).round(2)
    merged_df["COUNT"] = (merged_df[logfc_cols] != 0).sum(axis=1)
    merged_df["GO"] = merged_df.pop("GO")
    merged_df["GO"] = merged_df["GO"].replace("0", "-")
    return merged_df


def dict_to_excel(dictionary, comploc, writer):
    """
    Exports the dataframes stored in a dictionary to a Excel, including the strain names
    and the comparison and genomic location (chromosome or plasmids) in the sheet names.
    """
    for strain, df in dictionary.items():
        sheet_name = f"{strain}_{comploc}"
        if len(sheet_name) > 31:
            sheet_name = sheet_name[:31]
        df.to_excel(writer, sheet_name=sheet_name, index=False)


def apply_bold_format_to_iron_related_degs_and_store(df, sheet_name, writer):
    """
    Searches for DEGs involved in iron-related processes, excluding proteins forming clusters with iron or
    involved in the electron transport chain, formats these DEGs in bold, and returns these rows in a list. 
    """
    wanted_terms = ["Fe", "iron", "enterobactin", "ferrous", "heme", "siderophore", "ferritin"]
    exclude_terms = ["cluster", "cytochrome", "electron"]
    worksheet = writer.book.get_worksheet_by_name(sheet_name)
    highlighted_rows = []
    if worksheet:
        for idx, row in df.iterrows():
            row_text = " ".join(str(x) for x in row.astype(str).values)
            if any(term in row_text for term in wanted_terms) and all(term not in row_text for term in exclude_terms):
                worksheet.set_row(idx + 1, None, bold_format)
                row_with_source = row.copy()
                row_with_source["SourceSheet"] = sheet_name
                highlighted_rows.append(row_with_source)
    return highlighted_rows


if __name__ == "__main__":
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                    description="""
                                    DEGs TABLE PARSER

                                    This script reads the DE_results_*_filtered_padj_chromosome_GOannot.tsv
                                    and DE_results_*_filtered_padj_plasmids.tsv files located in the provided directories,
                                    outputs them in a single Excel in different sheets,
                                    and parses the DE_results_*_filtered_padj_chromosome_GOannot.tsv files
                                    to find common chromosomal DEGs between the strains per comparison,
                                    outputting the summary table in a different sheet per comparison.

                                    The script also outputs a summary table of the DEGs of pOXA-VIM across strains,
                                    and of genes involved in iron-related processes.
                                    """
                                    )

    parser.add_argument("-d", nargs="+", help="Directories to parse", type=str, dest="dir_names", required=True)
    parser.add_argument("-n", help="Base name for the output file", type=str, dest="out_name", required=True)

    args = parser.parse_args()
    dir_names = args.dir_names
    dir_names = [i.replace("/", "") for i in dir_names]
    out_name = args.out_name

    ### Reading files:
    # The name must be [comp_separated_by_-]_[genomic_location]. The chromosome files must end with _chr and plasmid with _plas

    files = {
        "pOXA-48-vs-PF_chr": "DE_results_pOXA-48-vs-PF_filtered_padj_chromosome_GOannot.tsv",
        "pVIM-1-vs-PF_chr": "DE_results_pVIM-1-vs-PF_filtered_padj_chromosome_GOannot.tsv",
        "pVIM-1-vs-pOXA-48_chr": "DE_results_pVIM-1-vs-pOXA-48_filtered_padj_chromosome_GOannot.tsv",
        "pOXA-48-vs-PF_plas": "DE_results_pOXA-48-vs-PF_filtered_padj_plasmids.tsv",
        "pVIM-1-vs-PF_plas": "DE_results_pVIM-1-vs-PF_filtered_padj_plasmids.tsv",
        "pVIM-1-vs-pOXA-48_plas": "DE_results_pVIM-1-vs-pOXA-48_filtered_padj_plasmids.tsv",
        "pVIM-1-vs-pOXA-48_pOXA-VIM": "DE_results_pVIM-1-vs-pOXA-48_filtered_padj.tsv"
    }

    dict_comploc_strain_df = {}
    for comploc, file in files.items():
        dict_strain_df = {}
        for directory in dir_names:
            if os.path.isdir(directory):
                dict_strain_df.update(table_to_dict(directory, file))
            else:
                print(f"{directory} is an invalid directory name")
        dict_comploc_strain_df[comploc] = dict_strain_df

    ### Preparing dataframes:

    dict_comploc_strain_df_sub = {}
    for comploc in files.keys():
        if comploc.endswith("chr") or comploc.endswith("pOXA-VIM"):
            dict_comploc_strain_df_sub[comploc] = prepare_dfs_in_dict_for_merging(dict_comploc_strain_df[comploc])

    ### Merging dataframes:

    dict_comploc_merged_dfs = {}
    for comploc, strain_df_sub in dict_comploc_strain_df_sub.items():
        dict_comploc_merged_dfs[comploc] = merge_dfs_in_dict_from_a_comparative_analysis(strain_df_sub)

    ### Exporting to excel:

    with pd.ExcelWriter(f"{out_name}.xlsx", engine="xlsxwriter") as writer:
        highlighted_rows_all = []
        for comploc in files.keys():
            try:
                if comploc.endswith("chr") or comploc.endswith("pOXA-VIM"):
                    sheet_name = f"summ_{comploc}"
                    merged_df_result = dict_comploc_merged_dfs[comploc]
                    merged_df_result.to_excel(writer, sheet_name=sheet_name, index=False)

                    workbook = writer.book
                    bold_format = workbook.add_format({"bold": True})

                    if comploc.endswith("chr"):
                        highlighted_rows_all.extend(apply_bold_format_to_iron_related_degs_and_store(merged_df_result, sheet_name, writer))

            except Exception as e:
                print(f"An unexpected error in comploc {comploc} occurred: {e}")

        if highlighted_rows_all:
            summary_df = pd.DataFrame(highlighted_rows_all)
            col = summary_df.pop("SourceSheet")
            summary_df.insert(0, col.name, col)
            summary_df.to_excel(writer, sheet_name="iron-related_DEGs", index=False)

        for comploc in files.keys():
            try:
                if comploc.endswith("chr") or comploc.endswith("plas"):
                    dict_df_result = dict_comploc_strain_df[comploc]
                    dict_to_excel(dict_df_result, comploc, writer)
            except Exception as e:
                print(f"An unexpected error in com_loc {comploc} occurred: {e}")

    print(f"Excel {out_name}.xlsx successfully created.")
