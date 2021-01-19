# -*- coding: utf-8 -*-
import pandas as pd
import logging

"""
Author: Snorre Sulheim
Created: 10.01.2018

# Description
This file adds doi annotations from ComplementaryData/annotations/reaction_notes_and_references.csv
"""

def add_doi_annotations(scoGEM, annotations_fn):
    df = pd.read_csv(annotations_fn, sep = ";", encoding = "ISO-8859-1")
    for index, row in df.iterrows():
        if pd.notna(row["DOI"]):
            r_id = row["Reaction ID"]
            doi_list = [x.strip() for x in row["DOI"].split(",")]
            r = scoGEM.reactions.get_by_id(r_id)
            r.annotation["doi"] = doi_list
            logging.info("Added doi annotation {1} to {0}".format(r_id, doi_list))

def add_gene_annotations(scoGEM, annotations_fn):
    """
    Add the following identifiers to all genes (if available):
    1) Uniprot ID
    2) GO id
    3) Pfam
    4) RefSeq

    Didn't find the correct MIRIAM identifier for PANTHER
    """
    df = pd.read_csv(annotations_fn, sep = ";", skiprows = [1])
    df.set_index("SCO_ID", inplace = True)
    for index, row in df.iterrows():
        gene_id = index.replace(".", "")
        try:
            gene = scoGEM.genes.get_by_id(gene_id)
        except KeyError:
            print("Gene {0} is not in the model".format(gene_id))
            continue
        # Uniprot
        gene.annotation["uniprot"] = row["Uniprot_Entry"]
        annotations = [row["Uniprot_Entry"]]

        # GO
        if pd.notna(row["GO_ID"]):
            GO_list = [x.strip() for x in row["GO_ID"].split(";") if len(x)]
            gene.annotation["go"] = GO_list
            annotations += GO_list

        # Pfam
        if pd.notna(row["Pfam"]):
            pfam_list = [x.strip() for x in row["Pfam"].split(";") if len(x)]
            gene.annotation["pfam"] = pfam_list
            annotations += pfam_list


        # RefSeq
        if pd.notna(row["RefSeq"]):
            refseq_list = [x.strip() for x in row["RefSeq"].split(";") if len(x)]
            gene.annotation["refseq"] = refseq_list
            annotations += refseq_list

        logging.info("{0} was given the following annotations: {1}".format(gene_id, ", ".join(annotations)))


if __name__ == '__main__':
    import cobra
    scoGEM_FN = "../../ModelFiles/xml/scoGEM.xml"
    scoGEM = cobra.io.read_sbml_model(scoGEM_FN)
    # annotations_fn = "../../ComplementaryData/annotations/reaction_notes_and_references.csv"
    # add_doi_annotations(scoGEM, annotations_fn)
    annotations_fn = "../../ComplementaryData/annotations/genes.csv"
    add_gene_annotations(scoGEM, annotations_fn)