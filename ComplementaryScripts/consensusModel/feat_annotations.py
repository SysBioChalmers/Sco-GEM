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

if __name__ == '__main__':
    import cobra
    scoGEM_FN = "../../ModelFiles/xml/scoGEM.xml"
    scoGEM = cobra.io.read_sbml_model(scoGEM_FN)
    annotations_fn = "../../ComplementaryData/annotations/reaction_notes_and_references.csv"
    add_doi_annotations(scoGEM, annotations_fn)