# -*- coding: utf-8 -*-
"""
This file reconstructs scoGEM, the genome-scale model for Streptomyces coelicolor A3(2).
Author: Snorre Sulheim
Created: 28.02.2019
Email: snorre.sulheim@sintef.no


# Description
The main function of this script is to update the curated annotation of subsystem and pathways.
However, all work related to subsystem / pathway annotation goes in here, e.g. creating 
csv-file used to export, querying KEGG or BioCyc for annotations etc
"""

import cobra
import pandas as pd 
from pathlib import Path

REPO_DIR = Path(__file__).parent.parent.parent

def export_reaction_subsystem_and_pathways(model, csv_fn):
    """
    Use pandas to write a csv-file which can be used to curate the subsystem
    annotation. The csv-file is written with the following columns:
    
    Reaction ID, Reaction name, KEGG ID, Biocyc-annotation, 
    KEGG Subsystem, KEGG pathway, Subsystem
    """
    annotation_keys = ["kegg.reaction", "biocyc", "kegg.pathway", "kegg.subsystem", "subsystem"]
    reactions_list = []
    for r in model.reactions:
        r_list = [r.id, r.name]

        for key in annotation_keys:
            try:
                r_list.append(r.annotation[key])
            except KeyError:
                r_list.append(None)
        reactions_list.append(r_list)
    df = pd.DataFrame(reactions_list)
    df.columns = ["Reaction ID", "Reaction name"] +  annotation_keys
    # Add empty column
    df["curated pathway"] = None 
    df["curated subsystem"] = None 
    print(df.head())
    df.to_csv(csv_fn, sep = ";")
    



if __name__ == '__main__':
    model_fn = REPO_DIR / "ModelFiles" / "xml" / "scoGEM.xml"
    model = cobra.io.read_sbml_model(str(model_fn))
    csv_fn = str(REPO_DIR / "ComplementaryData" / "curation" / "subsystem.csv")
    export_reaction_subsystem_and_pathways(model, csv_fn)