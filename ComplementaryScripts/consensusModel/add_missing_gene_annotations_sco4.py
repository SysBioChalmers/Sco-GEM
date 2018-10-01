#!/usr/bin/env python
"""
Author: Snorre Sulheim
Created: 19.09.2018


# Description
This file adds the missing gene annotations identified by Raven 2.0 in the Sco4 reconstruction.
This list of file is given in Tabls S7 in the supplementaryTables following the Sco4 publication

Rxn ID  NOTE                                                        Gene assocation
RMI     Required for growth using L-rhamnose as a carbon source     SCO0812
RMPA    Required for growth using L-rhamnose as a carbon source     SCO0813
DABTD   Required for growth using D-arabitol as a carbon source     SCO1901
PROD2   Required for growth using L-proline as a carbon source      SCO5519
THRPD   Required for growth (adenosylcobalamin biosynthesis)        SCO1859
ADCL    Required for growth (tetrahydrofolate biosynthesis)         SCO1546
CBLAT   To utilize extracellular cob(1)alamin                       SCO1851 or SCO5381
BPNT    Required for growth (acyl-carrier protein biosynthesis)     SCO5161
DTMPK   Required for growth (dTTP biosynthesis)                     SCO3542
OXPTNDH Required for growth using L-lysine as a carbon source       SCO1204 or SCO1612 or SCO1706 or SCO3420 or SCO3486 or SCO4780 or SCO5666 or SCO5679 or SCO6441 or SCO6793 or SCO7035 or SCO7139
GLNTRS  tRNA system                                                 SCO5547

# Note
The follwoing reactions has changed id in iKS1317
- DABTD -> ABTDG
- THRPD -> THRPDC


# From compare_iMK1208_content.csv comparison
We find additional reactions in Sco4 where the gene annotation is changed.
Most of the changes in gene annotation in compare_iMK1208_content.csv is that they have removed
s00001

PDH( SCO1269 and SCO1270 and SCO2183 and SCO0884 and SCO1268 ) or ( SCO1269 and SCO1270 and SCO2183 and SCO0884 and SCO7123 ) or ( SCO1269 and SCO1270 and SCO2183 and SCO0884 and SCO2181 ) or ( SCO1269 and SCO1270 and SCO2183 and SCO2180 and SCO1268 ) or ( SCO1269 and SCO1270 and SCO2183 and SCO2180 and SCO7123 ) or ( SCO1269 and SCO1270 and SCO2183 and SCO2180 and SCO2181 ) or ( SCO1269 and SCO1270 and SCO2183 and SCO4919 and SCO1268 ) or ( SCO1269 and SCO1270 and SCO2183 and SCO4919 and SCO7123 ) or ( SCO1269 and SCO1270 and SCO2183 and SCO4919 and SCO2181 ) or ( SCO1269 and SCO1270 and SCO2371 and SCO0884 and SCO1268 ) or ( SCO1269 and SCO1270 and SCO2371 and SCO0884 and SCO7123 ) or ( SCO1269 and SCO1270 and SCO2371 and SCO0884 and SCO2181 ) or ( SCO1269 and SCO1270 and SCO2371 and SCO2180 and SCO1268 ) or ( SCO1269 and SCO1270 and SCO2371 and SCO2180 and SCO7123 ) or ( SCO1269 and SCO1270 and SCO2371 and SCO2180 and SCO2181 ) or ( SCO1269 and SCO1270 and SCO2371 and SCO4919 and SCO1268 ) or ( SCO1269 and SCO1270 and SCO2371 and SCO4919 and SCO7123 ) or ( SCO1269 and SCO1270 and SCO2371 and SCO4919 and SCO2181 ) or ( SCO1269 and SCO1270 and SCO7124 and SCO0884 and SCO1268 ) or ( SCO1269 and SCO1270 and SCO7124 and SCO0884 and SCO7123 ) or ( SCO1269 and SCO1270 and SCO7124 and SCO0884 and SCO2181 ) or ( SCO1269 and SCO1270 and SCO7124 and SCO2180 and SCO1268 ) or ( SCO1269 and SCO1270 and SCO7124 and SCO2180 and SCO7123 ) or ( SCO1269 and SCO1270 and SCO7124 and SCO2180 and SCO2181 ) or ( SCO1269 and SCO1270 and SCO7124 and SCO4919 and SCO1268 ) or ( SCO1269 and SCO1270 and SCO7124 and SCO4919 and SCO7123 ) or ( SCO1269 and SCO1270 and SCO7124 and SCO4919 and SCO2181 )
CU2abc; SCO0164 or SCO6460
CBIabc: SCO2273 #Changed from CBIuabcto CBIabc
CBL1abc: SCO2273
GSNt2: SCO3915 or SCO0079
INSt2: SCO3915 or SCO0079 # Changed from INSt2r
GLNTRS: SCO5547

The final list of reactions where the gene annotation is changed is REACTION_GENE_DICT.
The uncommented reactions are already in iKS1317

"""
import logging
import pandas as pd
import cobra

COMPARE_iMK1208_CONTENT_FN = "../../ComplementaryData/curation/compare_iMK1208_content.csv"

REACTION_GENE_DICT = {
    # "RMI":     "SCO0812", # change is already in iKS1317
    "RMPA":    "SCO0813",
    "ABTDG":   "SCO1901",
    "PROD2":   "SCO5519",
    "THRPDC":   "SCO1859",
    "ADCL":    "SCO1546",
    # "CBLAT":   "SCO1851 or SCO5381", # change is already in iKS1317
    # "BPNT":    "SCO5161", 
    # "DTMPK":   "SCO3542", # change is already in iKS1317
    "OXPTNDH": "SCO1204 or SCO1612 or SCO1706 or SCO3420 or SCO3486 or SCO4780 or SCO5666 or SCO5679 or SCO6441 or SCO6793 or SCO7035 or SCO7139",
    "GLNTRS":  "SCO5547",
    "CU2abc":  "SCO0164 or SCO6460",
    "CBIabc": "SCO2273",
    "CBL1abc": "SCO2273",
    "GSNt2": "SCO3915 or SCO0079",
    "INSt2": "SCO3915 or SCO0079",
    "PDH":  "SCO1269 and SCO1270 and (SCO2183 or SCO2371 or SCO7124) and (SCO0884 or SCO04919 or SCO2180) and (SCO1268 or SCO7123 or SCO2181) and (SCO3815 or SCO3829)" #This is a union of Sco4 and iKS1317 annotations

}





def add_gene_annotations(scoGEM):
    for reaction_id, gene_string in REACTION_GENE_DICT.items():
        try:
            reaction = scoGEM.reactions.get_by_id(reaction_id)
        except KeyError:
            logging.warning("The gene annotation {0} could not be added to reaction {1}".format(gene_string, reaction))
        else:
            reaction.gene_reaction_rule = gene_string

def get_changed_genes_in_sco4():
    df = pd.read_csv(COMPARE_iMK1208_CONTENT_FN, sep = ";")
    df_changed_sco4 = df[df["Sco4"].str.contains("genes")]
    return df_changed_sco4

def print_new_sco4_annotations(iMK1208, iKS1317, sco4):
    df_changed_sco4 = get_changed_genes_in_sco4()
    print(df_changed_sco4)
    for i, row in df_changed_sco4.iterrows():
        reaction_id = row["iMK1208 ID"]
        r_sco4 = sco4.reactions.get_by_id(reaction_id)
        r_iMK1208 = iMK1208.reactions.get_by_id(reaction_id)
        try:
            r = iKS1317.reactions.get_by_id(reaction_id)
        except:
            print("missing", reaction_id)
            print("{0}; {1}; {2}; {3}; {4}".format(reaction_id, r_iMK1208.gene_reaction_rule, "None", r_sco4.gene_reaction_rule, True))
        else:
            print("{0}; {1}; {2}; {3}; {4}".format(reaction_id, r_iMK1208.gene_reaction_rule, r.gene_reaction_rule, r_sco4.gene_reaction_rule, True))

if __name__ == '__main__':
    sco4_fn = r"C:\\Users\\snorres\\git\\scoGEM\\ComplementaryData\\models\\Sco4.xml"
    iKS1317_fn = r"C:\\Users\\snorres\\git\\scoGEM\\ComplementaryData\\models\\iKS1317.xml"
    iMK1208_fn = r"C:\\Users\\snorres\\OneDrive - SINTEF\\SINTEF projects\\INBioPharm\\SCM\\sbml models\\Kim2014\\kim_with_kegg.xml"

    iMK1208 = cobra.io.read_sbml_model(iMK1208_fn)
    sco4 = cobra.io.read_sbml_model(sco4_fn)
    iKS1317 = cobra.io.read_sbml_model(iKS1317_fn)
    # print_new_sco4_annotations(iMK1208, iKS1317, sco4)
    print([x.id for x in sco4.reactions.PDH.genes])
    print([x.id for x in iKS1317.reactions.PDH.genes])
    print(iKS1317.reactions.PDH.gene_reaction_rule)
